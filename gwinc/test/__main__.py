import os
import sys
import glob
import shutil
import signal
import logging
import tempfile
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from collections.abc import Mapping
from PyPDF2 import PdfFileReader, PdfFileWriter

from .. import IFOS, load_budget
from ..io import load_hdf5

try:
    import inspiral_range
except ImportError:
    inspiral_range = None


logging.basicConfig(
    format='%(message)s',
    level=os.getenv('LOG_LEVEL', logging.INFO))


TOLERANCE = 1e-6
CACHE_LIMIT = 5


def test_path(*args):
    """Return path to package file."""
    return os.path.join(os.path.dirname(__file__), *args)


def git_ref_resolve_hash(git_ref):
    """Resolve a git reference into its hash string."""
    try:
        return subprocess.run(
            ['git', 'show', '-s', '--format=format:%H', git_ref],
            capture_output=True, universal_newlines=True,
        ).stdout
    except subprocess.CalledProcessError:
        return None


def write_ref_hash(ref_hash):
    """Write ref hash to reference file

    """
    with open(test_path('ref_hash'), 'w') as f:
        f.write('{}\n'.format(ref_hash))


def load_ref_hash():
    """Load the current reference git hash.

    """
    try:
        with open(test_path('ref_hash')) as f:
            return f.read().strip()
    except IOError:
        return None


def prune_cache_dir():
    """Prune all but the N most recently accessed caches.

    """
    cache_dir = test_path('cache')
    if not os.path.exists(cache_dir):
        return
    expired_paths = sorted(
        [os.path.join(cache_dir, path) for path in os.listdir(cache_dir)],
        key=lambda path: os.stat(path).st_atime,
    )[CACHE_LIMIT:]
    if not expired_paths:
        return
    logging.info("pruning {} old caches...".format(len(expired_paths)))
    for path in expired_paths:
        logging.debug("pruning {}...".format(path))
        shutil.rmtree(path)


def gen_cache_for_ref(ref_hash, path):
    """generate cache from git reference

    The ref_hash should be a git hash, and path should be the location
    of the generated cache.

    The included shell script is used to extract the gwinc code from
    the appropriate git commit, and invoke a new python instance to
    generate the noise curves.

    """
    logging.info("creating new cache from reference {}...".format(ref_hash))
    subprocess.run(
        [test_path('gen_cache.sh'), ref_hash, path],
        check=True,
    )


def load_cache(path):
    """load a cache from path

    returns a dictionary with 'ref_hash' and 'ifos' keys.

    """
    logging.info("loading cache {}...".format(path))
    cache = {}
    ref_hash_path = os.path.join(path, 'ref_hash')
    if os.path.exists(ref_hash_path):
        with open(ref_hash_path) as f:
            ref_hash = f.read().strip()
    else:
        ref_hash = None
    logging.debug("cache hash: {}".format(ref_hash))
    cache['ref_hash'] = ref_hash
    cache['ifos'] = {}
    for f in sorted(os.listdir(path)):
        name, ext = os.path.splitext(f)
        if ext != '.h5':
            continue
        cache['ifos'][name] = os.path.join(path, f)
    return cache


def walk_traces(traces, root=()):
    """recursively walk through traces

    yields (key_tuple, value), where key_tuple is a tuple of nested
    dict keys of the traces dict locating the given value.

    """
    for key, val in traces.items():
        r = root+(key,)
        if isinstance(val, Mapping):
            for k, v in walk_traces(val, root=r):
                yield k, v
            continue
        else:
            yield r, val


def zip_noises(tracesA, tracesB, skip):
    """zip matching noises from traces"""
    for keys, (noiseA, _) in walk_traces(tracesA):
        # keys is a tuple of nested keys for noise
        name = keys[-1]
        # skip keys we'll add at the end
        if name in ['Total', 'int73']:
            continue
        if skip and name in skip:
            logging.warning("SKIPPING TEST: '{}'".format(name))
            continue

        try:
            t = tracesB
            for key in keys:
                t = t[key]
            noiseB, style = t
        except (KeyError, TypeError):
            logging.warning("MISSING B TRACE: '{}'".format(name))
            continue

        yield name, noiseA, noiseB

    yield 'Total', tracesA['Total'][0], tracesB['Total'][0]

    if 'int73' in tracesA:
        yield 'int73', tracesA['int73'][0], tracesB['int73'][0]


def compare_traces(tracesA, tracesB, tolerance=TOLERANCE, skip=None):
    """Compare two gwinc trace dictionaries

    Noises listed in `skip` will not be compared.

    Returns a dictionary of noises that differ fractionally (on a
    point-by-point basis) by more that `tolerance` between `tracesA`
    and `tracesB`.

    """
    #name_width = max([len(n[0][-1]) for n in walk_traces(tracesA)])
    name_width = 15
    diffs = OrderedDict()
    for name, noiseA, noiseB in zip_noises(tracesA, tracesB, skip):
        logging.debug("comparing {}...".format(name))

        ampA = np.sqrt(noiseA)
        ampB = np.sqrt(noiseB)
        diff = ampB - ampA
        frac = abs(diff / ampA)

        if max(frac) < tolerance:
            continue

        logging.warning("EXCESSIVE DIFFERENCE: {:{w}} {:6.1f} ppm".format(
            name, max(frac)*1e6, w=name_width))

        diffs[name] = (noiseA, noiseB, frac)

    return diffs


def plot_diffs(freq, diffs, styleA, styleB):
    spec = (len(diffs)+1, 2)
    sharex = None
    for i, nname in enumerate(diffs):
        noiseA, noiseB, frac = diffs[nname]

        axl = plt.subplot2grid(spec, (i, 0), sharex=None)
        axl.loglog(freq, np.sqrt(noiseA), **styleA)
        axl.loglog(freq, np.sqrt(noiseB), **styleB)
        axl.grid()
        axl.legend(loc='upper right')
        axl.set_ylabel(nname)
        if i == 0:
            axl.set_title("noise value")

        if i == 0:
            sharex = axl
        axr = plt.subplot2grid(spec, (i, 1), sharex=sharex)
        axr.loglog(freq, frac)
        axr.grid()
        axr.axhline(y=max(frac), color='r', linestyle='--')
        axr.text(max(freq)+4000, max(frac), '{:.1f} ppm'.format(max(frac)*1e6),
                 horizontalalignment='left', verticalalignment='center',
                 color='red')
        if i == 0:
            axr.set_title("fractional difference")

    axl.set_xlabel("frequency [Hz]")
    axr.set_xlabel("frequency [Hz]")
    plt.subplots_adjust(top=0.8, right=0.85, wspace=0.3)

##################################################

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--tolerance', '-t',  type=float, default=TOLERANCE,
        help='fractional tolerance [{}]'.format(TOLERANCE))
    parser.add_argument(
        '--skip', '-k', metavar='NOISE', action='append',
        help='traces to skip in comparison (multiple may be specified)')
    parser.add_argument(
        '--git-ref', '-g', metavar='HASH',
        help='specify git ref to compare against')
    rgroup = parser.add_mutually_exclusive_group()
    rgroup.add_argument(
        '--update-ref', '-u', metavar='HASH', nargs='?', const='HEAD',
        help="update the stored reference git hash to HASH (or 'HEAD' if not specified) and exit")
    rgroup.add_argument(
        '--plot', '-p', action='store_true',
        help='plot differences')
    rgroup.add_argument(
        '--report', '-r', metavar='REPORT.pdf',
        help='create PDF report of test results (only created if differences found)')
    parser.add_argument(
        'ifo', metavar='IFO', nargs='*',
        help='specific ifos to test (default all)')
    args = parser.parse_args()

    # get the current hash of HEAD
    head_hash = git_ref_resolve_hash('HEAD')
    if not head_hash:
        logging.warning("could not determine git HEAD hash.")

    # update the reference if specified
    if args.update_ref:
        if args.update_ref == 'HEAD':
            if not head_hash:
                sys.exit("Could not update reference to head.")
            logging.info("updating reference to HEAD...")
            ref_hash = head_hash
        else:
            ref_hash = git_ref_resolve_hash(args.update_ref)
        logging.info("updating reference git hash: {}".format(ref_hash))
        write_ref_hash(ref_hash)
        sys.exit()

    # get the reference hash
    if args.git_ref:
        ref_hash = git_ref_resolve_hash(args.git_ref)
    else:
        ref_hash = load_ref_hash()
        if not ref_hash:
            pass
        try:
            with open(test_path('ref_hash')) as f:
                ref_hash = f.read().strip()
        except IOError:
            logging.warning("could not open reference")
            sys.exit("Unspecified reference git hash, could not run test.")

    logging.info("head hash: {}".format(head_hash))
    logging.info("ref  hash: {}".format(ref_hash))

    # don't bother test if hashes match
    if ref_hash == head_hash:
        logging.info("HEAD matches reference, not bothering to calculate.")
        logging.info("Use --git-ref to compare against an arbitrary git commit.")
        sys.exit()

    # load the cache
    cache_path = test_path('cache', ref_hash)
    if not os.path.exists(cache_path):
        prune_cache_dir()
        gen_cache_for_ref(ref_hash, cache_path)
    cache = load_cache(cache_path)

    if args.report:
        base, ext = os.path.splitext(args.report)
        if ext != '.pdf':
            parser.error("Test reports only support PDF format.")
        outdir = tempfile.TemporaryDirectory()

    if args.ifo:
        ifos = args.ifo
    else:
        ifos = IFOS

    style_cache = dict(label='reference', linestyle='-')
    style_head = dict(label='head', linestyle='--')

    fail = False

    # compare
    for name in ifos:
        logging.info("{} tests...".format(name))

        try:
            path = cache['ifos'][name]
        except KeyError:
            logging.warning("IFO {} not found in cache")
            fail |= True
            continue

        freq, traces_cache, attrs = load_hdf5(path)

        Budget = load_budget(name)
        traces_head = Budget(freq).run()

        if inspiral_range:
            total_cache = traces_cache['Total'][0]
            total_head = traces_head['Total'][0]
            range_func = inspiral_range.range
            H = inspiral_range.waveform.CBCWaveform(freq)
            fom_cache = range_func(freq, total_cache, H=H)
            traces_cache['int73'] = inspiral_range.int73(freq, total_cache)[1], None
            fom_head = range_func(freq, total_head, H=H)
            traces_head['int73'] = inspiral_range.int73(freq, total_head)[1], None
            fom_summary = """
inspiral {func} {m1}/{m2} Msol:
{label_cache}: {fom_cache:.2f} Mpc
{label_head}: {fom_head:.2f} Mpc
""".format(
                func=range_func.__name__,
                m1=H.params['m1'],
                m2=H.params['m2'],
                label_cache=style_cache['label'],
                fom_cache=fom_cache,
                label_head=style_head['label'],
                fom_head=fom_head,
            )
        else:
            fom_summary = ''

        diffs = compare_traces(traces_cache, traces_head, args.tolerance, args.skip)

        if not diffs:
            logging.info("{} tests pass.".format(name))
            continue

        logging.warning("{} tests FAIL".format(name))
        fail |= True
        if args.plot or args.report:
            plot_diffs(freq, diffs, style_cache, style_head)
            plt.suptitle('''{} {}/{} noise comparison
(noises that differ by more than {} ppm)
reference git hash: {}
{}'''.format(name, style_cache['label'], style_head['label'],
             args.tolerance*1e6, cache['ref_hash'], fom_summary))
            if args.report:
                pwidth = 10
                pheight = (len(diffs) * 5) + 2
                plt.gcf().set_size_inches(pwidth, pheight)
                plt.savefig(os.path.join(outdir.name, name+'.pdf'))
            else:
                plt.show()

    if not fail:
        logging.info("all tests pass.")
        return 0

    if args.report:
        logging.info("generating report {}...".format(args.report))
        pdf_writer = PdfFileWriter()
        for path in glob.glob(os.path.join(outdir.name, '*.pdf')):
            pdf_reader = PdfFileReader(path)
            for page in range(pdf_reader.getNumPages()):
                pdf_writer.addPage(pdf_reader.getPage(page))
        with open(args.report, 'wb') as f:
            pdf_writer.write(f)
        outdir.cleanup()

    logging.info("TESTS FAILED.")
    return 1

##################################################

if __name__ == '__main__':
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    sys.exit(main())
