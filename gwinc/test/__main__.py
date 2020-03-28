import os
import sys
import glob
import signal
import logging
import tempfile
import argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from collections.abc import Mapping
from PyPDF2 import PdfFileReader, PdfFileWriter

logging.basicConfig(format='%(message)s',
                    level=os.getenv('LOG_LEVEL', logging.INFO))

from .. import IFOS, load_budget
from ..io import load_hdf5, save_hdf5

try:
    import inspiral_range
except ImportError:
    inspiral_range = None


TOLERANCE = 1e-6
CACHE_PATH = os.path.join(os.path.dirname(__file__), 'cache')


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


def plot_diffs(freq, diffs, tolerance,
               name, labelA, labelB, fom_title='',
               save=None):
    spec = (len(diffs)+1, 2)
    sharex = None
    for i, nname in enumerate(diffs):
        noiseA, noiseB, frac = diffs[nname]

        axl = plt.subplot2grid(spec, (i, 0), sharex=None)
        axl.loglog(freq, np.sqrt(noiseA), label=labelA)
        axl.loglog(freq, np.sqrt(noiseB), label=labelB)
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

        plt.suptitle('''{} {}/{} noise comparison
(noises that differ by more than {} ppm)
{}'''.format(name, labelA, labelB, tolerance*1e6, fom_title))

    axl.set_xlabel("frequency [Hz]")
    axr.set_xlabel("frequency [Hz]")

    plt.subplots_adjust(top=0.8, right=0.85, wspace=0.3)
    if save:
        pwidth = 10
        pheight = (len(diffs) * 5) + 2
        plt.gcf().set_size_inches(pwidth, pheight)
        plt.savefig(save)
    else:
        plt.show()

##################################################

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tolerance', '-t',  type=float, default=TOLERANCE,
                        help='fractional tolerance [{}]'.format(TOLERANCE))
    parser.add_argument('--skip', '-k', metavar='NOISE', action='append',
                        help='traces to skip in comparison (multiple may be specified)')
    parser.add_argument('--cache', '-c', metavar='PATH', default=CACHE_PATH,
                        help='specify alternate IFO traces cache path')
    rgroup = parser.add_mutually_exclusive_group()
    rgroup.add_argument('--plot', '-p', action='store_true',
                        help='plot differences')
    rgroup.add_argument('--report', '-r', metavar='REPORT.pdf',
                        help='create PDF report of test results (only created if differences found)')
    rgroup.add_argument('--gen-cache', action='store_true',
                        help='update/create IFO traces cache directory')
    parser.add_argument('ifo', metavar='IFO', nargs='*',
                        help='specific ifos to test (default all)')
    args = parser.parse_args()

    if args.gen_cache:
        try:
            os.makedirs(args.cache)
        except FileExistsError:
            pass
        freq = np.logspace(np.log10(5), np.log10(6000), 3000)
        for name in IFOS:
            Budget = load_budget(name)
            traces = Budget(freq).run()
            path = os.path.join(args.cache, name+'.h5')
            save_hdf5(path, freq, traces)
        return
        
    if args.report:
        base, ext = os.path.splitext(args.report)
        if ext != '.pdf':
            parser.error("Test reports only support PDF format.")
        outdir = tempfile.TemporaryDirectory()
    
    # find all cached IFOs
    logging.info("loading cache {}...".format(args.cache))
    cached_ifos = {}
    for f in sorted(os.listdir(args.cache)):
        name, ext = os.path.splitext(f)
        if ext != '.h5':
            continue
        cached_ifos[name] = os.path.join(args.cache, f)

    # select
    if args.ifo:
        ifos = {name:cached_ifos[name] for name in args.ifo}
    else:
        ifos = cached_ifos

    labelA = 'cache'
    labelB = 'head'

    fail = False

    # compare
    for name, path in ifos.items():
        logging.info("{} tests...".format(name))

        freq, tracesA, attrs = load_hdf5(path)

        Budget = load_budget(name)
        tracesB = Budget(freq).run()

        if inspiral_range:
            totalA = tracesA['Total'][0]
            totalB = tracesB['Total'][0]
            range_func = inspiral_range.range
            H = inspiral_range.waveform.CBCWaveform(freq)
            fomA = range_func(freq, totalA, H=H)
            tracesA['int73'] = inspiral_range.int73(freq, totalA)[1], None
            fomB = range_func(freq, totalB, H=H)
            tracesB['int73'] = inspiral_range.int73(freq, totalB)[1], None
            fom_summary = """
inspiral {func} {m1}/{m2} Msol:
{labelA}: {fomA:.2f} Mpc
{labelB}: {fomB:.2f} Mpc
""".format(
                func=range_func.__name__,
                m1=H.params['m1'],
                m2=H.params['m2'],
                labelA=labelA,
                fomA=fomA,
                labelB=labelB,
                fomB=fomB,
            )
        else:
            fom_summary = ''

        diffs = compare_traces(tracesA, tracesB, args.tolerance, args.skip)

        if diffs:
            logging.warning("{} tests FAIL".format(name))
            fail |= True
            if args.plot or args.report:
                if args.report:
                    save = os.path.join(outdir.name, name+'.pdf')
                else:
                    save = None
                plot_diffs(
                    freq, diffs, args.tolerance,
                    name, labelA, labelB, fom_summary,
                    save=save,
                )
        else:
            logging.info("{} tests pass.".format(name))

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
