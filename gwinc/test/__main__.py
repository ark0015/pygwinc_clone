"""Beginnings of a test suite for pygwinc

"""
import os
import sys
import signal
import pickle
import subprocess
import hashlib
import numpy as np
import matplotlib.pyplot as plt
import argparse

import logging
logging.basicConfig(format='%(message)s',
                    level=logging.DEBUG)

from .. import load_ifo, gwinc


FLO = 5
FHI = 6000
NPOINTS = 3000
FRACTIONAL_TOLERANCE = 0.01
# comparisons to skip
SKIP = [
    # 'Seismic',
    # 'Suspension Thermal',
    # 'Newtonian Gravity',
    'Total',
    ]


def path_hash(path):
    """Calculate SHA1 hash of path, either directory or file"""
    path = os.path.expanduser(path)
    if os.path.isdir(path):
        d = path
        f = '.'
    else:
        d = os.path.dirname(path)
        f = os.path.basename(path)
    CWD = os.getcwd()
    os.chdir(d)
    cmd = 'find {} -type f ! -wholename "*/.*" -print0 | sort -z | xargs -0 sha1sum | sha1sum'.format(f)
    sha1sum_out = subprocess.check_output(cmd, shell=True)
    sha1sum = sha1sum_out.split()[0]
    os.chdir(CWD)
    return sha1sum



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plot', '-p', action='store_true', help='plot differences')
    parser.add_argument('--save', '-s', help='save plot to file')
    parser.add_argument('IFO', help='IFO name or description file')
    args = parser.parse_args()

    logging.info("loading IFO '{}'...".format(args.IFO))
    ifo = load_ifo(args.IFO)

    freq = np.logspace(np.log10(FLO), np.log10(FHI), NPOINTS)

    mdata_pkl = os.path.join(os.path.dirname(__file__), '{}.pkl'.format(args.IFO))
    ifo_hash = hashlib.sha1(ifo.to_txt()).hexdigest()
    gwinc_hash = path_hash(os.getenv('GWINCPATH'))

    mrecalc = True

    if os.path.exists(mdata_pkl):
        logging.info("loading MATLAB data {}...".format(mdata_pkl))
        with open(mdata_pkl, 'rb') as f:
            if sys.version_info.major > 2:
                mdata = pickle.load(f, encoding='latin1')
            else:
                mdata = pickle.load(f)

        # don't recalculcate MAT gwinc if ifo and gwinc dir hashes
        # haven't changed.
        if mdata['ifo_hash'] == ifo_hash and mdata['gwinc_hash'] == gwinc_hash:
            mrecalc = False
        if mdata['ifo_hash'] != ifo_hash:
            logging.info("ifo hash has changed: {}".format(ifo_hash))
        if mdata['gwinc_hash'] != gwinc_hash:
            logging.info("GWINC hash has changed: {}".format(gwinc_hash))

    if mrecalc:
        logging.info("calculating MATLAB noises...")
        from ..gwinc_matlab import gwinc_matlab
        mscore, mnoises, mifo = gwinc_matlab(freq, ifo)
        mdata = dict(score=mscore, noises=mnoises, ifo=mifo, ifo_hash=ifo_hash, gwinc_hash=gwinc_hash)
        with open(mdata_pkl, 'wb') as f:
            pickle.dump(mdata, f)

    score, noises, ifo = gwinc(freq, ifo)

    mnoises = mdata['noises']

    diffs = {}
    for name, noise in noises.items():
        if name == 'Freq':
            continue
        if name in SKIP:
            logging.warning("SKIPPING TEST: '{}'".format(name))
            continue

        try:
            mnoise = mnoises[name]
        except KeyError:
            continue
        # logging.info("compare: {}".format(name))

        diff = np.sqrt(mnoise) - np.sqrt(noise)
        frac = abs(diff / np.sqrt(noise))

        if max(frac) < FRACTIONAL_TOLERANCE:
            continue

        logging.warning("EXCESSIVE DIFFERENCE: '{}'".format(name))
        logging.warning("  max: {:e}, min: {:e}".format(max(frac), min(frac)))

        diffs[name] = frac


    if args.plot:
        spec = (len(diffs), 2)
        for i, name in enumerate(diffs):
            axl = plt.subplot2grid(spec, (i, 0))
            axl.loglog(freq, np.sqrt(noises[name]), label='pygwinc')
            axl.loglog(freq, np.sqrt(mnoises[name]), label='matlab')
            axl.grid()
            axl.legend(loc='upper right')
            # ax.set_title(name)
            axl.set_ylabel(name)
    
            axr = plt.subplot2grid(spec, (i, 1))
            axr.loglog(freq, diffs[name], label=name)
            axr.grid()
            # ax.set_title(name)

        axl.set_xlabel("frequency [Hz]")
        axr.set_xlabel("frequency [Hz]")
    
        plt.suptitle("noises that differ by more than 1% [(mat-py)/py]")
        if args.save:
            plt.gcf().set_size_inches(11, 20)
            plt.savefig(args.save)
        else:
            plt.show()

    if len(diffs) > 0:
        return 1
    return 0


if __name__ == '__main__':
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    sys.exit(main())
