"""Beginnings of a test suite for pygwinc

"""
import os
import sys
import signal
import pickle
import numpy as np
import matplotlib.pyplot as plt
import argparse

import logging
logging.basicConfig(level=logging.WARNING)

from .. import load_ifo, gwinc


IFO = 'aLIGO'
FLO = 5
FHI = 6000
NPOINTS = 3000
# comparisons to skip
SKIP = [
    'Seismic',
    'Suspension Thermal',
    'Newtonian Gravity',
    'Total',
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plot', '-p', action='store_true', help='plot differences')
    args = parser.parse_args()

    ifo = load_ifo(IFO)

    freq = np.logspace(np.log10(FLO), np.log10(FHI), NPOINTS)

    matdata = os.path.join(os.path.dirname(__file__), 'matlab.pkl')
    if os.path.exists(matdata):
        logging.info("using existing {}...".format(matdata))
        with open(matdata, 'rb') as f:
            if sys.version_info.major > 2:
                mdata = pickle.load(f, encoding='latin1')
            else:
                mdata = pickle.load(f)
    else:
        logging.info("calculating matlab noise...")
        gwincpath = os.path.join(os.path.dirname(__file__), 'gwinc')
        from ..gwinc_matlab import gwinc_matlab
        score, noises, ifo = gwinc_matlab(freq, ifo, gwincpath=gwincpath)
        mdata = dict(score=score, noises=noises, ifo=ifo)
        with open(matdata, 'wb') as f:
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

        if max(frac) < 0.01:
            continue

        logging.warning("Excessive difference: {}".format(name))
        logging.warning("  max: {:e}, min: {:e}".format(max(frac), min(frac)))

        diffs[name] = frac


    if args.plot:
        spec = (len(diffs), 2)
        for i, name in enumerate(diffs):
            ax = plt.subplot2grid(spec, (i, 0))
            ax.loglog(freq, np.sqrt(noises[name]), label='pygwinc')
            ax.loglog(freq, np.sqrt(mnoises[name]), label='matlab')
            ax.grid()
            ax.legend(loc='upper right')
            # ax.set_title(name)
            ax.set_ylabel(name)
    
            ax = plt.subplot2grid(spec, (i, 1))
            ax.loglog(freq, diffs[name], label=name)
            ax.grid()
            # ax.set_title(name)
    
        plt.suptitle("noises that differ by more than 1% [(mat-py)/py]")
        plt.show()


    if len(diffs) > 0:
        return 1
    return 0


if __name__ == '__main__':
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    sys.exit(main())
