"""Beginnings of a test suite for pygwinc

todo:

 * unittest
 * write test to compare matgwinc default IFOs and pygwinc default

"""
import os
import signal
import scipy
import scipy.io
import numpy as np
import matplotlib.pyplot as plt

import logging
logging.basicConfig(level=logging.INFO)

from .. import load_ifo, noise_calc
from ..plot import plot_noise
from ..ifo import Struct
# NOTE: gwinc needs to be imported before matlab for some very strange
# reason:
#
# Traceback (most recent call last):
#   File "./ifo2txt.py", line 9, in <module>
#     import gwinc
#   File "/home/jrollins/ligo/src/pygwinc/gwinc/__init__.py", line 9, in <module>
#     from . import plot
#   File "/home/jrollins/ligo/src/pygwinc/gwinc/plot.py", line 2, in <module>
#     import matplotlib.pyplot as plt
#   File "/usr/lib/python2.7/dist-packages/matplotlib/pyplot.py", line 29, in <module>
#     import matplotlib.colorbar
#   File "/usr/lib/python2.7/dist-packages/matplotlib/colorbar.py", line 32, in <module>
#     import matplotlib.artist as martist
#   File "/usr/lib/python2.7/dist-packages/matplotlib/artist.py", line 15, in <module>
#     from .transforms import (Bbox, IdentityTransform, TransformedBbox,
#   File "/usr/lib/python2.7/dist-packages/matplotlib/transforms.py", line 39, in <module>
#     from matplotlib._path import (affine_transform, count_bboxes_overlapping_bbox,
# ImportError: /opt/matlab/R2016b/extern/engines/python/dist/matlab/engine/glnxa64/../../../../../../../sys/os/glnxa64/libstdc++.so.6: version `CXXABI_1.3.9' not found (required by /usr/lib/python2.7/dist-packages/matplotlib/_path.x86_64-linux-gnu.so)
#
# no idea what that's about
import matlab.engine


GWINC_MAT_SCRIPT = '''
source = SourceModel();
source.BlackHole.Mass1 = 20; % [MSol]
source.BlackHole.Mass2 = 20; % [MSol]

flo = 5;
fhi = 6000;
npoints = 3000;
f = logspace(log10(flo), log10(fhi), npoints);

ifo = {ifoname}();
ifo = precompIFO(ifo);
[score, noises, ifo] = gwinc(f, [], ifo, source, 0);
'''

def run_gwinc(eng, ifoname, matpath):
    # This doesn't work I think because the ifo.Suspension.Stages is not a
    # "scalar structure":
    # ifo = eng.IFOModel()
    # ifo = eng.precompIFO(ifo)
    script = GWINC_MAT_SCRIPT.format(ifoname=ifoname)
    print()
    print(script)
    print()
    eng.eval(script, nargout=0)
    eng.save(matpath, 'score', 'noises', 'ifo', nargout=0)


def main():

    ifoname = 'IFOModel'
    gwincpath = os.path.join(os.path.dirname(__file__),
                             'gwinc')
    logging.info('gwincpath: {}'.format(gwincpath))
    matpath = os.path.join(os.path.dirname(__file__),
                           'gwinc.mat')
    logging.info('matpath: {}'.format(matpath))

    if not os.path.exists(matpath):
        logging.info("starting matlab engine...")
        eng = matlab.engine.start_matlab()
        eng.addpath(gwincpath)
        logging.info("running gwinc matlab script...")

        run_gwinc(eng, ifoname, matpath)
        eng.exit()

    else:
        logging.info("Using existing {}...".format(matpath))

    logging.info("loading {}...".format(matpath))
    mat_dict = scipy.io.loadmat(matpath,
                                squeeze_me=True,
                                struct_as_record=False)

    ifo = Struct.from_matstruct(mat_dict['ifo'])
    matnoises = Struct.from_matstruct(mat_dict['noises']).to_dict()
    f = matnoises['Freq']

    # FIXME: pygwinc expects some newer attributes, for newer
    # suspensionthermal calc that requires stage dilution, among other
    # things.  Should make sure it's fully backwards compatible.
    ifo.Seismic.Omicron = 1
    ifo.Suspension.FiberType = 0
    ifo.Suspension.Stage[0].Dilution = np.NaN
    ifo.Suspension.Stage[1].Dilution = 106
    ifo.Suspension.Stage[2].Dilution = 80
    ifo.Suspension.Stage[3].Dilution = 87

    logging.info("calculating noise...")
    noises = noise_calc(ifo, f)

    # plt.figure()
    # plot_noise(noises)
    # plt.title('pygwinc with matgwinc ifo')
    # plt.show(block=False)

    # mapping from pygwinc names to matgwinc names
    MAPPING = {
        'Quantum Vacuum': 'Quantum',
        'Seismic': 'Seismic',
        'Newtonian Gravity': 'Newtonian',
        'Suspension Thermal': 'SuspThermal',
        'Coating Brownian': 'MirrorThermal',
        # 'Coating Thermo-Optic': 
        # 'ITM Thermo-Refractive':
        # 'ITM Carrier Density':
        'Substrate Brownian': 'SubBrown',
        'Substrate Thermo-Elastic': 'SubTE',
        'Excess Gas': 'ResGas',
        }

    diffs = {}
    for name, noise in noises.iteritems():
        try:
            matnoise = matnoises[MAPPING[name]]
        except KeyError:
            continue
        logging.info("plot compare: {}".format(name))
        # FIXME: for CoatTO
        if isinstance(matnoise, dict):
            logging.info("  skip")
            continue
        diff = matnoise - noise
        logging.info("  max: {}, min: {}".format(max(diff), min(diff)))
        diffs[name] = diff

    fig, ax = plt.subplots(len(diffs), sharex=True)
    for i, name in enumerate(diffs):
        diff = diffs[name]
        ax[i].semilogx(f, diff, label=name)
        ax[i].legend(loc='lower right')
    ax[0].set_title('difference: (mat_gwinc - pygwinc)')
    plt.xlabel('frequency [Hz]')
    plt.show()


if __name__ == '__main__':
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    main()
