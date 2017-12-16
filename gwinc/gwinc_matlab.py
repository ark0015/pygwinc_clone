import os
import tempfile
import scipy.io
import scipy.constants
import numpy as np

import logging

from .ifo import Struct, dictlist2recarray
from .plot import plot_noise

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


def start_matlab():
    """Start a MATLAB engine"""
    logging.info("starting matlab engine...")
    eng = matlab.engine.start_matlab()
    return eng


def ifo_add_constants(ifo):
    """Add "constants" sub-Struct to ifo Struct

    This is required by MATLAB gwinc.

    """
    #ifo.Constants.E0      = 8.8541878176e-12                 # F / m; Permittivity of Free Space
    ifo.Constants.E0      = scipy.constants.epsilon_0
    #ifo.Constants.hbar    = 1.054572e-34                     # J - s; (Plancks constant) / (2 * pi)
    ifo.Constants.hbar    = scipy.constants.hbar
    #ifo.Constants.kB      = 1.380658e-23                     # J / K; Boltzman Constant
    ifo.Constants.kB      = scipy.constants.k
    #ifo.Constants.h       = ifo.Constants.hbar * 2 * pi      # J - s; Planks constant
    #ifo.Constants.R       = 8.31447215                       # J / (K * mol); Gas Constant
    #ifo.Constants.m_e     = 9.10938291e-31                   # kg; electron mass
    #ifo.Constants.c       = 2.99792458e8                     # m / s; speed of light in vacuum
    ifo.Constants.c       = scipy.constants.c
    ifo.Constants.yr      = 365.2422 * 86400                 # sec; Seconds in a year
    #ifo.Constants.M_earth = 5.972e24                         # mass of Earth [kg]
    ifo.Constants.R_earth = 6.3781e6                         # radius of Earth [m]
    #ifo.Constants.fs      = 16384                            # Sampling frequency (Hz)
    ifo.Constants.AU      = 149597870700                     # m; Astronomical unit, IAU 2012 Resolution B2
    ifo.Constants.parsec  = ifo.Constants.AU * (648000 / np.pi) # m, IAU 2015 Resolution B2
    ifo.Constants.Mpc     = ifo.Constants.parsec * 1e6       # m, IAU 2015 Resolution B2
    ifo.Constants.SolarMassParameter = 1.3271244e20          # m^3 / s^2; G * MSol, IAU 2015 Resolution B3
    #ifo.Constants.G       = 6.67408e-11                      # m^3 / (kg  s^2); Grav. const
    ifo.Constants.G       = scipy.constants.G
    #                                                          # http://arxiv.org/abs/1507.07956
    ifo.Constants.MSol    = ifo.Constants.SolarMassParameter / ifo.Constants.G # kg; Solar mass
    #ifo.Constants.g       = 9.806                            # m / s^2; grav. acceleration
    ifo.Constants.g       = scipy.constants.g
    #                                                          # http://physics.nist.gov/cuu/Constants/ 
    #ifo.Constants.H0      = 67110;                            # ms^( - 1); Hubble const.
    #                                                          # http://arxiv.org/pdf/1303.5076v3.pdf
    #ifo.Constants.omegaM  = 0.3175;                           # Mass density parameter 
    #                                                          # http://arxiv.org/pdf/1303.5076v3.pdf
    #ifo.Constants.omegaLambda = 1 - ifo.Constants.omegaM;     # Cosmological constant density parameter
    #                                                          # omegaK = 0 (flat universe) is assumed
    return ifo


NOISE_NAME_MAP = {
    'Quantum': 'Quantum Vacuum',
    'Newtonian': 'Newtonian Gravity',
    'CoatBrown': 'Coating Brownian',
    'CoatTO': 'Coating Thermo-Optic',
    'SubBrown': 'Substrate Brownian',
    'SubTE': 'Substrate Thermo-Elastic',
    'SuspThermal': 'Suspension Thermal',
    'ResGas': 'Excess Gas',
}
def _rename_noises(d):
    nd = {}
    for k,v in d.iteritems():
        try:
            nk = NOISE_NAME_MAP[k]
        except KeyError:
            nk = k
        if isinstance(v, dict):
            nd[nk] = _rename_noises(v)
        else:
            nd[nk] = v
    return nd


def gwinc_matlab(f, ifo, fig=False, title=None, gwincpath=None, eng=None, **kwargs):
    """Execute gwinc in MATLAB with the specified ifo model.

    This uses the python matlab.engine (see start_matlab()) to
    calculate noises with MATLAB gwinc at the specified path
    `gwincpath` (defaults to 'gwinc' in the current directory if not
    specified, or uses the GWINCPATH environment variable).

    returns [source, noises, ifo] as returned from matlab.

    """
    if not gwincpath:
        gwincpath = os.getenv('GWINCPATH', 'gwinc')
    if not os.path.exists(os.path.join(gwincpath, 'gwinc.m')):
        raise IOError("Invalid matlab gwinc path: '{}'".format(gwincpath))

    if not eng:
        eng = start_matlab()

    logging.info("adding gwinc path: {}".format(gwincpath))
    eng.addpath(gwincpath)

    # add Constants attribute to ifo structure
    ifo_add_constants(ifo)

    # this stupidity because you can't just give the engine a np.ndarray
    eng.workspace['f'] = f.tolist()
    eng.eval('f = cell2mat(f);', nargout=0)

    # similar stupidity prevents this (this time recarrays in the dict):
    #eng.workspace['ifo'] = ifo.to_dict(array=True)
    with tempfile.NamedTemporaryFile(suffix='.mat') as f:
        scipy.io.savemat(f, ifo.to_dict(array=True))
        eng.eval("ifo = load('{}');".format(f.name), nargout=0)

    eng.eval('[score, noises, ifo] = gwinc(f, [], ifo, [], 0);', nargout=0)

    with tempfile.NamedTemporaryFile(suffix='.mat') as f:
        eng.save(f.name, 'score', 'noises', 'ifo', nargout=0)
        data = scipy.io.loadmat(f, squeeze_me=True, struct_as_record=False)

    score = data['score']
    mnoises = Struct.from_matstruct(data['noises']).to_dict()
    ##### blow out 'MirrorThermal' sub-dict
    for n,d in mnoises['MirrorThermal'].iteritems():
        if n == 'Total':
            continue
        mnoises[n] = d
    del mnoises['MirrorThermal']
    #####
    noises = _rename_noises(mnoises)

    return score, noises, ifo
