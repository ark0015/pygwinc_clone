import os
import tempfile
import scipy.io
import scipy.constants
import numpy as np
import logging

from .struct import Struct

##################################################

MATLAB_ENGINE = None

class Matlab:
    def __init__(self, gwincpath=None):
        """Start a MATLAB engine for GWINC processing

         engine is provided, the GWINC path is added
        to it's path.

        """
        global MATLAB_ENGINE

        if MATLAB_ENGINE:
            return

        import matlab.engine

        if not gwincpath:
            gwincpath = os.path.expanduser(os.getenv('GWINCPATH', 'gwinc'))
        if not os.path.exists(os.path.join(gwincpath, 'gwinc.m')):
            raise IOError("Invalid MATLAB GWINC path: '{}'".format(gwincpath))

        logging.info("starting MATLAB engine...")
        MATLAB_ENGINE = matlab.engine.start_matlab()

        logging.info("gwinc path: {}".format(gwincpath))
        MATLAB_ENGINE.addpath(gwincpath)


    @property
    def eng(self):
        return MATLAB_ENGINE

    @property
    def workspace(self):
        return MATLAB_ENGINE.workspace

    def addpath(self, *args):
        return MATLAB_ENGINE.addpath(*args)

    def eval(self, *args, **kwargs):
        return MATLAB_ENGINE.eval(*args, **kwargs)


    def load_array(self, var, array):
        """Load numpy array into workspace as vector.

        `var` is name of workspace variable, and `array` is numpy
        ndarray.

        """
        # this stupidity because you can't just give the engine a np.ndarray
        MATLAB_ENGINE.workspace[var] = array.tolist()
        MATLAB_ENGINE.eval('{0} = cell2mat({0});'.format(var), nargout=0)


    def load_struct(self, var, struct):
        """Load pygwinc.Struct array into workspace as vector.

        `var` is name of workspace variable, and `struct` is
        pygwinc.Struct.

        """
        # similar stupidity prevents this (this time recarrays in the dict):
        #matlab.workspace['ifo'] = ifo.to_dict(array=True)
        with tempfile.NamedTemporaryFile(suffix='.mat') as f:
            scipy.io.savemat(f, struct.to_dict(array=True))
            MATLAB_ENGINE.eval("{} = load('{}');".format(var, f.name), nargout=0)


    def extract(self, *wvars):
        """Extract workspace variables from engine.

        Returns dict with wvars as keys.

        """
        assert len(wvars) > 0
        with tempfile.NamedTemporaryFile(suffix='.mat') as f:
            MATLAB_ENGINE.save(f.name, *wvars, nargout=0)
            data = scipy.io.loadmat(f, squeeze_me=True, struct_as_record=False)
        if len(wvars) == 1:
            return data[wvars[0]]
        else:
            return data

##################################################

def ifo_add_constants(ifo):
    """Add "constants" sub-Struct to ifo Struct

    This is required by MATLAB gwinc.

    """
    # permittivity of free space [F / m]
    #ifo.Constants.E0      = 8.8541878176e-12
    ifo.Constants.E0      = scipy.constants.epsilon_0
    # Plancks constant [J s]
    #ifo.Constants.hbar    = 1.054572e-34
    ifo.Constants.hbar    = scipy.constants.hbar
    # Planks constant [J s]
    #ifo.Constants.h       = ifo.Constants.hbar * 2 * pi
    # Boltzman constant [J / K]
    #ifo.Constants.kB      = 1.380658e-23
    ifo.Constants.kB      = scipy.constants.k
    # gas constant [J / (K * mol)]
    #ifo.Constants.R       = 8.31447215
    # electron mass [kg]
    #ifo.Constants.m_e     = 9.10938291e-31
    # speed of light in vacuum [m / s]
    #ifo.Constants.c       = 2.99792458e8
    ifo.Constants.c       = scipy.constants.c
    # seconds in a year [s]
    ifo.Constants.yr      = 365.2422 * 86400
    # mass of Earth [kg]
    #ifo.Constants.M_earth = 5.972e24
    # radius of Earth [m]
    ifo.Constants.R_earth = 6.3781e6
    # sampling frequency [Hz]
    #ifo.Constants.fs      = 16384
    # Astronomical unit, IAU 2012 Resolution B2 [m]
    ifo.Constants.AU      = 149597870700
    # IAU 2015 Resolution B2 [m]
    ifo.Constants.parsec  = ifo.Constants.AU * (648000 / np.pi)
    # IAU 2015 Resolution B2 [m]
    ifo.Constants.Mpc     = ifo.Constants.parsec * 1e6
    # IAU 2015 Resolution B3 [m^3 / s^2; G * MSol]
    ifo.Constants.SolarMassParameter = 1.3271244e20
    # gravitational const [m^3 / (kg  s^2)]
    #ifo.Constants.G       = 6.67408e-11
    ifo.Constants.G       = scipy.constants.G
    # solar mass [kg]
    # http://arxiv.org/abs/1507.07956
    ifo.Constants.MSol    = ifo.Constants.SolarMassParameter / ifo.Constants.G
    # gravitational acceleration [m / s^2]
    #ifo.Constants.g       = 9.806
    ifo.Constants.g       = scipy.constants.g
    # Hubble constant [ms^( - 1)]
    # http://physics.nist.gov/cuu/Constants/
    ifo.Constants.H0      = 67110
    # http://arxiv.org/pdf/1303.5076v3.pdf
    # Mass density parameter
    ifo.Constants.omegaM  = 0.3175
    # http://arxiv.org/pdf/1303.5076v3.pdf
    # Cosmological constant density parameter
    # omegaK = 0 (flat universe) is assumed
    ifo.Constants.omegaLambda = 1 - ifo.Constants.omegaM
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
    for k,v in d.items():
        try:
            nk = NOISE_NAME_MAP[k]
        except KeyError:
            nk = k
        if isinstance(v, dict):
            nd[nk] = _rename_noises(v)
        else:
            nd[nk] = v
    return nd


def gwinc_matlab(f, ifo, plot=False):
    """Execute gwinc in MATLAB with the specified ifo model.

    This uses the python matlab.engine (see Matlab class) to calculate
    noises with MATLAB gwinc (gwinc directory specified by GWINCPATH
    environment variable).

    Returns `score` (dict), `noises` (dict), and `ifo` (Struct) as
    returned from MATLAB.

    If `plot` is True will cause MATLAB to produce it's own plot for
    the noise budget.

    """
    matlab = Matlab()

    # add Constants attribute to ifo structure
    ifo_add_constants(ifo)

    matlab.load_array('f', f)
    matlab.load_struct('ifo', ifo)

    if plot:
        plot_flag = '3'
    else:
        plot_flag = '0'

    cmd = "[score, noises, ifo] = gwinc(f, [], ifo, SourceModel, {});".format(plot_flag)
    matlab.eval(cmd, nargout=0)

    data = matlab.extract('score', 'noises', 'ifo')

    score = data['score']
    mnoises = Struct.from_matstruct(data['noises']).to_dict()
    ##### blow out 'MirrorThermal' sub-dict
    for n,d in mnoises['MirrorThermal'].items():
        if n == 'Total':
            continue
        mnoises[n] = d
    del mnoises['MirrorThermal']
    #####
    noises = _rename_noises(mnoises)
    ifo = Struct.from_matstruct(data['ifo'])

    return score, noises, ifo
