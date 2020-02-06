from __future__ import division
import os
import logging
import importlib
import numpy as np

from .ifo import available_ifos
from .struct import Struct
from .precomp import precompIFO
from .plot import plot_noise
from .io import load_hdf5, save_hdf5


def _load_module(name_or_path):
    """Load module from name or path.

    Return loaded module and module path.

    """
    if os.path.exists(name_or_path):
        path = name_or_path.rstrip('/')
        modname = os.path.splitext(os.path.basename(path))[0]
        if os.path.isdir(path):
            path = os.path.join(path, '__init__.py')
        spec = importlib.util.spec_from_file_location(modname, path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    else:
        mod = importlib.import_module(name_or_path)
    try:
        path = mod.__path__[0]
    except AttributeError:
        path = mod.__file__
    return mod, path


def load_budget(name_or_path):
    """Load GWINC IFO Budget by name or from path.

    Named IFOs should correspond to one of the IFOs available in the
    gwinc package (see gwinc.available_ifos()).  If a path is provided
    it should either be a budget package (directory) or module (ending
    in .py), or an IFO struct definition (see gwinc.Struct).

    If a budget package path is provided, and the package includes an
    'ifo.yaml' file, that file will be loaded into a Struct and
    assigned as an attribute to the returned Budget class.

    If a bare struct is provided the base aLIGO budget definition will
    be assumed.

    Returns a Budget class with 'ifo', 'freq', and 'plot_style', ifo
    structure, frequency array, and plot style dictionary, with the
    last three being None if they are not defined in the budget.

    """
    ifo = None

    if os.path.exists(name_or_path):
        path = name_or_path.rstrip('/')
        bname, ext = os.path.splitext(os.path.basename(path))

        if ext in Struct.STRUCT_EXT:
            logging.info("loading struct {}...".format(path))
            ifo = Struct.from_file(path)
            bname = 'aLIGO'
            modname = 'gwinc.ifo.aLIGO'
            logging.info("loading budget {}...".format(modname))

        else:
            modname = path
            logging.info("loading module path {}...".format(modname))

    else:
        if name_or_path not in available_ifos():
            raise RuntimeError("Unknonw IFO '{}' (available IFOs: {}).".format(
                name_or_path,
                available_ifos(),
            ))
        bname = name_or_path
        modname = 'gwinc.ifo.'+name_or_path
        logging.info("loading module {}...".format(modname))

    mod, modpath = _load_module(modname)

    Budget = getattr(mod, bname)
    ifopath = os.path.join(modpath, 'ifo.yaml')
    if not ifo and ifopath:
        ifo = Struct.from_file(ifopath)
    Budget.ifo = ifo

    return Budget


def gwinc(freq, ifo, source=None, plot=False, PRfixed=True):
    """Calculate strain noise budget for a specified interferometer model.

    Argument `freq` is the frequency array for which the noises will
    be calculated, and `ifoin` is the IFO model (see the `load_ifo()`
    function).

    If `source` structure provided, so evaluates the sensitivity of
    the detector to several potential gravitational wave
    sources.

    If `plot` is True a plot of the budget will be created.

    Returns tuple of (score, noises, ifo)

    """
    # assume generic aLIGO configuration
    # FIXME: how do we allow adding arbitrary addtional noise sources
    # from just ifo description, without having to specify full budget
    Budget = load_budget('aLIGO')
    ifo = precompIFO(freq, ifo, PRfixed)
    traces = Budget(freq, ifo=ifo).calc_trace()
    plot_style = getattr(Budget, 'plot_style', {})

    # construct matgwinc-compatible noises structure
    noises = {}
    for name, (data, style) in traces.items():
        noises[style.get('label', name)] = data
    noises['Freq'] = freq

    pbs = ifo.gwinc.pbs
    parm = ifo.gwinc.parm
    finesse = ifo.gwinc.finesse
    prfactor = ifo.gwinc.prfactor
    if ifo.Laser.Power * prfactor != pbs:
        pass
        #warning(['Thermal lensing limits input power to ' num2str(pbs/prfactor, 3) ' W']);

    # report astrophysical scores if so desired
    score = None
    if source:
        score = int73(freq, noises['Total'], ifo, source)
        score.Omega = intStoch(freq, noises['Total'], 0, ifo, source)

    # --------------------------------------------------------
    # output graphics

    if plot:
        logging.info('Laser Power:            %7.2f Watt' % ifo.Laser.Power)
        logging.info('SRM Detuning:           %7.2f degree' % (ifo.Optics.SRM.Tunephase*180/pi))
        logging.info('SRM transmission:       %9.4f' % ifo.Optics.SRM.Transmittance)
        logging.info('ITM transmission:       %9.4f' % ifo.Optics.ITM.Transmittance)
        logging.info('PRM transmission:       %9.4f' % ifo.Optics.PRM.Transmittance)
        logging.info('Finesse:                %7.2f' % finesse)
        logging.info('Power Recycling Gain:   %7.2f' % prfactor)
        logging.info('Arm Power:              %7.2f kW' % (parm/1000))
        logging.info('Power on BS:            %7.2f W' % pbs)

        # coating and substrate thermal load on the ITM
        PowAbsITM = (pbs/2) * \
                    np.hstack([(finesse*2/np.pi) * ifo.Optics.ITM.CoatingAbsorption,
                               (2 * ifo.Materials.MassThickness) * ifo.Optics.ITM.SubstrateAbsorption])

        logging.info('Thermal load on ITM:    %8.3f W' % sum(PowAbsITM))
        logging.info('Thermal load on BS:     %8.3f W' %
                     (ifo.Materials.MassThickness*ifo.Optics.SubstrateAbsorption*pbs))
        if (ifo.Laser.Power*prfactor != pbs):
            logging.info('Lensing limited input power: %7.2f W' % (pbs/prfactor))

        if source:
            logging.info('BNS Inspiral Range:     ' + str(score.effr0ns) + ' Mpc/ z = ' + str(score.zHorizonNS))
            logging.info('BBH Inspiral Range:     ' + str(score.effr0bh) + ' Mpc/ z = ' + str(score.zHorizonBH))
            logging.info('Stochastic Omega: %4.1g Universes' % score.Omega)

        plot_noise(ifo, traces, **plot_style)

    return score, noises, ifo
