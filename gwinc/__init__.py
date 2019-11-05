import os
import logging
import importlib

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
