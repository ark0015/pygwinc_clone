import os
import fnmatch

from ..struct import Struct
from ..gwinc_matlab import Matlab


def available_ifos():
    """List available included IFO files"""
    ifos = []
    for f in os.listdir(os.path.dirname(__file__)):
        if fnmatch.fnmatch(f, '*.yaml'):
            ifos.append(f.split('.')[0])
    return ifos


def load_ifo(name_or_path):
    """Load IFO by name or from file.

    Named IFOs should correspond to the basename of .yaml IFO
    definition files included with pygwinc (see available_ifos()
    above).

    When specifying by path files may be either .yaml, .mat or .m.
    For .m files, the file is expected to include either an object or
    function that corresponds to the basename of the file.  The MATLAB
    engine will be invoked to execute the .m code and extract the
    resultant IFO data.

    """
    if os.path.exists(name_or_path):
        path = name_or_path
    else:
        path = os.path.join(os.path.dirname(__file__),
                            name_or_path+'.yaml')

    (root, ext) = os.path.splitext(path)

    if ext == '.m':
        matlab = Matlab()
        matlab.addpath(os.path.dirname(path))
        func_name = os.path.basename(root)
        matlab.eval("ifo = {};".format(func_name), nargout=0)
        ifo = matlab.extract('ifo')
        return Struct.from_matstruct(ifo)

    else:
        return Struct.from_file(path)
