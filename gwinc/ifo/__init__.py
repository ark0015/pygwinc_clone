import os
import fnmatch

from ..struct import Struct

def available_ifos():
    """List available included IFO files"""
    ifos = []
    for f in os.listdir(os.path.dirname(__file__)):
        if fnmatch.fnmatch(f, '*.yaml'):
            ifos.append(f.split('.')[0])
    return ifos


def load_ifo(name_or_path):
    """Load IFO by name or from file.

    IFO names will correspond to basename of included .yaml IFO
    definition file.

    When specifying path may be either .yaml or .mat.

    """
    if os.path.exists(name_or_path):
        path = name_or_path
    else:
        path = os.path.join(os.path.dirname(__file__),
                            name_or_path+'.yaml')
    s = Struct.from_file(path)
    return s
