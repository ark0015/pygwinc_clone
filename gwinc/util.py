import datetime
import h5py

from . import struct


def load_hdf5(path):
    """Load IFO and noises from HDF5 file.

    Returns (name, ifo, noises) tuple load from file.

    """
    with h5py.File(path, 'r') as f:
        name = f.attrs['name']
        ifo = struct.Struct.from_yaml(f.attrs['IFO'])
        noises = {}
        for name, data in f['/traces'].items():
            noises[name] = data.value
        return name, ifo, noises


def save_hdf5(name, ifo, noises, path):
    """Save IFO and noises to HDF5 file.

    """
    with h5py.File(path, 'w') as f:
        f.attrs['schema'] = 'GWINC noise budget'
        # FIXME: add GWINC version
        f.attrs['date'] = datetime.datetime.now().isoformat()
        f.attrs['name'] = name
        f.attrs['IFO'] = ifo.to_yaml()
        tgrp = f.create_group('traces')
        for noise, data in noises.items():
            tgrp.create_dataset(noise, data=data)
