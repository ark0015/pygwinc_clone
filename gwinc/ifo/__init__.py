import os
import re
import yaml
import scipy
import scipy.io
import numpy as np
from scipy.io.matlab.mio5_params import mat_struct


# HACK: fix loading number in scientific notation
#
# https://stackoverflow.com/questions/30458977/yaml-loads-5e-6-as-string-and-not-a-number
#
# An apparent bug in python-yaml prevents it from regognizing
# scientific notation as a float.  The following is a modified version
# of the parser that recognize scientific notation appropriately.
#
loader = yaml.SafeLoader
loader.add_implicit_resolver(
    u'tag:yaml.org,2002:float',
    re.compile(u'''^(?:
     [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
    |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
    |\\.[0-9_]+(?:[eE][-+][0-9]+)?
    |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
    |[-+]?\\.(?:inf|Inf|INF)
    |\\.(?:nan|NaN|NAN))$''', re.X),
    list(u'-+0123456789.'))


def dictlist2recarray(l):
    def dtype(v):
        if isinstance(v, int):
            return float
        else:
            return type(v)
    # get dtypes from first element dict
    dtypes = [(k, dtype(v)) for k,v in l[0].iteritems()]
    values = [tuple(el.values()) for el in l]
    out = np.array(values, dtype=dtypes)
    return out.view(np.recarray)


class Struct(object):
    """Matlab struct-like object

    This is a simple implementation of a MATLAB struct-like object
    that stores values as attributes of a simple class: and allows assigning to
    attributes recursively, e.g.:

    >>> s = Struct()
    >>> s.a = 4
    >>> s.b = Struct()
    >>> s.b.c = 8

    Various classmethods allow creating one of these objects from YAML
    file, a nested dict, or a MATLAB struct object.

    """

    # FIXME: This would be a way to allow setting nested struct
    # attributes, e.g.:
    #
    # >>> s = Struct()
    # >>> s.a.b.c = 4
    #
    # Usage of __getattr__ like this is dangerous and creates
    # non-intuitive behavior (i.e. an empty struct is returned with
    # accessing attributes that don't exist).  Is there a way to
    # accomplish this without that adverse side affect?
    #
    # def __getattr__(self, name):
    #     if name not in self.__dict__:
    #         self.__dict__[name] = Struct()
    #     return self.__dict__[name]

    ##########

    def __contains__(self, item):
        return item in self.__dict__

    def to_dict(self, array=False):
        """Return nested dictionary representation of Struct.

        If `array` is True any lists encountered will be turned into
        numpy arrays, and lists of Structs will be turned into record
        arrays.  This is need to convert to structure arrays in
        matlab.

        """
        d = {}
        for k,v in self.__dict__.iteritems():
            if isinstance(v, Struct):
                d[k] = v.to_dict(array=array)
            else:
                if isinstance(v, list):
                    try:
                        # this should fail if the elements of v are
                        # not Struct
                        # FIXME: need cleaner way to do this
                        v = [i.to_dict(array=array) for i in v]
                        if array:
                            v = dictlist2recarray(v)
                    except AttributeError:
                        if array:
                            v = np.array(v)
                elif isinstance(v, int):
                    v = float(v)
                d[k] = v
        return d

    def to_yaml(self, path=None):
        """Return YAML representation of Struct as nested dict.

        Or write Struct to YAML file if file 'path' argument
        specified.

        """
        y = yaml.dump(self.to_dict(), default_flow_style=False)
        if path:
            with open(path, 'w') as f:
                f.write(y)
        else:
            return y

    # def __repr__(self):
    #     return self.to_yaml().strip('\n')

    def __str__(self):
        return '<GWINC Struct: {}>'.format(self.__dict__.keys())

    def __iter__(self):
        return iter(self.__dict__)

    def walk(self):
        """Iterate over all leaves in the struct tree.

        """
        for k,v in self.__dict__.iteritems():
            if type(v) is Struct:
                for sk,sv in v.walk():
                    yield k+'.'+sk, sv
            else:
                try:
                    for i,vv in enumerate(v):
                        for sk,sv in vv.walk():
                            yield '{}[{}].{}'.format(k,i,sk), sv
                except (AttributeError, TypeError):
                    yield k, v

    @classmethod
    def from_dict(cls, d):
        """Create Struct from nested dict.

        """
        c = cls()
        for k,v in d.iteritems():
            if type(v) == dict:
                c.__dict__[k] = Struct.from_dict(v)
            else:
                try:
                    c.__dict__[k] = map(Struct.from_dict, v)
                except (AttributeError, TypeError):
                    c.__dict__[k] = v
        return c

    @classmethod
    def from_matstruct(cls, s):
        """Create Struct from scipy.io.matlab mat_struct object.

        """
        c = cls()
        try:
            s = s['ifo']
        except:
            pass
        for k,v in s.__dict__.iteritems():
            if k in ['_fieldnames']:
                # skip these fields
                pass
            elif type(v) is mat_struct:
                c.__dict__[k] = Struct.from_matstruct(v)
            else:
                # handle lists of Structs
                try:
                    c.__dict__[k] = map(Struct.from_matstruct, v)
                except:
                    c.__dict__[k] = v
        return c
                

    @classmethod
    def from_file(cls, path):
        """Load Struct from .yaml or GWINC .mat file.

        File type will be determined by extension.

        """
        (root, ext) = os.path.splitext(path)
        with open(path, 'r') as f:
            if ext in ['.yaml', '.yml']:
                d = yaml.load(f, Loader=loader)
                return cls.from_dict(d)
            elif ext == '.mat':
                s = scipy.io.loadmat(f, squeeze_me=True, struct_as_record=False)
                return cls.from_matstruct(s)
            else:
                raise IOError("Unknown file type: {}".format(ext))



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


##################################################


if __name__ == '__main__':
    import sys
    ifo = load_ifo(sys.argv[1])
    # print(ifo.to_yaml())
    print(ifo.to_dict())
