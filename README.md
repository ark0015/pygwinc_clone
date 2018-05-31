[![pipeline status](https://git.ligo.org/gwinc/pygwinc/badges/master/pipeline.svg)](https://git.ligo.org/gwinc/pygwinc/commits/master)

# Python port of GW Interferometer Noise Calculator

![gwinc](https://gwinc.docs.ligo.org/pygwinc/aLIGO.png)


## basic usage

`pygwinc` creates noise budgets based on detector descriptions
provided in either .yml or .mat files (see below).  Once the detector
description is loaded, the noise budget can be calculated and plotted:
```python
>>> import gwinc
>>> import numpy as np
>>> freq = np.logspace(1, 3, 1000)
>>> ifo = gwinc.load_ifo('aLIGO')
>>> ifo = gwinc.precompIFO(ifo)
>>> noises = gwinc.noise_calc(ifo, freq)
>>> gwinc.plot_noise(ifo, noises)
```
Or the `gwinc` convenience function can be used to handle it all:
```
>>> score, data, ifo = gwinc.gwinc(freq, ifo, plot=True)
```


## command line interface

You can make gwinc plots directly from the command line by executing
the package directly:
```shell
$ python3 -m gwinc -h
usage: gwinc [-h] [--flo FLO] [--fhi FHI] [--npoints NPOINTS] [--title TITLE]
             [--matlab] [--fom FOM] [--dump | --save SAVE | --interactive]
             [IFO]

Plot GWINC noise budget for specified IFO.

If the inspiral_range package is installed, various figures of merit
can be calculated for the resultant spectrum with the --fom argument,
e.g.:

  gwinc --fom horizon ...
  gwinc --fom range:m1=20,m2=20 ...

See documentation for inspiral_range package for details.

positional arguments:
  IFO                   IFO name or description file path (.yaml or .mat)

optional arguments:
  -h, --help            show this help message and exit
  --flo FLO, -fl FLO    lower frequency bound in Hz [5]
  --fhi FHI, --fh FHI   upper frequency bound in Hz [6000]
  --npoints NPOINTS, -n NPOINTS
                        number of frequency points [3000]
  --title TITLE, -t TITLE
                        plot title
  --matlab, -m          use MATLAB gwinc engine to calculate noises
  --fom FOM             calculate inspiral range for resultant spectrum
                        ('func[:param=val,param=val]')
  --dump, -d            print IFO parameters to stdout and exit
  --save SAVE, -s SAVE  save figure to file
  --interactive, -i     open interactive shell when plotting
```


## detector description files

`pygwinc` can load detector descriptions in different formats: the new
YAML .yaml format, the original MATLAB gwinc .mat format, or even from
a MATLAB .m file.  `pygwinc` includes .yaml detector descriptions for
various detectors:

* [aLIGO.yaml](https://git.ligo.org/gwinc/pygwinc/blob/master/gwinc/ifo/aLIGO.yaml)
* [A+.yaml](https://git.ligo.org/gwinc/pygwinc/blob/master/gwinc/ifo/A+.yaml)
* [Voyager.yaml](https://git.ligo.org/gwinc/pygwinc/blob/master/gwinc/ifo/Voyager.yaml)


## comparison with MATLAB gwinc

`pygwinc` includes the ability use MATLAB gwinc directly via the
MATLAB python interface (see the CLI '--matlab' option above).  This
also allows for easy direct comparison between the pygwinc and
matgwinc noise budgets.

If you have a local checkout of matgwinc (at e.g. /path/to/gwinc) and
a local installation of MATLAB and it's python interface (at
e.g. /opt/matlab/python/lib/python3.6/site-packages) you can run the
comparison as so:
```shell
$ export GWINCPATH=/path/to/matgwinc
$ export PYTHONPATH=/opt/matlab/python/lib/python3.6/site-packages
$ python3 -m gwinc.test -p aLIGO
```
This will produce a summary page of the various noise spectra that
differ between matgwinc and pygwinc.

Latest comparison plots from continuous integration:

* [aLIGO comparison](https://gwinc.docs.ligo.org/pygwinc/aLIGO_test.png)
* [A+ comparison](https://gwinc.docs.ligo.org/pygwinc/A+_test.png)
