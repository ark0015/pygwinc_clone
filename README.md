[![pipeline status](https://git.ligo.org/gwinc/pygwinc/badges/master/pipeline.svg)](https://git.ligo.org/gwinc/pygwinc/commits/master)

# Python port of GW Interferometer Noise Calculator

![gwinc](http://gwinc.docs.ligo.org/pygwinc/aLIGO.png)

## tests

To compare pygwinc vs. matlab gwinc, run the included test command
(requires a local installation of matlab and its python interface):

*   checkout, or link an existing checkout of, matlab gwinc in
    the test directory:

        $ cd pygwinc/gwinc/test
        $ ln -s ~/ligo/src/iscmodeling/gwinc

*   Execute the test command (specifying path to matlab python
    interface if needed):

        $ PYTHONPATH=/opt/matlab/python/lib/python2.7/site-packages python -m gwinc.test

This will produce difference plots of noises calculated from matlab
gwinc and pygwinc.

## basic usage

`pygwinc` creates noise budgets based on detector descriptions
provided in either .yml or .mat files (see below).  Once the detector
description is loaded, the noise budget can be calculated with the
`gwinc` command:
```python
>>> import gwinc
>>> import numpy as np
>>> ifo = gwinc.load_ifo('aLIGO')
>>> freq = np.logspace(1, 3, 1000)
>>> score, data, ifo = gwinc.gwinc(freq, ifo)
```
A convenience function to plot the resulting noise budget is included:
```
>>> gwinc.plot_noise(data)
```

## command line interface

You can make gwinc plots directly from the command line by executing
the package directly:
```shell
~/ligo/src/gwinc $ python -m gwinc -h
usage: gwinc [-h] [--flo FLO] [--fhi FHI] [--npoints NPOINTS] [--title TITLE]
             [--matlab] [--fom FOM] [--dump | --save SAVE | --interactive]
             [IFO]

Plot GWINC noise budget for specified IFO

If the inspiral_range package is installed, various figures of merit
can be calculated for the resultant spectrum with the --fom argument,
e.g.:

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
                        ('func:param=val,param=val')
  --dump, -d            print IFO parameters to stdout and exit
  --save SAVE, -s SAVE  save figure to file
  --interactive, -i     open interactive shell when plotting
```

## detector description files

`pygwinc` can load detector descriptions in two different format: the
original MATLAB gwinc .mat format, or the new YAML .yaml format.
`pygwinc` includes two .yaml detector descriptions:

    * gwinc/ifo/aLIGO.yaml
    * gwinc/ifo/Voyager.yaml
