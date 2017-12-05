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
