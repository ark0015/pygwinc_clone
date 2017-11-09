# Python port of GW Interferometer Noise Calculator

![gwinc](https://git.ligo.org/gwinc/pygwinc/uploads/f324c0ce516db509a4028858faed3891/gwinc.png)

## tests

To compare pygwinc vs. MATLAB gwinc, run the included test command
(requires a local installation of MATLAB and it's python interface):

*   checkout, or create a link existing checkout of, matlab gwinc in
    the test directory:

        $ cd pygwin/gwinc/test
        $ ln -s ~/ligo/src/iscmodeling/gwinc

*   Execute the test command (specifying path to MATLAB python
    interface if needed):

        $ PYTHONPATH=/opt/matlab/python/lib/python2.7/site-packages python -m gwinc.test

This will produce difference plots of noises calculated from MATLAB
gwinc and pygwinc.
