[![pipeline status](https://git.ligo.org/gwinc/pygwinc/badges/master/pipeline.svg)](https://git.ligo.org/gwinc/pygwinc/commits/master)

# Python Gravitational Wave Interferometer Noise Calculator

![gwinc](https://gwinc.docs.ligo.org/pygwinc/aLIGO.png)

`pygwinc` is a multi-faceted tool for processing and plotting noise
budgets for ground-based gravitational wave detectors.  It's primary
feature is a collection of mostly analytic noise calculation functions
for various sources of noise affecting detectors (see [noise
functions](#noise-functions) below):

* quantum noise
* mirror coating thermal noise
* mirror substrate thermal noise
* suspension fiber thermal noise
* seismic noise
* Newtonian/gravity-gradient noise
* residual gas noise

`pygwinc` is also a generalized noise budgeting tool, allowing users
to create arbitrary noise budgets (for any experiment, not just
ground-based GW detectors) using measured or analytically calculated
data.  See the [budget interface](#budget-interface) section below.

`pygwinc` includes canonical budgets for various well-known current
and future detectors:

* [aLIGO](https://gwinc.docs.ligo.org/pygwinc/aLIGO.png)
* [A+](https://gwinc.docs.ligo.org/pygwinc/Aplus.png)
* [Voyager](https://gwinc.docs.ligo.org/pygwinc/Voyager.png)
* [Cosmic Explorer 1](https://gwinc.docs.ligo.org/pygwinc/CE1.png)
* [Cosmic Explorer 2](https://gwinc.docs.ligo.org/pygwinc/CE2.png)

See [IFO.md](IFO.md) for the latest CI-generated plots and hdf5 cached
data.

For calculating various common "inspiral range" figures of merit
please see the
[`inspiral_range`](https://git.ligo.org/gwinc/inspiral-range) package.


## usage

### command line interface

`pygwinc` provides a command line interface that can be used to
calculate and plot noise budgets for the various canonical IFOs
described above, save/plot hdf5 trace data, and dump budget IFO
parameters:
```shell
$ python3 -m gwinc aLIGO
```

You can play with IFO parameters and see the effects on the budget by
dumping the pre-defined parameters to a [YAML-formatted parameter
file](#yaml-parameter-files), editing the parameter file, and
re-calculating the noise budget:
```shell
$ python3 -m gwinc --yaml aLIGO > my_aLIGO.yaml
$ edit my_aLIGO.yaml
$ python3 -m gwinc -d my_aLIGO.yaml aLIGO
                             aLIGO    my_aLIGO.yaml
Materials.Coating.Philown    5e-05            3e-05
$ python3 -m gwinc my_aLIGO.yaml
```

Stand-alone YAML files will always assume the nominal ['aLIGO' budget
description](gwinc/ifo/aLIGO).

[Custom budgets](#budget-interface) may also be processed by providing
the path to the budget module/package:
```shell
$ python3 -m gwinc path/to/mybudget
```

See command help for more info:
```shell
$ python3 -m gwinc -h
```


### basic python library usage

For custom code, parameter optimization, etc. all functionality can be
accessed through the `gwinc` library interface:
```python
>>> import gwinc
>>> import numpy as np
>>> freq = np.logspace(1, 3, 1000)
>>> Budget = gwinc.load_budget('aLIGO')
>>> traces = Budget(freq).run()
>>> fig = gwinc.plot_noise(freq, traces)
>>> fig.show()
```

The `load_budget()` function takes most of the same inputs as the
command line interface (e.g. IFO names, budget module paths, YAML
parameter files), and returns the un-instantiated Budget class defined
in the specified budget module.


## noise functions

`pygwinc` noise functions are available in the `gwinc.noise` package.
This package includes multiple modules for the different types of
noises, e.g. `suspensionthermal`, `coatingthermal`, `quantum`, etc.)

The various noise functions need many different parameters to
calculate their noise outputs.  The parameters are expected to be in
the form of object attributes, e.g.:
```python
def coating_brownian(f, materials, wavelength, wBeam, dOpt):
    ...
    # extract substructures
    sub = materials.Substrate
    ...
    # substrate properties
    Ysub = sub.MirrorY
```

The `materials` input argument in this case is expected to be an
object with a `Substrate` attribute, which itself has a `MirrorY`
attribute.


### `gwinc.Struct` objects

To make all this easier `pygwinc` provides a `Struct` class that can
hold parameters in attributes and additionally acts like a dictionary.
`Struct`s can be created from dictionaries, or loaded from various
file formats (see below).


### YAML parameter files

The easiest way to store all budget parameters is in a YAML file.
YAML files can be loaded directly into `gwinc.Struct` objects via
the `Struct.from_file()` class method:
```python
from gwinc import Struct
ifo = Struct.from_file('/path/to/ifo.yaml')
```

YAML parameter files can also be given to the `load_budget()` function
as described above, in which case the base 'aLIGO' budget structure
will be assumed and returned, with the YAML Struct inserted in the
`Budget.ifo` class attribute.

Here are the included ifo.yaml files for all the canonical IFOs:

* [aLIGO.yaml](gwinc/ifo/aLIGO/ifo.yaml)
* [Aplus.yaml](gwinc/ifo/Aplus/ifo.yaml)
* [Voyager.yaml](gwinc/ifo/Voyager/ifo.yaml)
* [CE1.yaml](gwinc/ifo/CE1/ifo.yaml)
* [CE2.yaml](gwinc/ifo/CE2/ifo.yaml)

The `Struct.from_file()` class method can also load MATLAB structs
defined in .mat files, for compatibility with
[matgwinc](https://git.ligo.org/gwinc/matgwinc), and MATLAB .m files,
although the later requires the use of the [python MATLAB
engine](https://www.mathworks.com/help/matlab/matlab-engine-for-python.html).


## budget interface

The basic structure of a `pygwinc` noise budget is a "budget module".
A budget modules is a standard python module (single `.py` file) or
package (directory containing `__inti__.py` file) which contains at
least one `nb.Budget` class definition named after the module name.

The `gwinc.nb` package provides various `BudgetItem` classes that can
inherited to define the various components of a budget:

* `nb.Noise`: BudgetItem describing a noise source
* `nb.Calibration`: BudgetItem describing a noise calibration
* `nb.Budget`: BudgetItem describing a group of noises

Here's an example of a budget module name `MyBudget`:
```shell
$ find MyBudget
MyBudget/
MyBudget/__init__.py
MyBudget/ifo.yaml
$
```

```python
# MyBudget/__init__.py

import numpy as np
from gwinc import nb
from gwinc import noise


class SuspensionThermal(nb.Noise):
    """Suspension thermal noise"""
    style = dict(
        label='Suspension Thermal',
        color='#ad900d',
        linestyle='--',
    )

    def calc(self):
        n = noise.suspensionthermal.suspension_thermal(
            self.freq, self.ifo.sus)
        return 2 * n


class MeasuredNoise(nb.Noise):
    style = dict(
        label='Measured Noise',
        color='#838209',
        linestyle='-',
    )

    def load(self):
        psd, freq = np.loadtxt('/path/to/measured/psd.txt')
        self.data = self.interpolate(freq, psd)

    def calc(self):
        return self.data


class MyCalibration(nb.Calibration):
    def calc(self):
        return np.ones_like(self.freq) * 1234


class MyBudget(nb.Budget):
    noises = [
        SuspensionThermal,
        MeasuredNoise,
    ]
    
    calibrations = [
        MyCalibration,
    ]
```

This budget would then be loaded with the `gwinc.load_budget()`
function, and processed with the `Budget.run()` method:
```python
Budget = load_budget('/path/to/MyBudget')
budget = Budget(freq)
traces = budget.run()
```

Other than the necessary `frequency` Budget init argument, any
additional keyword arguments are assigned as class attributes to the
budget object, and to all of it's constituent sub
noises/calibrations/budgets.

Note that the `SuspensionThermal` Noise class above uses the
`suspension_thermal` analytic noise calculation function, which takes
a "suspension" Struct as input argument.  In this case, this
suspension Struct is extracted from the `self.ifo` Struct at
`self.ifo.sus`.

If a budget module defined as a package includes an `ifo.yaml`
[parameter file](#parameter-files) in the package directory, the
`load_budget()` function will automatically load the YAML data into a
`gwinc.Struct` and include it as an `Budget.ifo` attribute in the
returned `Budget` class.  This would provide the `self.ifo` needed in
the `SuspensionThermal` Noise class above and is therefore a
convenient way to provide parameter structures in budget packages.
Otherwise it would need to be created/loaded in some other way and
passed to the budget at instantiation, e.g.:
```python
Budget = load_budget('/path/to/MyBudget')
ifo = Struct.from_file('/path/to/MyBudget.ifo')
budget = Budget(freq, ifo=ifo)
traces = budget.run()
```

The IFOs included in `gwinc.ifo` provide examples of the use of the
budget interface:

* [aLIGO](gwinc/ifo/aLIGO)
* [Aplus](gwinc/ifo/Aplus)
* [Voyager](gwinc/ifo/Voyager)
* [CE1](master/gwinc/ifo/CE1)
* [CE2](master/gwinc/ifo/CE2)


### BudgetItem load/update/calc methods

The Noise/Calibration/Budget `BudgetItem`s have three core methods
that can be overridden by the user to handle arbitrary data
processing.  This is useful for creating budgets from "live" dynamic
noise measurements and the like:

* `load()`: initial loading of static data
* `update(**kwargs)`: update data/attributes
* `calc()`: return final data array

See the built-in documentation for more info (e.g. `pydoc3
gwinc.nb.BudgetItem`)
  
  
<!-- ## comparison with MATLAB gwinc -->

<!-- `pygwinc` includes the ability use MATLAB gwinc directly via the -->
<!-- MATLAB python interface (see the CLI '--matlab' option above).  This -->
<!-- also allows for easy direct comparison between the pygwinc and -->
<!-- matgwinc noise budgets. -->

<!-- If you have a local checkout of matgwinc (at e.g. /path/to/gwinc) and -->
<!-- a local installation of MATLAB and it's python interface (at -->
<!-- e.g. /opt/matlab/python/lib/python3.6/site-packages) you can run the -->
<!-- comparison as so: -->
<!-- ```shell -->
<!-- $ export GWINCPATH=/path/to/matgwinc -->
<!-- $ export PYTHONPATH=/opt/matlab/python/lib/python3.6/site-packages -->
<!-- $ python3 -m gwinc.test -p aLIGO -->
<!-- ``` -->
<!-- This will produce a summary page of the various noise spectra that -->
<!-- differ between matgwinc and pygwinc. -->

<!-- Latest comparison plots from continuous integration: -->

<!-- * [aLIGO comparison](https://gwinc.docs.ligo.org/pygwinc/aLIGO_test.png) -->
<!-- * [A+ comparison](https://gwinc.docs.ligo.org/pygwinc/A+_test.png) -->
