from gwinc.ifo.noises import *


class CE(nb.Budget):

    name = 'Cosmic Explorer'

    noises = [
        QuantumVacuum,
        Seismic,
        Newtonian,
        SuspensionThermal,
        CoatingBrownian,
        CoatingThermoOptic,
        ITMThermoRefractive,
        ITMCarrierDensity,
        SubstrateBrownian,
        SubstrateThermoElastic,
        ExcessGas,
    ]
