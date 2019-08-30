from gwinc.ifo.noises import *


class CE(nb.Budget):

    name = 'Cosmic Explorer'

    noises = [
        QuantumVacuum,
        Seismic,
        Newtonian,
        AtmosphericInfrasound,
        SuspensionThermal,
        CoatingBrownian,
        CoatingThermoOptic,
        ITMThermoRefractive,
        ITMCarrierDensity,
        SubstrateBrownian,
        SubstrateThermoElastic,
        ExcessGas,
    ]
