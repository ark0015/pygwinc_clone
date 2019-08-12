from gwinc.ifo.noises import *


class CE2(nb.Budget):

    name = 'Cosmic Explorer 2'

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
