from .. import nb
from .. import noise


class QuantumVacuum(nb.Noise):
    """Quantum Vacuum

    """
    style = dict(
        label='Quantum',
        color='#ad03de',
    )

    def calc(self):
        return noise.quantum.shotrad(self.freq, self.ifo)


class Seismic(nb.Noise):
    """Seismic

    """
    style = dict(
        label='Seismic',
        color='#653700',
    )

    def calc(self):
        return noise.seismic.seismic(self.freq, self.ifo)


class Newtonian(nb.Noise):
    """Newtonian Gravity

    """
    style = dict(
        label='Newtonian Gravity',
        color='#15b01a',
    )

    def calc(self):
        return noise.newtonian.gravg(self.freq, self.ifo)


class SuspensionThermal(nb.Noise):
    """Suspension Thermal

    """
    style = dict(
        label='Suspension Thermal',
        color='#0d75f8',
    )

    def calc(self):
        return noise.suspensionthermal.susptherm(self.freq, self.ifo)


class CoatingBrownian(nb.Noise):
    """Coating Brownian

    """
    style = dict(
        label='Coating Brownian',
        color='#fe0002',
    )

    def calc(self):
        return noise.coatingthermal.coatbrownian(self.freq, self.ifo)


class CoatingThermoOptic(nb.Noise):
    """Coating Thermo-Optic

    """
    style = dict(
        label='Coating Thermo-Optic',
        color='#02ccfe',
        linestyle='--',
    )

    def calc(self):
        return noise.coatingthermal.thermooptic(self.freq, self.ifo)


class ITMThermoRefractive(nb.Noise):
    """ITM Thermo-Refractive

    """
    style = dict(
        label='ITM Thermo-Refractive',
        color='#448ee4',
        linestyle='--',
    )

    def calc(self):
        return noise.substratethermal.thermorefractiveITM(self.freq, self.ifo)


class ITMCarrierDensity(nb.Noise):
    """ITM Carrier Density

    """
    style = dict(
        label='ITM Carrier Density',
        color='#929591',
        linestyle='--',
    )

    def calc(self):
        return noise.substratethermal.carrierdensity(self.freq, self.ifo)


class SubstrateBrownian(nb.Noise):
    """Substrate Brownian

    """
    style = dict(
        label='Substrate Brownian',
        color='#fb7d07',
        linestyle='--',
    )

    def calc(self):
        return noise.substratethermal.subbrownian(self.freq, self.ifo)


class SubstrateThermoElastic(nb.Noise):
    """Substrate Thermo-Elastic

    """
    style = dict(
        label='Substrate Thermo-Elastic',
        color='#f5bf03',
        linestyle='--',
    )

    def calc(self):
        return noise.substratethermal.subtherm(self.freq, self.ifo)


class ExcessGas(nb.Noise):
    """Excess Gas

    """
    style = dict(
        label='Excess Gas',
        color='#ad900d',
        linestyle='--',
    )

    def calc(self):
        return noise.residualgas.gas(self.freq, self.ifo)
