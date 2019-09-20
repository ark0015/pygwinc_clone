import numpy as np

from .. import nb
from .. import noise

class QuantumVacuum(nb.Noise):
    """Quantum Vacuum

    """
    style = dict(
        label='Quantum Vacuum',
        color='#ad03de',
    )

    def calc(self):
        return noise.quantum.shotrad(self.freq, self.ifo)


class Seismic(nb.Noise):
    """Seismic

    """
    style = dict(
        label='Seismic',
        color='#855700',
    )

    def calc(self):
        if 'PlatformMotion' in self.ifo.Seismic:
            if self.ifo.Seismic.PlatformMotion == 'BSC':
                nt, nr = noise.seismic.seismic_BSC_ISI(self.freq)
            elif self.ifo.Seismic.PlatformMotion == '6D':
                nt, nr = noise.seismic.seismic_BSC_ISI_6D(self.freq)
            else:
                nt, nr = noise.seismic.seismic_BSC_ISI(self.freq)
        else:
            nt, nr = noise.seismic.seismic_BSC_ISI(self.freq)
        n, nh, nv = noise.seismic.seismic_suspension_fitered(
            self.ifo.Suspension, nt)
        return n * 4 * self.ifo.gwinc.dhdl_sqr


class Newtonian(nb.Noise):
    """Newtonian Gravity

    """
    style = dict(
        label='Newtonian Gravity',
        color='#15b01a',
    )

    def calc(self):
        return noise.newtonian.gravg(self.freq, self.ifo)


class NewtonianRayleigh(nb.Noise):
    """Newtonian Gravity, Rayleigh waves

    """
    style = dict(
        label='Newtonian (Rayleigh waves)',
        color='#1b2431',
    )

    def calc(self):
        return noise.newtonian.gravg_rayleigh(self.freq, self.ifo)


class NewtonianBody(nb.Noise):
    """Newtonian Gravity, body waves

    """
    style = dict(
        label='Newtonian (body waves)',
        color='#85a3b2',
    )

    def calc(self):
        return (noise.newtonian.gravg_pwave(self.freq, self.ifo)
               + noise.newtonian.gravg_swave(self.freq, self.ifo))


class NewtonianInfrasound(nb.Noise):
    """Newtonian Gravity, infrasound

    """
    style = dict(
        label='Newtonian (infrasound)',
        color='#ffa62b',
    )

    def calc(self):
        return noise.newtonian.atmois(self.freq, self.ifo)


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
        wavelength = self.ifo.Laser.Wavelength
        materials = self.ifo.Materials
        wBeam_ITM = self.ifo.Optics.ITM.BeamRadius
        wBeam_ETM = self.ifo.Optics.ETM.BeamRadius
        dOpt_ITM = self.ifo.Optics.ITM.CoatLayerOpticalThickness
        dOpt_ETM = self.ifo.Optics.ETM.CoatLayerOpticalThickness
        nITM = noise.coatingthermal.coating_brownian(
            self.freq, materials, wavelength, wBeam_ITM, dOpt_ITM)
        nETM = noise.coatingthermal.coating_brownian(
            self.freq, materials, wavelength, wBeam_ETM, dOpt_ETM)
        return (nITM + nETM) * 2 * self.ifo.gwinc.dhdl_sqr


class CoatingThermoOptic(nb.Noise):
    """Coating Thermo-Optic

    """
    style = dict(
        label='Coating Thermo-Optic',
        color='#02ccfe',
        linestyle='--',
    )

    def calc(self):
        wavelength = self.ifo.Laser.Wavelength
        materials = self.ifo.Materials
        wBeam_ITM = self.ifo.Optics.ITM.BeamRadius
        wBeam_ETM = self.ifo.Optics.ETM.BeamRadius
        dOpt_ITM = self.ifo.Optics.ITM.CoatLayerOpticalThickness
        dOpt_ETM = self.ifo.Optics.ETM.CoatLayerOpticalThickness
        nITM, junk1, junk2, junk3 = noise.coatingthermal.coating_thermooptic(
            self.freq, materials, wavelength, wBeam_ITM, dOpt_ITM[:])
        nETM, junk1, junk2, junk3 = noise.coatingthermal.coating_thermooptic(
            self.freq, materials, wavelength, wBeam_ETM, dOpt_ETM[:])
        return (nITM + nETM) * 2 * self.ifo.gwinc.dhdl_sqr


class ITMThermoRefractive(nb.Noise):
    """ITM Thermo-Refractive

    """
    style = dict(
        label='ITM Thermo-Refractive',
        color='#448ee4',
        linestyle='--',
    )

    def calc(self):
        gPhase = self.ifo.gwinc.finesse * 2/np.pi
        n = noise.substratethermal.substrate_thermorefractive(
            self.freq, self.ifo.Materials, self.ifo.Optics.ITM.BeamRadius)
        return n * 2 / (gPhase)**2 * self.ifo.gwinc.dhdl_sqr


class ITMCarrierDensity(nb.Noise):
    """ITM Carrier Density

    """
    style = dict(
        label='ITM Carrier Density',
        color='#929591',
        linestyle='--',
    )

    def calc(self):
        gPhase = self.ifo.gwinc.finesse * 2/np.pi
        n = noise.substratethermal.substrate_carrierdensity(
            self.freq, self.ifo.Materials, self.ifo.Optics.ITM.BeamRadius)
        return n * 2 / (gPhase)**2 * self.ifo.gwinc.dhdl_sqr


class SubstrateBrownian(nb.Noise):
    """Substrate Brownian

    """
    style = dict(
        label='Substrate Brownian',
        color='#fb7d07',
        linestyle='--',
    )

    def calc(self):
        nITM = noise.substratethermal.substrate_brownian(
            self.freq, self.ifo.Materials, self.ifo.Optics.ITM.BeamRadius)
        nETM = noise.substratethermal.substrate_brownian(
            self.freq, self.ifo.Materials, self.ifo.Optics.ETM.BeamRadius)
        return (nITM + nETM) * 2 * self.ifo.gwinc.dhdl_sqr


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
        color='#add00d',
        linestyle='--',
    )

    def calc(self):
        return noise.residualgas.gas(self.freq, self.ifo)
