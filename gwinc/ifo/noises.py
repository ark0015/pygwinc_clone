import numpy as np
from numpy import pi, sin, exp, sqrt

from .. import const
from .. import nb
from .. import noise

##################################################

def dhdl(f, armlen):
    """Strain to length conversion for noise power spetra

    This takes into account the GW wavelength and is only important
    when this is comparable to the detector arm length.

    From R. Schilling, CQG 14 (1997) 1513-1519, equation 5,
    with n = 1, nu = 0.05, ignoring overall phase and cos(nu)^2.
    A small value of nu is used instead of zero to avoid infinities.

    Returns the square of the dh/dL function, and the same divided by
    the arm length squared.

    """
    c = const.c
    nu_small = 15*pi/180
    omega_arm = pi * f * armlen / c
    omega_arm_f = (1 - sin(nu_small)) * pi * f * armlen / c
    omega_arm_b = (1 + sin(nu_small)) * pi * f * armlen / c
    sinc_sqr = 4 / abs(sin(omega_arm_f) * exp(-1j * omega_arm) / omega_arm_f
                       + sin(omega_arm_b) * exp(1j * omega_arm) / omega_arm_b)**2
    # keep DC value equal to 1
    sinc_sqr /= sinc_sqr[0]
    dhdl_sqr = sinc_sqr / armlen**2
    return dhdl_sqr, sinc_sqr

##################################################

class Strain(nb.Calibration):
    def calc(self):
        dhdl_sqr, sinc_sqr = dhdl(self.freq, self.ifo.Infrastructure.Length)
        return dhdl_sqr

##################################################

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
        return n * 4


class Newtonian(nb.Noise):
    """Newtonian Gravity

    """
    style = dict(
        label='Newtonian Gravity',
        color='#15b01a',
    )

    def calc(self):
        n = noise.newtonian.gravg(self.freq, self.ifo.Seismic)
        return n * 4


class NewtonianRayleigh(nb.Noise):
    """Newtonian Gravity, Rayleigh waves

    """
    style = dict(
        label='Newtonian (Rayleigh waves)',
        color='#1b2431',
    )

    def calc(self):
        n = noise.newtonian.gravg_rayleigh(self.freq, self.ifo.Seismic)
        return n * 2


class NewtonianBody(nb.Noise):
    """Newtonian Gravity, body waves

    """
    style = dict(
        label='Newtonian (body waves)',
        color='#85a3b2',
    )

    def calc(self):
        np = noise.newtonian.gravg_pwave(self.freq, self.ifo.Seismic)
        ns = noise.newtonian.gravg_swave(self.freq, self.ifo.Seismic)
        return (np + ns) * 4


class NewtonianInfrasound(nb.Noise):
    """Newtonian Gravity, infrasound

    """
    style = dict(
        label='Newtonian (infrasound)',
        color='#ffa62b',
    )

    def calc(self):
        n = noise.newtonian.atmois(self.freq, self.ifo.Atmospheric, self.ifo.Seismic)
        return n * 2


class SuspensionThermal(nb.Noise):
    """Suspension Thermal

    """
    style = dict(
        label='Suspension Thermal',
        color='#0d75f8',
    )

    def calc(self):
        n = noise.suspensionthermal.suspension_thermal(self.freq, self.ifo.Suspension)
        return n * 4


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
        return (nITM + nETM) * 2


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
        return (nITM + nETM) * 2


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
        return n * 2 / (gPhase)**2


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
        return n * 2 / (gPhase)**2


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
        return (nITM + nETM) * 2


class SubstrateThermoElastic(nb.Noise):
    """Substrate Thermo-Elastic

    """
    style = dict(
        label='Substrate Thermo-Elastic',
        color='#f5bf03',
        linestyle='--',
    )

    def calc(self):
        nITM = noise.substratethermal.substrate_thermoelastic(
            self.freq, self.ifo.Materials, self.ifo.Optics.ITM.BeamRadius)
        nETM = noise.substratethermal.substrate_thermoelastic(
            self.freq, self.ifo.Materials, self.ifo.Optics.ETM.BeamRadius)
        return (nITM + nETM) * 2


class ExcessGas(nb.Noise):
    """Excess Gas

    """
    style = dict(
        label='Excess Gas',
        color='#add00d',
        linestyle='--',
    )

    def calc(self):
        n = noise.residualgas.residual_gas_cavity(self.freq, self.ifo)
        # FIXME HACK: it's unclear if a phase noise in the arms like
        # the excess gas noise should get the same dhdL strain
        # calibration as the other displacement noises.  However, we
        # would like to use the one Strain calibration for all noises,
        # so we need to divide by the sinc_sqr here to undo the
        # application of the dhdl in the Strain calibration.  But this
        # is ultimately a superfluous extra calculation with the only
        # to provide some simplication in the Budget definition, so
        # should be re-evaluated at some point.
        dhdl_sqr, sinc_sqr = dhdl(self.freq, self.ifo.Infrastructure.Length)
        return n * 2 / sinc_sqr
