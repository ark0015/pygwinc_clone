'''Functions to calculate substrate thermal noise

'''

from __future__ import division, print_function
from numpy import exp, inf, pi, sqrt
import numpy as np
import scipy.special
import scipy.integrate

from .. import const
from ..const import BESSEL_ZEROS as zeta
from ..const import J0M as j0m


def substrate_carrierdensity(f, materials, wBeam, exact=False):
    """Substrate thermal displacement noise spectrum from charge carrier density fluctuations

    For semiconductor substrates.

    :f: frequency array in Hz
    :materials: gwinc optic materials structure
    :wBeam: beam radius (at 1 / e^2 power)
    :exact: whether to use adiabatic approximation or exact calculation (False)

    :returns: displacement noise power spectrum at :f:, in meters

    """
    H = materials.MassThickness
    diffElec = materials.Substrate.ElectronDiffusion
    diffHole = materials.Substrate.HoleDiffusion
    mElec = materials.Substrate.ElectronEffMass
    mHole = materials.Substrate.HoleEffMass
    cdDens = materials.Substrate.CarrierDensity
    gammaElec = materials.Substrate.ElectronIndexGamma
    gammaHole = materials.Substrate.HoleIndexGamma
    r0 = wBeam/np.sqrt(2)
    omega = 2*pi*f

    if exact:
        def integrand(k, om, D):
            return D * k**3 * exp(-k**2 * wBeam**2/4) / (D**2 * k**4 + om**2)
    
        integralElec = np.array([scipy.integrate.quad(lambda k: integrand(k, om, diffElec), 0, inf)[0] for om in omega])
        integralHole = np.array([scipy.integrate.quad(lambda k: integrand(k, om, diffHole), 0, inf)[0] for om in omega])

        # From P1400084 Heinert et al. Eq. 15
        #psdCD = @(gamma,m,int) 2*(3/pi^7)^(1/3)*kBT*H*gamma^2*m/hbar^2*cdDens^(1/3)*int; %units are meters
        # FIXME: why the unused argument here?
        def psdCD(gamma, m, int_):
            return 2/pi * H * gamma**2 * cdDens * int_

        psdElec = psdCD(gammaElec, mElec, integralElec)
        psdHole = psdCD(gammaHole, mHole, integralHole)
    else:
        psdElec = 4*H*gammaElec**2*cdDens*diffElec/(pi*r0**4*omega**2)
        psdHole = 4*H*gammaHole**2*cdDens*diffHole/(pi*r0**4*omega**2)

    return psdElec + psdHole


def substrate_thermorefractive(f, materials, wBeam, exact=False):
    """Substrate thermal displacement noise spectrum from thermorefractive fluctuations

    For semiconductor substrates.

    :f: frequency array in Hz
    :materials: gwinc optic materials structure
    :wBeam: beam radius (at 1 / e^2 power)
    :exact: whether to use adiabatic approximation or exact calculation (False)

    :returns: displacement noise power spectrum at :f:, in meters

    """
    H = materials.MassThickness
    kBT = const.kB * materials.Substrate.Temp
    Temp = materials.Substrate.Temp
    rho = materials.Substrate.MassDensity
    beta = materials.Substrate.dndT
    C = materials.Substrate.MassCM
    kappa = materials.Substrate.MassKappa
    r0 = wBeam/np.sqrt(2)
    omega = 2*pi*f

    if exact:
        def integrand(k, om, D):
            return D * k**3 * exp(-k**2 * wBeam**2/4) / (D**2 * k**4 + om**2)

        inte = np.array([scipy.integrate.quad(lambda k: integrand(k, om, kappa/(rho*C)), 0, inf)[0] for om in omega])

        # From P1400084 Heinert et al. Eq. 15
        #psdCD = @(gamma,m,int) 2*(3/pi^7)^(1/3)*kBT*H*gamma^2*m/hbar^2*cdDens^(1/3)*int; %units are meters
        psdTR = lambda int_: 2/pi * H * beta**2 * kBT * Temp / (rho*C) * int_;

        psd = psdTR(inte)
        psd = 2/pi * H * beta**2 * kBT * Temp / (rho*C) * inte

    else:
        psd = 4*H*beta**2*kappa*kBT*Temp/(pi*r0**4*omega**2*(rho*C)**2)

    return psd


def substrate_brownian(f, materials, wBeam):
    """Substrate thermal displacement noise spectrum due to substrate mechanical loss

    :f: frequency array in Hz
    :materials: gwinc optic materials structure
    :wBeam: beam radius (at 1 / e^2 power)

    :returns: displacement noise power spectrum at :f:, in meters

    """
    Y = materials.Substrate.MirrorY
    sigma = materials.Substrate.MirrorSigma
    c2 = materials.Substrate.c2
    n = materials.Substrate.MechanicalLossExponent
    alphas = materials.Substrate.Alphas
    kBT = const.kB * materials.Substrate.Temp

    cftm, aftm = substrate_brownian_FiniteCorr(materials, wBeam)

    # Bulk substrate contribution
    phibulk = c2 * f**n
    cbulk = 8 * kBT * aftm * phibulk / (2 * pi * f)

    # Surface loss contribution
    # csurf = alphas/(Y*pi*wBeam^2)
    csurf = alphas*(1-2*sigma)/((1-sigma)*Y*pi*wBeam**2)
    csurf *= 8 * kBT / (2 * pi * f)

    return csurf + cbulk


def substrate_brownian_FiniteCorr(materials, wBeam):
    """Substrate brownian noise finite-size test mass correction

    :materials: gwinc optic materials structure
    :wBeam: beam radius (at 1 / e^2 power)

    :returns: correction factors tuple:
    cftm = finite mirror correction factor
    aftm = amplitude coefficient for thermal noise:
           thermal noise contribution to displacement noise is
           S_x(f) = (8 * kB * T / (2*pi*f)) * Phi(f) * aftm

    Equation references to Bondu, et al. Physics Letters A 246 (1998)
    227-236 (hereafter BHV) and Liu and Thorne gr-qc/0002055 (hereafter LT)

    """
    a = materials.MassRadius
    h = materials.MassThickness
    Y = materials.Substrate.MirrorY
    sigma = materials.Substrate.MirrorSigma

    # LT uses e-folding of power
    r0 = wBeam / sqrt(2)
    km = zeta/a

    Qm = exp(-2*km*h) # LT eq. 35a

    Um = (1-Qm)*(1+Qm)+4*h*km*Qm
    Um = Um/((1-Qm)**2-4*(km*h)**2*Qm) # LT 53 (BHV eq. btwn 29 & 30)

    x = exp(-(zeta*r0/a)**2/4)
    s = sum(x/(zeta**2*j0m)) # LT 57

    x2 = x*x
    U0 = sum(Um*x2/(zeta*j0m**2))
    U0 = U0*(1-sigma)*(1+sigma)/(pi*a*Y) # LT 56 (BHV eq. 3)

    p0 = 1/(pi*a**2) # LT 28
    DeltaU = (pi*h**2*p0)**2
    DeltaU = DeltaU + 12*pi*h**2*p0*sigma*s
    DeltaU = DeltaU + 72*(1-sigma)*s**2
    DeltaU = DeltaU*a**2/(6*pi*h**3*Y) # LT 54

    # LT 58 (eq. following BHV 31)
    aftm = DeltaU + U0

    # amplitude coef for infinite TM, LT 59
    # factored out: (8 * kB * T * Phi) / (2 * pi * f)
    aitm = (1 - sigma**2) / (2 * sqrt(2 * pi) * Y * r0)

    # finite mirror correction
    cftm = aftm / aitm

    return cftm, aftm


def subtherm(f, ifo):
    """Noise from thermoelastic fluctuations in mirror

    """
    wITM = ifo.Optics.ITM.BeamRadius
    wETM = ifo.Optics.ETM.BeamRadius
    sigma = ifo.Materials.Substrate.MirrorSigma

    L = ifo.Infrastructure.Length
    kBT = const.kB * ifo.Materials.Substrate.Temp

    rho = ifo.Materials.Substrate.MassDensity
    kappa = ifo.Materials.Substrate.MassKappa # thermal conductivity
    alpha = ifo.Materials.Substrate.MassAlpha # thermal expansion
    CM = ifo.Materials.Substrate.MassCM # heat capacity @ constant mass
    Temp = ifo.Materials.Substrate.Temp # temperature

    S0 = 8*(1+sigma)**2*kappa*alpha**2*Temp*kBT # note kBT has factor Temp
    S0 = S0/(sqrt(2*pi)*(CM*rho)**2)
    SITM = S0/(wITM/sqrt(2))**3 # LT 18 less factor 1/omega^2
    SETM = S0/(wETM/sqrt(2))**3 # LT 18 less factor 1/omega^2

    # Corrections for finite test masses:
    SITM = SITM * subthermFiniteCorr(ifo, 'ITM')
    SETM = SETM * subthermFiniteCorr(ifo, 'ETM')

    # 2 mirrors each type, factor omega^2, dimensionless strain
    n = 2 * (SITM + SETM)/(2*pi*f)**2 * ifo.gwinc.dhdl_sqr

    return n


def subthermFiniteCorr(ifo, opticName):
    """Finite size test mass correction to noise amplitude coefficient

    (Liu & Thorne gr-qc/0002055 equation 46)
    
    Equation references to Bondu, et al. Physics Letters A 246 (1998)
    227-236 (hereafter BHV) or Liu and Thorne gr-qc/0002055 (hereafter LT)

    """
    # extract some numbers
    a = ifo.Materials.MassRadius
    h = ifo.Materials.MassThickness
    w = ifo.Optics[opticName].BeamRadius
    sigma = ifo.Materials.Substrate.MirrorSigma

    # do the work
    r0 = w/sqrt(2) # LT uses power e-folding
    km = zeta/a

    Qm = exp(-2*km*h) # LT eq. 35a

    pm = exp(-(km*r0)**2/4)/(pi*(a*j0m)**2) # LT 37

    c0 = 6*(a/h)**2*sum(j0m*pm/zeta**2) # LT 32
    c1 = -2*c0/h # LT 32
    p0 = 1/(pi*a**2) # LT 28
    c1 = c1 + p0/(2*h) # LT 40

    coeff = (1-Qm)*((1-Qm)*(1+Qm)+8*h*km*Qm)
    coeff = coeff + 4*(h*km)**2*Qm*(1+Qm)
    coeff = coeff*km*(pm*j0m)**2*(1-Qm)
    coeff = coeff/((1-Qm)**2-4*(h*km)**2*Qm)**2
    coeff = sum(coeff) + h*c1**2/(1+sigma)**2
    coeff = coeff*(sqrt(2*pi)*r0)**3*a**2 # LT 46
    return coeff
