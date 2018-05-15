from __future__ import division
from numpy import pi, sqrt
from scipy.io import loadmat
import scipy.special
import logging

from .struct import Struct
from .noise.coatingthermal import getCoatDopt


def precompIFO(ifo, PRfixed=True):
    """Add precomputed data to the IFO model.

    To prevent recomputation of these precomputed data, if the
    ifo argument contains ifo.gwinc.PRfixed, and this matches
    the argument PRfixed, no changes are made.

    """

    if 'gwinc' not in ifo:
        ifo.gwinc = Struct()

    ifo.gwinc.PRfixed = PRfixed

    ##############################
    # derived temp

    if 'Temp' not in ifo.Materials.Substrate:
        ifo.Materials.Substrate.Temp = ifo.Constants.Temp

    ##############################
    # suspension vertical-horizontal coupling

    if 'VHCoupling' not in ifo.Suspension:
        ifo.Suspension.VHCoupling = Struct()
        ifo.Suspension.VHCoupling.theta = ifo.Infrastructure.Length / ifo.Constants.R_earth

    ##############################
    # optics values

    # calculate optics' parameters
    ifo.Materials.MirrorVolume = pi*ifo.Materials.MassRadius**2 * \
                                 ifo.Materials.MassThickness
    ifo.Materials.MirrorMass = ifo.Materials.MirrorVolume* \
                               ifo.Materials.Substrate.MassDensity
    ifo.Optics.ITM.Thickness = ifo.Materials.MassThickness

    # coating layer optical thicknesses - mevans 2 May 2008
    if 'CoatLayerOpticalThickness' not in ifo.Optics.ITM:
        T = ifo.Optics.ITM.Transmittance
        dL = ifo.Optics.ITM.CoatingThicknessLown
        dCap = ifo.Optics.ITM.CoatingThicknessCap
        ifo.Optics.ITM.CoatLayerOpticalThickness = getCoatDopt(ifo, T, dL, dCap=dCap)
        T = ifo.Optics.ETM.Transmittance
        dL = ifo.Optics.ETM.CoatingThicknessLown
        dCap = ifo.Optics.ETM.CoatingThicknessCap
        ifo.Optics.ETM.CoatLayerOpticalThickness = getCoatDopt(ifo, T, dL, dCap=dCap)

    ##############################
    # beam parameters

    armlen = ifo.Infrastructure.Length

    # g-factors
    g1 = 1 - armlen / ifo.Optics.Curvature.ITM
    g2 = 1 - armlen / ifo.Optics.Curvature.ETM
    gcav = sqrt(g1 * g2 * (1 - g1 * g2))
    gden = g1 - 2 * g1 * g2 + g2

    if (g1 * g2 * (1 - g1 * g2)) <= 0:
        raise Exception('Unstable arm cavity g-factors.  Change ifo.Optics.Curvature')
    elif gcav < 1e-3:
        log.warning('Nearly unstable arm cavity g-factors.  Reconsider ifo.Optics.Curvature')

    ws = sqrt(armlen * ifo.Laser.Wavelength / pi)
    w1 = ws * sqrt(abs(g2) / gcav)
    w2 = ws * sqrt(abs(g1) / gcav)

    # waist size
    w0 = ws * sqrt(gcav / abs(gden))
    zr = pi * w0**2 / ifo.Laser.Wavelength
    z1 = armlen * g2 * (1 - g1) / gden
    z2 = armlen * g1 * (1 - g2) / gden

    ifo.Optics.ITM.BeamRadius = w1
    ifo.Optics.ETM.BeamRadius = w2

    ##############################
    # calc power and IFO parameters

    pbs, parm, finesse, prfactor, Tpr = precompPower(ifo, PRfixed)

    ifo.gwinc.pbs = pbs
    ifo.gwinc.parm = parm
    ifo.gwinc.finesse = finesse
    ifo.gwinc.prfactor = prfactor
    ifo.gwinc.gITM = g1
    ifo.gwinc.gETM = g2
    ifo.gwinc.BeamWaist = w0
    ifo.gwinc.BeamRayleighRange = zr
    ifo.gwinc.BeamWaistToITM = z1
    ifo.gwinc.BeamWaistToETM = z2
    ifo.Optics.PRM.Transmittance = Tpr

    ##############################
    # calc quantum parameters

    fSQL, fGammaIFO, fGammaArm = precompQuantum(ifo)

    ifo.gwinc.fSQL = fSQL
    ifo.gwinc.fGammaIFO = fGammaIFO
    ifo.gwinc.fGammaArm = fGammaArm

    ##############################
    # precompute bessels zeros, needed in coat and substrate thermal

    besselzeros = scipy.special.jn_zeros(1, 300)
    ifo.Constants.BesselZeros = besselzeros

    ##############################
    # saved seismic spectrum

    if 'darmSeiSusFile' in ifo.Seismic and ifo.Seismic.darmSeiSusFile:
        darmsei = loadmat(ifo.Seismic.darmSeiSusFile)
        ifo.Seismic.darmseis_f = darmsei['darmseis_f'][0]
        ifo.Seismic.darmseis_x = darmsei['darmseis_x'][0]

    return ifo


def precompPower(ifo, PRfixed=True):
    """Compute power on beamsplitter, finesse, and power recycling factor.

    """
    # constants
    c       = scipy.constants.c
    pin     = ifo.Laser.Power
    lambda_ = ifo.Laser.Wavelength
    t1      = sqrt(ifo.Optics.ITM.Transmittance)
    r1      = sqrt(1 - ifo.Optics.ITM.Transmittance)
    t2      = sqrt(ifo.Optics.ETM.Transmittance)
    r2      = sqrt(1 - ifo.Optics.ETM.Transmittance)
    t3      = sqrt(ifo.Optics.SRM.Transmittance)
    t5      = sqrt(ifo.Optics.PRM.Transmittance)
    r5      = sqrt(1 - ifo.Optics.PRM.Transmittance)
    wl      = 2*pi * c/lambda_
    lrec    = ifo.Optics.SRM.CavityLength
    effic   = ifo.Optics.PhotoDetectorEfficiency
    loss    = ifo.Optics.Loss                          # single TM loss
    bsloss  = ifo.Optics.BSLoss
    acoat   = ifo.Optics.ITM.CoatingAbsorption
    pcrit   = ifo.Optics.pcrit

    # Finesse, effective number of bounces in cavity, power recycling factor
    finesse = 2*pi / (t1**2 + 2*loss)        # arm cavity finesse
    neff    = 2 * finesse / pi

    # Arm cavity reflectivity with finite loss
    garm = t1 / (1 - r1*r2*sqrt(1-2*loss))  # amplitude gain wrt input field
    rarm = r1 - t1 * r2 * sqrt(1-2*loss) * garm

    if PRfixed:
        Tpr = ifo.Optics.PRM.Transmittance  # use given value
    else:
        #prfactor = 1/(2*loss * neff + bsloss);         % power recycling factor
        Tpr = 1-(rarm*sqrt(1-bsloss))**2 # optimal recycling mirror transmission
        t5 = sqrt(Tpr)
        r5 = sqrt(1 - Tpr)
    prfactor = t5**2 / (1 + r5 * rarm * sqrt(1-bsloss))**2

    pbs  = pin * prfactor          # BS power from input power
    parm = pbs * garm**2 / 2       # arm power from BS power

    asub = 1.3*2*ifo.Optics.ITM.Thickness*ifo.Optics.SubstrateAbsorption
    pbsl = 2*pcrit/(asub+acoat*neff) # bs power limited by thermal lensing

    #pbs = min([pbs, pbsl]);
    if pbs > pbsl:
        logging.warning('P_BS exceeds BS Thermal limit!')

    return pbs, parm, finesse, prfactor, Tpr


def precompQuantum(ifo):
    """Compute quantum noise parameters.

    """

    # physical constants
    hbar = scipy.constants.hbar # J s
    c    = scipy.constants.c    # m / s

    # IFO parameters
    lambda_= ifo.Laser.Wavelength
    Titm   = ifo.Optics.ITM.Transmittance
    Tsrm   = ifo.Optics.SRM.Transmittance
    m      = ifo.Materials.MirrorMass
    L      = ifo.Infrastructure.Length
    Lsrc   = ifo.Optics.SRM.CavityLength

    # power on BS (W) computed in precompBSPower
    Pbs    = ifo.gwinc.pbs
    Parm   = ifo.gwinc.parm

    # derived parameters
    w0 = 2 * pi * c / lambda_      # carrier frequency (rad/s)
    gammaArm = Titm * c / (4 * L)  # arm cavity pole (rad/s)
    fGammaArm = gammaArm / (2*pi)
    rSR = sqrt(1 - Tsrm)

    # fSQL as defined in D&D paper (eq 33 in P1400018 and/or PRD paper)
    tSR = sqrt(Tsrm)
    fSQL = (1/(2*pi))*(8/c)*sqrt((Parm*w0)/(m*Titm))*(tSR/(1+rSR))

    # gammaIFO in Hz
    fGammaIFO = fGammaArm * ((1 + rSR) / (1 - rSR))

    return fSQL, fGammaIFO, fGammaArm
