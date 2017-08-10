from __future__ import division, print_function
import scipy.constants
import scipy.special
import numpy as np
from numpy import pi, sum, zeros, exp, imag, sqrt, sin, cos, sinh, cosh


def coatbrownian(f, ifo):
    """effect of optical coating on Brownian thermal noise
    returns strain noise power spectrum in 1 / Hz
    
    Added by G Harry 8/3/02 from work by Nakagawa, Gretarsson, et al.
    Expanded to reduce approximation, GMH 8/03
    Modified to return strain noise, PF 4/07
    Modified to accept coating with non-quater-wave layers, mevans 25 Apr 2008"""

    Length = ifo.Infrastructure.Length
    wBeam_ITM = ifo.Optics.ITM.BeamRadius
    wBeam_ETM = ifo.Optics.ETM.BeamRadius
    dOpt_ITM = ifo.Optics.ITM.CoatLayerOpticalThickness
    dOpt_ETM = ifo.Optics.ETM.CoatLayerOpticalThickness

    # compute Brownian noise for specified coating structure
    SbrITM  = getCoatBrownian(f, ifo, wBeam_ITM, dOpt_ITM)
    SbrETM  = getCoatBrownian(f, ifo, wBeam_ETM, dOpt_ETM)

    n = 2 * (SbrITM + SbrETM) / Length**2

    return n


def getCoatBrownian(f, ifo, wBeam, dOpt):
    """returns the coating brownian noise for a given collection of
    coating layers.  The layers are assumed to be alernating low-n
    high-n layers, with low-n first.
    
    f = frequency vector in Hz
    ifo = parameter struct from IFOmodel.m
    
    opticName = name of the Optic struct to use for wBeam and dOpt
    wBeam = ifoArg.Optics.(opticName).BeamRadius
    dOpt = ifoArg.Optics.(opticName).CoatLayerOpticalThickness
    
    wBeam = beam radius (at 1 / e^2 power)
    dOpt = coating layer thickness vector (Nlayer x 1)
         = the optical thickness, normalized by lambda, of each coating layer.
    
    SbrZ = Brownian noise spectra for one mirror in m^2 / Hz
    
    adapted from bench62
    based on Harry et al., Class Quant Grav 24 (2007) 405-415"""

    # Constants
    subTemp = ifo.Materials.Substrate.Temp
    kBT    = scipy.constants.k * subTemp
    lambda_ = ifo.Laser.Wavelength

    Ysub     = ifo.Materials.Substrate.MirrorY
    sigmasub = ifo.Materials.Substrate.MirrorSigma
    
    Yhighn     = ifo.Materials.Coating.Yhighn
    sigmahighn = ifo.Materials.Coating.Sigmahighn
    phihighn   = ifo.Materials.Coating.Phihighn
    nH         = ifo.Materials.Coating.Indexhighn
    
    Ylown     = ifo.Materials.Coating.Ylown
    sigmalown = ifo.Materials.Coating.Sigmalown
    philown   = ifo.Materials.Coating.Philown
    nL        = ifo.Materials.Coating.Indexlown

    # compute thickness of each material in the coating
    dlown  = sum(dOpt[::2])  * lambda_ / nL
    dhighn = sum(dOpt[1::2]) * lambda_ / nH
    dCoat  = dlown + dhighn
    
    # for debugging, this is a rough but direct estimate
    Llown = dlown * philown
    Lhighn = dhighn * phihighn
    Lsum = Llown + Lhighn
    Lall = [Lsum, Llown, Lhighn]
    
    ################# this part is directly from bench62 #################
    Yperp = dCoat/(dhighn/Yhighn+dlown/Ylown)
    phiperp = Yperp/dCoat*(dlown*philown/Ylown + dhighn*phihighn/Yhighn)
    Ypara = 1/dCoat*(Yhighn*dhighn + Ylown*dlown)
    phipara = 1/(dCoat*Ypara)*(Ylown*philown*dlown + Yhighn*phihighn*dhighn)
    
    # This is a kludge, the real formula is very complicted but this
    # average works really well
    sigma1 = 1/2*(sigmahighn+sigmalown)

    # This is exact
    sigma2 = (sigmahighn*Yhighn*dhighn+sigmalown*Ylown*dlown)/(Yhighn*dhighn+Ylown*dlown)

    # Brownian contribution to coating thermal noise, low Poisson ratio limit
    # cITM = dITM/(pi*wITM^2)*(Ypara/Ysub^2*phipara+phiperp/Yperp);
    # cETM = dCoat/(pi*wETM^2)*(Ypara/Ysub^2*phipara+phiperp/Yperp);

    # Brownian contribution to coating thermal noise, full formula
    cbnoise = dCoat * (1-sigmasub**2)/(pi*wBeam**2) * (
        (1/(Yperp*(1-sigmasub**2)) - 2*sigma2**2*Ypara/(Yperp**2*(1-sigmasub**2)*(1-sigma1)))*phiperp
        + Ypara*sigma2*(1-2*sigmasub)/(Yperp*Ysub*(1-sigma1)*(1-sigmasub))*(phipara-phiperp)
        + Ypara*(1+sigmasub)*(1-2*sigmasub)**2/(Ysub**2*(1-sigma1**2)*(1-sigmasub))*phipara)

    # noise power spectrum
    SbrZ = 4 * kBT * cbnoise / (2 * pi * f)

    return SbrZ


def thermooptic(f, ifo):
    """effect of thermoelastic and thermorefractive fluctuations in the mirror dielectric coating
    returns strain noise power spectrum in 1 / Hz
    
    Added by G Harry 8/27/03 from work by Fejer, Rowan, Braginsky, et al
    thermoelastic and thermorefractive effects combined coherently in 2006
    
    Reduced approximations and combined TE and TR
    effects with correct sign, mevans 25 Apr 2008"""
  
    Length = ifo.Infrastructure.Length
    wBeam_ITM = ifo.Optics.ITM.BeamRadius
    wBeam_ETM = ifo.Optics.ETM.BeamRadius
    dOpt_ITM = ifo.Optics.ITM.CoatLayerOpticalThickness
    dOpt_ETM = ifo.Optics.ETM.CoatLayerOpticalThickness

    # compute ThermoOptic noise for specified coating structure
    StoITM, junk1, junk2, junk3  = getCoatThermoOptic(f, ifo, wBeam_ITM, dOpt_ITM[:])
    StoETM, junk1, junk2, junk3  = getCoatThermoOptic(f, ifo, wBeam_ETM, dOpt_ETM[:])
  
    n = 2 * (StoITM + StoETM) / Length**2
    return n


def getCoatThermoOptic(f, ifo, wBeam, dOpt):
    """return power spectra of coating thermo-optic noise for a single optic
    
    f = frequency vector in Hz
    ifo = parameter struct from IFOmodel.m
    opticName = name of the Optic struct to use for wBeam and dOpt
    wBeam = ifoArg.Optics.(opticName).BeamRadius
    dOpt = ifoArg.Optics.(opticName).CoatLayerOpticalThickness
    
    wBeam = beam radius (at 1 / e^2 power)
    dOpt   = optical thickness / lambda of each layer
           = geometrical thickness * refractive index / lambda
    
    StoZ = power spectra of apparent mirror position in m^2 / Hz
    SteZ = thermo-elastic componenet of StoZ
    StrZ = thermo-refractive componenet of StoZ
    T = coating power transmission, made available as a cross-check"""

    # compute coefficients
    dTO, dTR, dTE, T, junk = getCoatTOPos(ifo, wBeam, dOpt)

    # compute correction factors
    gTO = getCoatThickCorr(f, ifo, dOpt, dTE, dTR)
    gTE = getCoatThickCorr(f, ifo, dOpt, dTE, 0)
    gTR = getCoatThickCorr(f, ifo, dOpt, 0, dTR)

    # compute thermal source spectrum
    SsurfT, junk = getCoatThermal(f, ifo, wBeam)

    StoZ = SsurfT * gTO * dTO**2
    SteZ = SsurfT * gTE * dTE**2
    StrZ = SsurfT * gTR * dTR**2
    return (StoZ, SteZ, StrZ, T)


def getCoatTOPos(ifo, wBeam, dOpt):
    """returns mirror position derivative wrt thermal fluctuations
    
    ifo  = parameter struct from IFOmodel.m
    opticName = name of the Optic struct to use for wBeam and dOpt
    wBeam = ifoArg.Optics.(opticName).BeamRadius
    dOpt = ifoArg.Optics.(opticName).CoatLayerOpticalThickness
    
    wBeam = beam radius, for finite mirror correction (at 1 / e^2 power)
    dOpt = optical thickness / lambda of each layer
         = geometrical thickness * refractive index / lambda

    dTO = total thermo-optic dz/dT
    dTR = thermo-refractive dz/dT
    dTE = thermo-elastic dz/dT
    
    compute thermal fluctuations with getCoatThermal
    (see also T080101)"""
  
    # parameters
    lambda_ = ifo.Laser.Wavelength
    nS = ifo.Materials.Substrate.RefractiveIndex
  
    # compute refractive index, effective alpha and beta
    nLayer, aLayer, bLayer, dLayer, sLayer = getCoatLayers(ifo, dOpt)

    # compute coating average parameters
    dc, Cc, Kc, aSub = getCoatAvg(ifo, dOpt)
  
    # compute reflectivity and parameters
    dphi_dT, dphi_TE, dphi_TR, rCoat = getCoatTOPhase(1, nS, nLayer, dOpt, aLayer, bLayer, sLayer)
    R = abs(rCoat)**2
    T = 1 - R
  
    # for debugging
    #disp(sprintf('R = %.3f, T = %.0f ppm', R, 1e6 * T))
  
    # convert from phase to meters, subtracting substrate
    dTR = dphi_TR * lambda_ / (4 * pi)
    dTE = dphi_TE * lambda_ / (4 * pi) - aSub * dc

    # mirror finite size correction
    Cfsm = getCoatFiniteCorr(ifo, wBeam, dOpt)
    dTE = dTE * Cfsm
  
    # add TE and TR effects (sign is already included)
    dTO = dTE + dTR
  
    return dTO, dTR, dTE, T, R


def getCoatThickCorr(f, ifo, dOpt, dTE, dTR):
    """finite coating thickness correction
    Uses correction factor from T080101, "Thick Coating Correction" (Evans)
    
    (see getCoatThermoOptic for example usage)"""

    ##############################################
    # For comparison in the bTR = 0 limit, the
    # equation from Fejer (PRD70, 2004)
    # modified so that gFC -> 1 as xi -> 0
    #  gTC = (2 ./ (R * xi.^2)) .* (sh - s + R .* (ch - c)) ./ ...
    #    (ch + c + 2 * R * sh + R^2 * (ch - c));
    # which in the limit of xi << 1 becomes
    #  gTC = 1 - xi * (R - 1 / (3 * R));

    # parameter extraction
    pS = ifo.Materials.Substrate
    Cs = pS.MassCM * pS.MassDensity
    Ks = pS.MassKappa
  
    # compute coating average parameters
    dc, Cc, Kc, junk = getCoatAvg(ifo, dOpt)
  
    # R and xi (from T080101, Thick Coating Correction)
    w = 2 * pi * f
    R = sqrt(Cc * Kc / (Cs * Ks))
    xi = dc * sqrt(2 * w * Cc / Kc)
  
    # trig functions of xi
    s = sin(xi)
    c = cos(xi)
    sh = sinh(xi)
    ch = cosh(xi)
  
    # pR and pE (dTR = -\bar{\beta} lambda, dTE = \Delta \bar{\alpha} d)
    pR = dTR / (dTR + dTE)
    pE = dTE / (dTR + dTE)
  
    # various parts of gTC
    g0 = 2 * (sh - s) + 2 * R * (ch - c)
    g1 = 8 * sin(xi / 2) * (R * cosh(xi / 2) + sinh(xi / 2))
    g2 = (1 + R**2) * sh + (1 - R**2) * s + 2 * R * ch
    gD = (1 + R**2) * ch + (1 - R**2) * c + 2 * R * sh

    # and finally, the correction factor
    gTC = (pE**2 * g0 + pE * pR * xi * g1 + pR**2 * xi**2 * g2) / (R * xi**2 * gD)
    return gTC


def getCoatThermal(f, ifo, wBeam):
    """returns the thermal noise spectra for a surface layer
    
    f = frequency vector in Hz
    ifo = parameter struct from IFOmodel.m
    
    wBeam = beam radius (at 1 / e^2 power)
          = usually ifo.Optics.ITM.BeamRadius or ifo.Optics.ETM.BeamRadius
    
    SsurfT = power spectra of thermal fluctuations in K^2 / Hz
    rdel = thermal diffusion length at each frequency in m"""
  
    # use substrate temperature
    subTemp = ifo.Materials.Substrate.Temp

    kBT2 = scipy.constants.k * subTemp**2
  
    pS = ifo.Materials.Substrate
    C_S = pS.MassCM * pS.MassDensity
    K_S = pS.MassKappa

    # omega
    w = 2 * pi * f
  
    # thermal diffusion length
    rdel = sqrt(2 * K_S / (C_S * w))
  
    # noise equation
    SsurfT = 4 * kBT2 / (pi * w * C_S * rdel * wBeam**2)
    return SsurfT, rdel


def getCoatLayers(ifo, dOpt):
    """get layer vectors for refractive index, effective alpha and beta
    and geometrical thickness
    (see getCoatTOPos for example usage)
    
    ifo    = parameter struct from IFOmodel.m
    dOpt   = optical thickness / lambda of each layer
           = geometrical thickness * refractive index / lambda
    
    nLayer = refractive index of each layer, ordered input to output (N x 1)
    aLayer = change in geometrical thickness with temperature
           = the effective thermal expansion coeffient of the coating layer
    bLayer = change in refractive index with temperature
           = dn/dT
    dLayer = geometrical thicness of each layer
    sLayer = Yamamoto thermo-refractive correction
           = alpha * (1 + sigma) / (1 - sigma)"""
   
    # coating parameters
    lambda_ = ifo.Laser.Wavelength
    
    pS = ifo.Materials.Substrate
    pC = ifo.Materials.Coating
  
    Y_S = pS.MirrorY
    sigS = pS.MirrorSigma
  
    alphaL = pC.Alphalown
    betaL = pC.Betalown
    Y_L = pC.Ylown
    sigL = pC.Sigmalown
    nL = pC.Indexlown
  
    alphaH = pC.Alphahighn
    betaH = pC.Betahighn
    Y_H = pC.Yhighn
    sigH = pC.Sigmahighn
    nH = pC.Indexhighn
  
    Nlayer = len(dOpt)
  
    # compute effective alpha
    def getExpansionRatio(Y_C, sigC, Y_S, sigS):
        ##############################################
        # Y_C and sigC are for the coating material (can also be substrate)
        # Y_S and sigS are for the substrate material
        #
  
        ce = ((1 + sigS) / (1 - sigC)) \
             * ( ((1 + sigC) / (1 + sigS)) + (1 - 2 * sigS) * Y_C / Y_S )
        return ce

    aLayer = zeros((Nlayer, 1))
    aLayer[::2] = alphaL * getExpansionRatio(Y_L, sigL, Y_S, sigS)
    aLayer[1::2] = alphaH * getExpansionRatio(Y_H, sigH, Y_S, sigS)
  
    # and beta
    bLayer = zeros((Nlayer, 1))
    bLayer[::2] = betaL
    bLayer[1::2] = betaH

    # and refractive index
    nLayer = zeros((Nlayer, 1))
    nLayer[::2] = nL
    nLayer[1::2] = nH
  
    # and geometrical thickness
    dLayer = lambda_ * dOpt / nLayer

    # and sigma correction
    sLayer = zeros((Nlayer, 1))
    sLayer[::2] = alphaL * (1 + sigL) / (1 - sigL)
    sLayer[1::2] = alphaH * (1 + sigH) / (1 - sigH)

    return nLayer, aLayer, bLayer, dLayer, sLayer


def getCoatAvg(ifo, dOpt):
    """get coating average properties
    (see getCoatTOPos for example usage)
    
    ifo  = parameter struct from IFOmodel.m
    dOpt = optical thickness / lambda of each layer
         = geometrical thickness * refractive index / lambda
    
    dc = total thickness (meters)
    Cc = heat capacity
    Kc = thermal diffusivity
    aSub = effective substrate thermal expansion (weighted by heat capacity)"""

  
    # coating parameters
    pS = ifo.Materials.Substrate
    pC = ifo.Materials.Coating
  
    alphaS = pS.MassAlpha
    C_S = pS.MassCM * pS.MassDensity
    sigS = pS.MirrorSigma
  
    C_L = pC.CVlown
    K_L = pC.ThermalDiffusivitylown

    C_H = pC.CVhighn
    K_H = pC.ThermalDiffusivityhighn
  
    # compute refractive index, effective alpha and beta
    junk1, junk2, junk3, dLayer, junk4 = getCoatLayers(ifo, dOpt)
  
    # heat capacity
    dc = sum(dLayer)
    dL = sum(dLayer[::2])
    dH = sum(dLayer[1::2])
    Cc = (C_L * dL + C_H * dH) / dc
  
    # thermal diffusivity
    KinvL = 1 / K_L
    KinvH = 1 / K_H
    Kc = dc / (KinvL * dL + KinvH * dH)
  
    # effective substrate thermal expansion
    aSub = 2 * alphaS * (1 + sigS) * Cc / C_S

    return dc, Cc, Kc, aSub


def getCoatTOPhase(nIn, nOut, nLayer, dOpt, aLayer, bLayer, sLayer):
    """returns coating reflection phase derivatives w.r.t. temperature
    
    nIn = refractive index of input medium (e.g., vacuum = 1)
    nOut = refractive index of output medium (e.g., SiO2 = 1.45231 @ 1064nm)
    nLayer = refractive index of each layer, ordered input to output (N x 1)
    dOpt   = optical thickness / lambda of each layer
           = geometrical thickness * refractive index / lambda
    aLayer = change in geometrical thickness with temperature
           = the effective thermal expansion coeffient of the coating layer
    bLayer = change in refractive index with temperature
           = dn/dT 
           = dd/dT - n * a
    
    dphi_dT = total thermo-optic phase derivative with respect to temperature
            = dphi_TE + dphi_TR
    dphi_TE = thermo-elastic phase derivative (dphi / dT)
    dphi_TR = thermo-refractive phase derivative (dphi / dT)
    rCoat = amplitude reflectivity of coating (complex)
    
    Note about aLayer: on a SiO2 substrate,
    a_Ta2O5 ~ 3.5 * alpha_Ta2O5
    a_SiO2 ~ 2.3 * alpha_SiO2
    
    see getCoatTOPos for more information
    (see also T080101)"""

    # vector of all refractive indexs
    nAll = np.vstack((nIn, nLayer, nOut))
  
    # reflectivity of each interface
    r = (nAll[:-1] - nAll[1:]) / (nAll[:-1] + nAll[1:])
  
    # combine reflectivities
    rbar = zeros(r.shape, dtype=complex)
    ephi = zeros(r.shape, dtype=complex)

    ephi[-1] = exp(-4j * pi * dOpt[-1])
    rbar[-1] = ephi[-1] * r[-1]
    for n in range(len(dOpt), 0, -1):
        # one-way phase in this layer
        if n > 1:
            ephi[n-1] = exp(-4j * pi * dOpt[n - 2])
        else:
            ephi[n-1] = 1
    
        # accumulate reflecitivity
        rbar[n-1] = ephi[n-1] * (r[n-1] + rbar[n]) / (1 + r[n-1] * rbar[n])
  
    # reflectivity derivatives
    dr_dphi = zeros(dOpt.shape, dtype=complex)
    for n in range(len(dOpt), 0, -1):
        dr_dphi[n-1] = -1j * rbar[n]
        for m in range(n, 0, -1):
            dr_dphi[n-1] = dr_dphi[n-1] * ephi[m-1] * (1 - r[m-1]**2) / (1 + r[m-1] * rbar[m])**2
  

    # geometrical distances
    dGeo = dOpt / nLayer

    # phase derivatives
    dphi_dd = 4 * pi * imag(dr_dphi / rbar[0])

    # thermo-refractive coupling
    dphi_TR = sum(dphi_dd * (bLayer + sLayer * nLayer) * dGeo)

    # thermo-elastic
    dphi_TE = 4 * pi * sum(aLayer * dGeo)

    # total
    dphi_dT = dphi_TR + dphi_TE
  
    # coating reflectivity
    rCoat = rbar[0]

    return dphi_dT, dphi_TE, dphi_TR, rCoat


def getCoatFiniteCorr(ifo, wBeam, dOpt):
    """finite mirror size correction
    Uses correction factor from PLA 2003 vol 312 pg 244-255
    "Thermodynamical fluctuations in optical mirror coatings"
    by V. B. Braginsky and S. P. Vyatchanin
    http://arxiv.org/abs/cond-mat/0302617
    
    ifo = parameter struct from IFOmodel.m
    opticName = name of the Optic struct to use for wBeam and dOpt
    wBeam = ifoArg.Optics.(opticName).BeamRadius
    dOpt = ifoArg.Optics.(opticName).CoatLayerOpticalThickness
    
    wBeam = beam radius (at 1 / e^2 power)
    dOpt   = optical thickness / lambda of each layer
           = geometrical thickness * refractive index / lambda
    
    (see getCoatTOPos for example usage)
    
    version 1 by Sam Wald, 2008"""

    # parameter extraction
    R = ifo.Materials.MassRadius      #substrate radius
    H = ifo.Materials.MassThickness   #substrate thickness
    lambda_ = ifo.Laser.Wavelength
    zeta = ifo.Constants.BesselZeros  # zeros of 1st order bessel function (J1)

    alphaS = ifo.Materials.Substrate.MassAlpha
    C_S = ifo.Materials.Substrate.MassCM * ifo.Materials.Substrate.MassDensity
    Y_S = ifo.Materials.Substrate.MirrorY
    sigS = ifo.Materials.Substrate.MirrorSigma

    alphaL = ifo.Materials.Coating.Alphalown
    C_L = ifo.Materials.Coating.CVlown
    Y_L = ifo.Materials.Coating.Ylown
    sigL = ifo.Materials.Coating.Sigmalown
    nL = ifo.Materials.Coating.Indexlown

    alphaH = ifo.Materials.Coating.Alphahighn
    C_H = ifo.Materials.Coating.CVhighn
    Y_H = ifo.Materials.Coating.Yhighn
    sigH = ifo.Materials.Coating.Sigmahighn
    nH = ifo.Materials.Coating.Indexhighn

    # coating sums
    dL = lambda_ * sum(dOpt[::2]) / nL
    dH = lambda_ * sum(dOpt[1::2]) / nH
    dc = dH + dL

    # AVERAGE SPECIFIC HEAT (simple volume average for coating)
    Cf = (C_L * dL + C_H * dH) / dc
    Cr = Cf / C_S

    # COATING AVERAGE VALUE X = ALPHAF*(1+POISSONf)/(1-POISSONf) avg
    xxL = alphaL * (1 + sigL) / (1 - sigL)
    xxH = alphaH * (1 + sigH) / (1 - sigH)
    Xf = (xxL * dL + xxH * dH) / dc
    Xr = Xf / alphaS

    # COATING AVERAGE VALUE Y = ALPHAF* YOUNGSF/(1-POISSONF) avg
    yyL = alphaL * Y_L / (1 - sigL)
    yyH = alphaH * Y_H / (1 - sigH)
    Yf = (yyL * dL + yyH * dH) / dc
    Yr = Yf / (alphaS * Y_S)

    # COATING AVERAGE VALUE Z = 1/(1-POISSONF) avg
    zzL = 1 / (1 - sigL)
    zzH = 1 / (1 - sigH)
    Zf = (zzL * dL + zzH * dH) / dc

    #################################### FINITE SIZE CORRECTION CALCULATION

    # beam size parameter used by Braginsky
    r0 = wBeam / sqrt(2)

    # values of J0 at zeros of J1
    j0m = scipy.special.jn(0, zeta)

    # between eq 77 and 78
    km = zeta / R
    Qm = exp(-2 * km * H)
    pm = exp(-km**2 * r0**2 / 4) / j0m;  # left out factor of pi * R^2 in denominator

    # eq 88
    Lm = Xr - Zf * (1 + sigS) + (Yr * (1 - 2 * sigS) + Zf - 2 * Cr) * \
         (1 + sigS) * (1 - Qm)**2 / ((1 - Qm)**2 - 4 * km**2 * H**2 * Qm)

    # eq 90 and 91
    S1 = (12 * R**2 / H**2) * sum(pm / zeta**2)
    S2 = sum(pm**2 * Lm**2)
    P = (Xr - 2 * sigS * Yr - Cr + S1 * (Cr - Yr * (1 - sigS)))**2 + S2

    # eq 60 and 70
    LAMBDA = -Cr + (Xr / (1 + sigS) + Yr * (1 - 2 * sigS)) / 2

    # eq 92
    Cfsm = sqrt((r0**2 * P) / (2 * R**2 * (1 + sigS)**2 * LAMBDA**2))
    return Cfsm
