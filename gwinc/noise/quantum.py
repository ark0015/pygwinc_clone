from __future__ import division
from numpy import pi, sqrt, arctan, sin, cos, exp, size, ones, zeros, log10, conj, sum
import numpy as np
import scipy.constants
import logging
from ..struct import Struct


def shotrad(f, ifo):
    """Quantum noise

    corresponding author: mevans
    modifications for resonant delay lines: Stefan Ballmer

    """

    try:
        sqzType = ifo.Squeezer.Type
    except AttributeError:
        sqzType = None

    #if sqzType == '2mode':
    #    from .quantum_2mode import shotrad as shotrad_2mode
    #    return shotrad_2mode(f, ifo)

    # deal with multiple bounces, required for resonant delay lines
    # Stefan Ballmer 2012
    if 'NFolded' in ifo.Infrastructure:
        if ifo.Infrastructure.travellingWave:
            ifo.Materials.MirrorMass=ifo.Materials.MirrorMass/ifo.Infrastructure.NFolded**2
        else:
            NN=ifo.Infrastructure.NFolded
            Ndispl=2 * NN-1
            ifo.Materials.MirrorMass=ifo.Materials.MirrorMass/Ndispl**2

    #####################################################
    # Call IFO Quantum Model
    #####################################################

    if 'Type' not in ifo.Optics:
        fname = shotradSignalRecycled
    else:
        namespace = globals()
        fname = namespace['shotrad' + ifo.Optics.Type]
    coeff, Mifo, Msig, Mn = fname(f, ifo)

    # check for consistent dimensions
    Nfield = Msig.shape[0]
    Nfreq = len(f)
    if any(np.array([Mifo.shape[0], Mifo.shape[1], Mn.shape[0]]) != Nfield) or \
       any(np.array([Mifo.shape[2], Msig.shape[2], Mn.shape[2]]) != Nfreq):
        logging.debug(Mifo.shape)
        logging.debug(Msig.shape)
        logging.debug(Mn.shape)
        logging.debug(Nfield, Nfreq)
        raise Exception('Inconsistent matrix sizes returned by %s' % str(fname))

    # deal with non-standard number of fields
    if Nfield != 2:
        if Nfield == 4:
            n = shotrad4(f, ifo, coeff, Mifo, Msig, Mn)
            return n
        else:
            raise Exception("shotrad doesn't know what to do with %d fields" % Nfield)

    #####################################################
    # Input Squeezing
    #####################################################

    # ------------------------------------------- equation 63 BnC PRD 2004
    #>>>>>>>>    QUANTUM NOISE POWER SPECTRAL DENSITY WITH SQZ [BnC PRD 2004, 62]
    #<<<<<<<<<<<<<<<<< Modified to include losses (KM)
    #<<<<<<<<<<<<<<<<< Modified to include frequency dependent squeezing angle (LB)
    # useful numbers
    #TODO, adjust Struct to allow deep defaulted access
    # Homodyne Readout phase
    eta_orig = ifo.Optics.get('Quadrature', Struct()).get('dc', None)

    ifoRead = ifo.get('Squeezer', Struct()).get('Readout', None)
    if ifoRead is None:
        eta = eta_orig
        if eta_orig is None:
            raise Exception("must add Quadrature.dc or Readout...")
    elif ifoRead.Type == 'DC':
        eta = np.sign(ifoRead.fringe_side) * np.arccos((ifoRead.defect_PWR_W / ifoRead.readout_PWR_W)**.5)
    elif ifoRead.Type == 'Homodyne':
        eta = ifoRead.Angle
    else:
        raise Exception("Unknown Readout Type")

    if eta_orig is not None:
        # logging.warn((
        #     'Quadrature.dc is redundant with '
        #     'Squeezer.Readout and is deprecated.'
        # ))
        if eta_orig != eta:
            raise Exception("Quadrature.dc inconsistent with Readout eta")

    lambda_PD = 1 - ifo.Optics.PhotoDetectorEfficiency  # PD losses

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # determine squeezer type, if any
    # and extract common parameters
    if 'Squeezer' not in ifo:
        sqzType = 'None'
    else:
        sqzType = ifo.Squeezer.get('Type', 'Freq Independent')

    # extract common parameters
    if sqzType == 'None':
        SQZ_DB = 0                               # Squeeing in dB
        alpha = 0                                # Squeeze angle
        lambda_in = 0                            # Loss to squeezing before injection [Power]
        ANTISQZ_DB = 0                           # Anti squeezing in db
        etaRMS = 0
    else:
        SQZ_DB = ifo.Squeezer.AmplitudedB        # Squeeing in dB
        lambda_in = ifo.Squeezer.InjectionLoss   # Loss to squeezing before injection [Power]
        alpha = ifo.Squeezer.SQZAngle        # Freq Indep Squeeze angle
        ANTISQZ_DB = ifo.Squeezer.get('AntiAmplitudedB', SQZ_DB)  # Anti squeezing in db
        etaRMS = ifo.Squeezer.get('LOAngleRMS', 0)  # quadrature noise

    # switch on squeezing type for other input squeezing modifications
    if sqzType == 'None':
        pass

    elif sqzType == 'Freq Independent':
        logging.debug('You are injecting %g dB of frequency independent squeezing' % SQZ_DB)

    elif sqzType == 'Optimal':
        # compute optimal squeezing angle
        alpha = sqzOptimalSqueezeAngle(Mifo, eta)

        logging.debug('You are injecting %g dB of squeezing with optimal frequency dependent squeezing angle' % SQZ_DB)

    elif sqzType == 'OptimalOptimal':
        # compute optimal squeezing angle, assuming optimal readout phase
        R = SQZ_DB / (20 * log10(exp(1)))
        MnPD = sqzInjectionLoss(Mn, lambda_PD)
        MsigPD = Msig * sqrt(1 - lambda_PD)
        alpha = sqzOptimalSqueezeAngle(Mifo, [], [R, lambda_in], MsigPD, MnPD)

        logging.debug('You are injecting %g dB of squeezing with optimal FD squeezing angle, for optimal readout phase' % SQZ_DB)

    elif sqzType == 'Freq Dependent':
        logging.debug('You are injecting %g dB of squeezing with frequency dependent squeezing angle' % SQZ_DB)

    else:
        raise Exception('ifo.Squeezer.Type must be None, Freq Independent, Optimal, or Frequency Dependent, not "%s"' % sqzType)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Define input matrix of Squeezing
    R = SQZ_DB / (20 * log10(exp(1)))                 # Squeeze factor
    R_anti = ANTISQZ_DB / (20 * log10(exp(1)))        # Squeeze factor
    Msqz = np.array([[exp(-R), 0], [0, exp(R_anti)]])

    # expand to Nfreq
    Msqz = np.transpose(np.tile(Msqz, (Nfreq,1,1)), axes=(1,2,0))

    # add input rotation
    MsqzRot = make2x2TF(cos(alpha), -sin(alpha), sin(alpha), cos(alpha))
    Msqz = getProdTF(MsqzRot, Msqz)

    # cheat to test optimal squeezing agle code
    #   if strcmp(sqzType, 'Optimal') || strcmp(sqzType, 'OptimalOptimal')
    #     Msqz = [exp(-R) 0; 0 exp(-R)];
    #   end

    # Include losses (lambda_in=ifo.Squeezer.InjectionLoss)
    Msqz = sqzInjectionLoss(Msqz, lambda_in)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Inject squeezed field into the IFO via some filter cavities
    if sqzType == 'Freq Dependent' and 'FilterCavity' in ifo.Squeezer:
        logging.debug('  Applying %d input filter cavities' % np.atleast_1d(ifo.Squeezer.FilterCavity).size)
        Mr, Msqz = sqzFilterCavityChain(f, np.atleast_1d(ifo.Squeezer.FilterCavity), Msqz)

    #####################################################
    # IFO Transfer and Output Filter Cavities
    #####################################################

    # apply the IFO dependent squeezing matrix to get
    #   the total noise due to quantum fluctuations coming in from the AS port
    Msqz = getProdTF(Mifo, Msqz)

    # add this to the other noises Mn
    Mnoise = np.hstack((Msqz, Mn))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # pass IFO output through some filter cavities
    if 'OutputFilter' in ifo:
        if ifo.OutputFilter.Type == 'None':
            # do nothing, say nothing
            pass

        elif ifo.OutputFilter.Type == 'Chain':
            logging.debug('  Applying %d output filter cavities' % np.atleast_1d(ifo.OutputFilter.FilterCavity).size)

            Mr, Mnoise = sqzFilterCavityChain(f, np.atleast_1d(ifo.OutputFilter.FilterCavity), Mnoise)
            Msig = getProdTF(Mr, Msig)
            #  Mnoise = getProdTF(Mn, Mnoise);

        elif ifo.OutputFilter.Type == 'Optimal':
            logging.debug('  Optimal output filtering!')

            # compute optimal angle, including upcoming PD losses
            MnPD = sqzInjectionLoss(Mnoise, lambda_PD)
            raise NotImplementedError("Cannot do optimal phase yet")
            zeta = sqzOptimalReadoutPhase(Msig, MnPD)

            # rotate by that angle, less the homodyne angle
            #zeta_opt = eta;
            cs = cos(zeta - eta)
            sn = sin(zeta - eta)
            Mrot = make2x2TF(cs, -sn, sn, cs)
            Mnoise = getProdTF(Mrot, Mnoise)
            Msig = getProdTF(Mrot, Msig)

        else:
            raise Exception('ifo.OutputFilter.Type must be None, Chain or Optimal, not "%s"' % ifo.OutputFilter.Type)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # add PD efficiency
    Mnoise = sqzInjectionLoss(Mnoise, lambda_PD)
    Msig = Msig * sqrt(1 - lambda_PD)

    # and compute the final noise
    def HDnoise(eta):
        vHD = np.array([[sin(eta), cos(eta)]])
        n = coeff * np.squeeze(sum(abs(getProdTF(vHD, Mnoise))**2, axis=1)) / \
            np.squeeze(sum(abs(getProdTF(vHD, Msig))**2, axis=1))
        return n

    if etaRMS <= 0:
        n = HDnoise(eta)
    else:
        # include quadrature noise (average of +- the RMS)
        n = (HDnoise(eta+etaRMS) + HDnoise(eta-etaRMS)) / 2

    # the above is the same as
    #    n = coeff * (vHD * Msqz * Msqz' * vHD') / (vHD * Msig * Msig' * vHD')
    #  where ' is the conjugate transpose operation.  Which is also
    #    n = coeff * sym(vHD * Msqz) / sym(vHD * Msig)
    #  where is the symmeterization operation
    #    sym(M) = real(M * M')
    #
    # it is also the same as taking the sum of the squared directly
    #   n = zeros(1, numel(f));
    #   for k = 1:numel(f)
    #     n(k) = coeff(k) * sum(abs((vHD * Msqz(:,:,k))).^2) ./ ...
    #       sum(abs((vHD * Msig(:,:,k))).^2);
    #   end

    return n * ifo.gwinc.sinc_sqr


def shotradSignalRecycled(f, ifo):
    """Quantum noise model for signal recycled IFO

    See shotrad for more info.

    All references to Buonanno & Chen PRD 64 042006 (2001) (hereafter BnC)
    Updated to include losses DEC 2006 Kirk McKenzie using BnC notation
    Updated to include squeezing April 2009 KM
    Updated April 2010 KM, LB

    moved out of shotrad May 2010, mevans
    output is used in shotrad to compute final noise as
      n = coeff * (vHD * Msqz * Msqz' * vHD') / (vHD * Md * Md' * vHD')
    where
      Msqz = [Mc MsqueezeInput, Mn]

    coeff = frequency dependent overall noise coefficient (Nx1)
    Mifo = IFO input-output relation for the AS port
    Msig = signal transfer to the AS port
    Mnoise = noise fields produced by losses in the IFO at the AS port

    """
    # f                                           % Signal Freq. [Hz]
    lambda_ = ifo.Laser.Wavelength               # Laser Wavelength [m]
    hbar    = scipy.constants.hbar               # Plancks Constant [Js]
    c       = scipy.constants.c                  # SOL [m/s]
    Omega   = 2*pi*f                             # [BnC, table 1] Signal angular frequency [rads/s]
    omega_0 = 2*pi*c/lambda_                     # [BnC, table 1] Laser angular frequency [rads/s]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    L       = ifo.Infrastructure.Length          # Length of arm cavities [m]
    l       = ifo.Optics.SRM.CavityLength        # SRC Length [m]
    T       = ifo.Optics.ITM.Transmittance       # ITM Transmittance [Power]
    m       = ifo.Materials.MirrorMass           # Mirror mass [kg]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    bsloss  = ifo.Optics.BSLoss                  # BS Loss [Power]
    mismatch = 1 - ifo.Optics.coupling           # Mismatch
    mismatch = mismatch + ifo.TCS.SRCloss        # Mismatch

    # BSloss + mismatch has been incorporated into a SRC Loss
    lambda_SR = mismatch + bsloss                # SR cavity loss [Power]

    tau     = sqrt(ifo.Optics.SRM.Transmittance) # SRM Transmittance [amplitude]
    rho     = sqrt(1 - tau**2 - lambda_SR)       # SRM Reflectivity [amplitude]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ds      = ifo.Optics.SRM.Tunephase           # SRC Detunning
    phi     = (pi-ds)/2                          # [BnC, between 2.14 & 2.15] SR Detuning

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    lambda_arm = ifo.Optics.Loss*2               # [BnC, after 5.2] Round Trip loss in arm [Power]
    gamma_ac = T*c/(4*L)                         # [KLMTV-PRD2001] Arm cavity half bandwidth [1/s]
    epsilon = lambda_arm/(2*gamma_ac*L/c)        # [BnC, after 5.2] Loss coefficent for arm cavity

    I_0     = ifo.gwinc.pbs                      # [BnC, Table 1] Power at BS (Power*prfactor) [W]
    I_SQL   = (m*L**2*gamma_ac**4)/(4*omega_0)   # [BnC, 2.14] Power to reach free mass SQL
    Kappa   = 2*((I_0/I_SQL)*gamma_ac**4)/ \
              (Omega**2*(gamma_ac**2+Omega**2))  # [BnC 2.13] Effective Radiation Pressure Coupling
    beta    = arctan(Omega/gamma_ac)             # [BnC, after 2.11] Phase shift of GW SB in arm
    h_SQL   = sqrt(8*hbar/(m*(Omega*L)**2))      # [BnC, 2.12] SQL Strain


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Coefficients [BnC, Equations 5.8 to 5.12]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    exp_1jbeta = exp(1j*beta)
    cos_beta = exp_1jbeta.real
    invexp_1jbeta = 1/exp_1jbeta
    exp_2jbeta = exp_1jbeta**2
    cos_2beta = exp_2jbeta.real
    invexp_2jbeta = 1/exp_2jbeta
    exp_4jbeta = exp_2jbeta**2
    C11_L   = ( (1+rho**2) * ( cos(2*phi) + Kappa/2 * sin(2*phi) ) -
                2*rho*cos_2beta - 1/4*epsilon * ( -2 * (1+exp_2jbeta)**2 * rho + 4 * (1+rho**2) *
                                                    cos_beta**2*cos(2*phi) + ( 3+exp_2jbeta ) *
                                                    Kappa * (1+rho**2) * sin(2*phi) ) +
                lambda_SR * ( exp_2jbeta*rho-1/2 * (1+rho**2) * ( cos(2*phi)+Kappa/2 * sin(2*phi) ) ) )

    C22_L   = C11_L

    C12_L   = tau**2 * ( - ( sin(2*phi) + Kappa*sin(phi)**2 )+
                         1/2*epsilon*sin(phi) * ( (3+exp_2jbeta) * Kappa * sin(phi) + 4*cos_beta**2 * cos(phi)) +
                         1/2*lambda_SR * ( sin(2*phi)+Kappa*sin(phi)**2) )

    C21_L   = tau**2 * ( (sin(2*phi)-Kappa*cos(phi)**2 ) +
                         1/2*epsilon*cos(phi) * ( (3+exp_2jbeta )*Kappa*sin(phi) - 4*cos_beta**2*sin(phi) ) +
                         1/2*lambda_SR * ( -sin(2*phi) + Kappa*cos(phi)**2) )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    D1_L    = ( - (1+rho*exp_2jbeta ) * sin(phi) +
                1/4*epsilon * ( 3+rho+2*rho*exp_4jbeta + exp_2jbeta*(1+5*rho) ) * sin(phi)+
                1/2*lambda_SR * exp_2jbeta * rho * sin(phi) )

    D2_L    = ( - (-1+rho*exp_2jbeta ) * cos(phi) +
                1/4*epsilon * ( -3+rho+2*rho*exp_4jbeta + exp_2jbeta * (-1+5*rho) ) * cos(phi)+
                1/2*lambda_SR * exp_2jbeta * rho * cos(phi) )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    P11     = 0.5 * sqrt(lambda_SR) * tau * \
              ( -2*rho*exp_2jbeta+2*cos(2*phi)+Kappa*sin(2*phi) )
    P22     = P11
    P12     = -sqrt(lambda_SR)*tau*sin(phi)*(2*cos(phi)+Kappa*sin(phi) )
    P21     =  sqrt(lambda_SR)*tau*cos(phi)*(2*sin(phi)-Kappa*cos(phi) )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # this was the PD noise source, but that belongs outside of this function
    #   I have used the equation for Q11 to properly normalize the other noises
    #   as well as the input-output relation Mc and the signal matrix Md

    Q11     = 1 / \
              ( invexp_2jbeta+rho**2*exp_2jbeta-rho*(2*cos(2*phi)+Kappa*sin(2*phi)) +
                1/2*epsilon*rho * (invexp_2jbeta*cos(2*phi)+exp_2jbeta*
                                   ( -2*rho-2*rho*cos_2beta+cos(2*phi)+Kappa*sin(2*phi) ) +
                                   2*cos(2*phi)+3*Kappa*sin(2*phi))-1/2*lambda_SR*rho *
                ( 2*rho*exp_2jbeta-2*cos(2*phi)-Kappa*sin(2*phi) ) )
    Q22     = Q11
    Q12     = 0
    Q21     = 0

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    N11     = sqrt(epsilon/2)*tau *(Kappa*(1+rho*exp_2jbeta)*sin(phi)+
                                    2*cos_beta*(invexp_1jbeta*cos(phi)-rho*exp_1jbeta*(cos(phi)+Kappa*sin(phi))))
    N22     = -sqrt(2*epsilon)*tau*(-invexp_1jbeta+rho*exp_1jbeta)*cos_beta*cos(phi)
    N12     = -sqrt(2*epsilon)*tau*(invexp_1jbeta+rho*exp_1jbeta)*cos_beta*sin(phi);
    N21     = sqrt(2*epsilon)*tau*(-Kappa*(1+rho)*cos(phi)+
                                   2*cos_beta*(invexp_1jbeta+rho*exp_1jbeta)*cos_beta*sin(phi))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # overall coefficient
    coeff = h_SQL**2/(2*Kappa*tau**2)

    # make normalization matrix
    Mq = make2x2TF(Q11, Q12, Q21, Q22)

    # 3D transfer matrices from vectors for each element
    Mifo = getProdTF(Mq, make2x2TF(C11_L, C12_L, C21_L, C22_L))
    Msig = getProdTF(Mq, np.array([D1_L, D2_L]).reshape(2, 1, f.size))

    # put all output noises together
    Mp = make2x2TF(P11, P12, P21, P22)
    Mn = make2x2TF(N11, N12, N21, N22)
    Mnoise = getProdTF(Mq, np.hstack((Mn, Mp)))

    return coeff, Mifo, Msig, Mnoise


def make2x2TF(A11, A12, A21, A22):
    """
    Create transfer matrix with 2x2xnF.
    """
    #now using numpy with broadcasting
    A11, A12, A21, A22 = np.broadcast_arrays(A11, A12, A21, A22)
    M3 = np.array([[A11, A12], [A21, A22]])
    return M3.reshape(2, 2, -1)


def getProdTF(lhs, rhs):
    """Compute the product of M Nout x Nin x Naf frequency dependent transfer matrices

    See also getTF.

    NOTE: To perform more complicated operations on transfer
          matrices, see LTI object FRD ("help frd").  This
          function is the same as: freqresp(frd(lhs) * frd(rhs), f)

    """
    # check matrix size
    if lhs.shape[1] != rhs.shape[0]:
        raise Exception('Matrix size mismatch size(lhs, 2) = %d != %d = size(rhs, 1)' % (lhs.shape[1], rhs.shape[0]))
    N = lhs.shape[0]
    M = rhs.shape[1]
    if len(lhs.shape) == 3:
        lhs = np.transpose(lhs, axes=(2, 0, 1))
    if len(rhs.shape) == 3:
        rhs = np.transpose(rhs, axes=(2, 0, 1))

    # compute product
    if len(lhs.shape) < 3 or lhs.shape[0] == 1:
        rslt = np.matmul(lhs, rhs)
    elif len(rhs.shape) < 3 or rhs.shape[0] == 1:
        rslt = np.matmul(lhs, rhs)
    elif lhs.shape[0] == rhs.shape[0]:
        rslt = np.matmul(lhs, rhs)
    else:
        raise Exception('Matrix size mismatch lhs.shape[2] = %d != %d = rhs.shape[2]' % (lhs.shape[2], rhs.shape[2]))

    if len(rslt.shape) == 3:
        rslt = np.transpose(rslt, axes=(1, 2, 0))
    return rslt


def sqzInjectionLoss(Min, L):
    """Injection losses for squeezed field

    lambda_in is defined as ifo.Squeezer.InjectionLoss

    """
    eye2 = np.eye(Min.shape[0], Min.shape[1])
    Meye = np.transpose(np.tile(eye2, (Min.shape[2],1,1)), axes=(1,2,0))

    Mout = np.hstack((Min * sqrt(1 - L), Meye * sqrt(L)))
    return Mout


def sqzFilterCavityChain(f, params, Mn):
    """Transfer relation for a chain of filter cavities

    Noise added by cavity losses are also output.

    f = frequency vector [Hz]
    param.fdetune = detuning [Hz]
    param.L = cavity length
    param.Ti = input mirror trasmission [Power]
    param.Li = input mirror loss
    param.Te = end mirror trasmission
    param.Le = end mirror loss
    param.Rot = phase rotation after cavity

    Mn0 = input noise
    Mc = input to output transfer
    Mn = filtered input noise, plus noise due to cavity losses

    Note:
        [Mc, Mn] = sqzFilterCavityChain(f, params, Mn0)
      is the same as
        [Mc, Mn] = sqzFilterCavityChain(f, params);
        Mn = [getProdTF(Mc, Mn0), Mn];

    corresponding author: mevans

    """
    # make an identity TF
    Mc = make2x2TF(ones(f.shape), 0, 0, 1)

    # loop through the filter cavites
    for k in range(params.size):
        # extract parameters for this filter cavity
        Lf = params[k].L
        fdetune = params[k].fdetune
        Ti = params[k].Ti
        Te = params[k].Te
        Lrt = params[k].Lrt
        theta = params[k].Rot

        # compute new Mn
        Mr, Mt, Mn = sqzFilterCavity(f, Lf, Ti, Te, Lrt, fdetune, Mn)

        # apply phase rotation after filter cavity
        Mrot = np.array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])
        Mn = getProdTF(Mrot, Mn)

        # update Mc
        Mc = getProdTF(Mrot, getProdTF(Mr, Mc))

    return Mc, Mn


def sqzFilterCavity(f, Lcav, Ti, Te, Lrt, fdetune, MinR, MinT=1):
    """Reflection/transmission matrix for filter cavity

    Function which gives the reflection matrix for vacuum fluctuations
    entering the input mirror and the transmission matrix for vacuum
    fluctuations entering the end mirror of one filter cavity.  The
    input parameters are the cavity parameters and the 2X2 matrix of
    the incoming fields in the two-photon formalism.

    (R_alpha x S_r) for a freq independent squeezed field.
    f = vector frequency in Hz
    Lf = length of the filter cavity
    Ti = transmission and losses of the input mirror
    Te = transmission and losses of the end mirror
    Lrt = round-trip losses in the cavity (mirror transmissoins not included)
    fdetune: detuning frequency of the filter cavity [Hz]
    MinR: squeezed field injected from the input mirror of the filter cavity (a1,a2 basis)
         if this argument is empty, it is assumed that the user will use Mr,
         so no noise field is added to Mnoise.  If no argument is given, or
         the scalar 1 is given, an Mr unsqueezed input is assumed and Mr is
         concatenated into Mnoise.
    MinT: squeezed field injected from the back of the filter cavity
         with MinR, this argument can be omitted or set to 1 to indicate
         and unsqueezed input. [] can be used to avoid adding a noise
         term to Mnoise.

    corresponding authors: LisaB, mevans

    """

    # reflectivities
    Ri = 1 - Ti
    Re = 1 - Te

    ri = sqrt(Ri)
    re = sqrt(Re)
    rr = ri * re * sqrt(1 - Lrt)  # include round-trip losses

    # Phases for positive and negative audio sidebands
    c = scipy.constants.c
    omega = 2 * pi * f
    wf = 2 * pi * fdetune
    Phi_p = 2 * (omega-wf)* Lcav / c
    Phi_m = 2 * (-omega-wf)* Lcav / c

    ephi_p = exp(1j * Phi_p)
    ephi_m = exp(1j * Phi_m)

    # cavity gains
    g_p = 1 / ( 1 - rr * ephi_p)
    g_m = 1 / ( 1 - rr * ephi_m)

    # Reflectivity for vacuum flactuation entering the cavity from
    # the input mirror (check sign)
    r_p = ri - re * Ti * ephi_p * g_p
    r_m = ri - re * Ti * ephi_m * g_m


    # Transmissivity for vacuum flactuation entering the cavity from
    # the back mirror (check sign)
    t_p = sqrt(Ti * Te * ephi_p) * g_p
    t_m = sqrt(Ti * Te * ephi_m) * g_m

    # Transmissivity for vacuum flactuation entering the cavity from
    # the losses in the cavity
    l_p = sqrt(Ti * Lrt * ephi_p) * g_p
    l_m = sqrt(Ti * Lrt * ephi_m) * g_m

    # Relfection matrix for vacuum fluctuations entering from the input
    # mirror in the A+, (a-)* basis
    Mr_temp = make2x2TF(r_p, 0, 0, conj(r_m))

    # Transmission matrix for vacuum fluctuations entering from the end mirror
    Mt_temp = make2x2TF(t_p, 0, 0, conj(t_m))

    # Transmission matrix for vacuum fluctuations entering from the end mirror
    Ml_temp = make2x2TF(l_p, 0, 0, conj(l_m))

    # Apply matrix which changes from two-photon basis to a+ and (a-)*
    Mbasis = np.array([[1, 1j], [1, -1j]])

    Mr = getProdTF(np.linalg.inv(Mbasis), getProdTF(Mr_temp, Mbasis))
    Mt = getProdTF(np.linalg.inv(Mbasis), getProdTF(Mt_temp, Mbasis))
    Ml = getProdTF(np.linalg.inv(Mbasis), getProdTF(Ml_temp, Mbasis))

    ###### output

    # reflected fields
    if MinR == []:
        Mnoise = zeros((2, 0, f.size))
    else:
        if np.isscalar(MinR):
            Mnoise = Mr * MinR
        else:
            Mnoise = getProdTF(Mr, MinR)

    # transmitted fields
    if MinT != [] and Te > 0:
        if np.isscalar(MinT) and MinT == 1:
            Mnoise = np.hstack((Mnoise, Mt))
        else:
            Mnoise = np.hstack((Mnoise, getProdTF(Mt, MinT)))

    # loss fields
    if Lrt > 0:
        Mnoise = np.hstack((Mnoise, Ml))

    return Mr, Mt, Mnoise
