from __future__ import division, print_function
import numpy as np
from numpy import log10, pi, sqrt
import copy

from .util import precompIFO
from . import noise
from . import plot


def noise_calc(ifo, f):
    """Calculate all IFO noises and return as dict

    """
    noises = {}

    # quad cases, for backward compatability
    if ifo.Suspension.Type in (0, 1, 2):
        hForce, vForce, hTable, vTable = noise.suspensionthermal.suspQuad(f, ifo)
    elif ifo.Suspension.Type == 'Quad':
        ifo.Suspension.Type = 0
        hForce, vForce, hTable, vTable = noise.suspensionthermal.suspQuad(f, ifo)
    else:
        fname = noise.suspensionthermal.__dict__['susp' + str(ifo.Suspension.Type)]
        hForce, vForce, hTable, vTable = fname(f, ifo)

    # if the suspension code supports different temps for the stages
    try:
        # full TF (conventional)
        ifo.Suspension.hForce = hForce.fullylossy
        ifo.Suspension.vForce = vForce.fullylossy
        # TFs with each stage lossy
        ifo.Suspension.hForce_singlylossy = hForce.singlylossy
        ifo.Suspension.vForce_singlylossy = vForce.singlylossy
    except:
        ifo.Suspension.hForce = hForce
        ifo.Suspension.vForce = vForce

    ifo.Suspension.hTable = hTable
    ifo.Suspension.vTable = vTable

    if 'Temp' not in ifo.Materials.Substrate:
        ifo.Materials.Substrate.Temp = ifo.Constants.Temp

    noises['Quantum Vacuum'] = noise.quantum.shotrad(f, ifo)
    noises['Suspension Thermal']  = noise.suspensionthermal.suspR(f, ifo)
    noises['Excess Gas']  = noise.residualgas.gas(f, ifo)
    noises['Substrate Brownian']  = noise.substratethermal.subbrownian(f, ifo)
    noises['Coating Brownian']  = noise.coatingthermal.coatbrownian(f, ifo)
    noises['Substrate Thermo-Elastic']  = noise.substratethermal.subtherm(f, ifo)
    noises['Newtonian Gravity']  = noise.newtonian.gravg(f, ifo)
    noises['Seismic'] = noise.seismic.seismic(f, ifo)
    noises['Coating Thermo-Optic'] = noise.coatingthermal.thermooptic(f, ifo)

    # calc semiconductor noise sources
    if 'isSemiConductor' in ifo.Materials.Substrate.__dict__ and ifo.Materials.Substrate.isSemiConductor:
        noises['ITM Thermo-Refractive'] = noise.substratethermal.thermorefractiveITM(f, ifo)
        noises['ITM Carrier Density'] = noise.substratethermal.carrierdensity(f, ifo)

    # adjust noise curves for multiple bounces - for resonant delay lines
    # Brownian noise scales as Neff (correlation-corrected spot number)
    # Displacement noises scale as N^2
    # Thermo-optic noise scales as N (incoherent between spots)
    if 'NFolded' in ifo.Infrastructure.__dict__:
        print('FOLDED')
        if ifo.Infrastructure.travellingWave:
            N = ifo.Infrastructure.NFolded
            sep_w = ifo.Infrastructure.DelayLineSpotSeparation
            Neff = getBrownianCorrelationFactor(N,sep_w) # not yet ported to python
            noises['Suspension Thermal'] *= N**2
            noises['Substrate Brownian'] *= Neff
            noises['Coating Brownian'] *= Neff
            noises['Substrate Thermo-Elastic'] *= N**2
            noises['Newtonian Gravity'] *= N**2
            noises['Seismic'] *= N**2
            noises['Coating Thermo-Optic'] *= N
        else:
            N = ifo.Infrastructure.NFolded
            sep_w = ifo.Infrastructure.DelayLineSpotSeparation
            Neff = getBrownianCorrelationFactorStandingWave(N,sep_w) # not yet ported to python
            Ndispl = (2*N-1)
            N_TO = 4*N-3 # use naive counting, needs improvement
            noises['Suspension Thermal'] *= Ndispl**2
            noises['Substrate Brownian'] *= Neff
            noises['Coating Brownian'] *= Neff
            noises['Substrate Thermo-Elastic'] *= Ndispl**2
            noises['Newtonian Gravity'] *= Ndispl**2
            noises['Seismic'] *= Ndispl**2
            noises['Coating Thermo-Optic'] *= N_TO
        #noises['Mirror Thermal'] = noises['Substrate Brownian'] + noises['Coating Brownian'] + \
        #                           noises['Substrate Thermo-Elastic'] + noises['Coating Thermo-Optic'] # total mirror thermal 

    noises['Total'] = sum(noises[curve] for curve in noises)
    noises['Freq'] = f

    return noises



def gwinc(flo, fhi, ifoin, source=None, fig=False):
    """Calculates strain noise due to various noise sources, for a
    specified set of interferometer parameters. Also evaluates the
    sensitivity of the interferometer to the detection of several potential 
    gravitational wave sources. Usage:
    
         VARARGOUT = GWINC(FLO,FHI,IFO,SOURCE,VARARGIN)
    
         FLO, FHI = minimum and maximum frequencies between which
                     calculations are made
         IFO       = structure containing interferometer parameters
         SOURCE    = structure containing source parameters
    
    Optional input arguments (the last 4 override IFO parameters):
         VARARGIN{1}: PLOT_FLAG set to 4 for score, only calculating shotrad
                                       3 for score and plots
                                       2 for score only
                                       1 to make plots but no score
                                       else 0 (DEF)
         VARARGIN{2}: LASER POWER -> ifo.Laser.Power
         VARARGIN{3}: SRC PHASE   -> ifo.Optics.SRM.Tunephase
         VARARGIN{4}: SRM TRANS   -> ifo.Optics.SRM.Transmittance
         VARARGIN{5}: ITM TRANS   -> ifo.Optics.ITM.Transmittance
         VARARGIN{6}: PRM TRANS   -> ifo.Optics.PRM.Transmittance
    
    Optional output arguments
         VARARGOUT{1}: SCORE  structure containing source sensitivities
         VARARGOUT{2}: NOISE  structure containing noise terms
    
    Ex.1    [score,noise] = gwinc(5,5000,IFOModel,SourceModel,1)"""

    ifo = copy.deepcopy(ifoin)
    # -------------------------------------------------------
    # parse arguments

    #makescore = 0;
    modeSR = 0
    PRfixed = 0

    # Parse varargin to decide to make plots or scores or both or neither

    # Stick it into the IFO so that it gets passed around
    ifo.modeSR = modeSR

    # Adjust these parameters as command line arguments

    # --------------------------------------------------------
    # add some precomputed info to the ifo struct
    ifo = precompIFO(ifo, PRfixed)

    pbs      = ifo.gwinc.pbs
    finesse  = ifo.gwinc.finesse
    prfactor = ifo.gwinc.prfactor
    if ifo.Laser.Power * prfactor != pbs:
        pass
        #warning(['Thermal lensing limits input power to ' num2str(pbs/prfactor, 3) ' W']);

    # Frequency grid on which everything is calculated
    if 'freq' in ifo.__dict__:
        f = ifo.freq
    else:
        f = np.logspace(log10(flo), log10(fhi), 250)

    noises = noise_calc(ifo, f)

    n = noises['Total']

    ifo.nse = noises
    retval = (n, noises)

    # Report astrophysical scores if so desired
    if source:
        sss = int73(f, n, ifo, source)
        sss.Omega = intStoch(f, n, 0, ifo, source)
        retval = (n, noises, sss)

    # --------------------------------------------------------
    # output graphics

    if fig:
        # Report input parameters
        if ifo.Optics.Type == 'DualCarrier_new':     #include the case for Dual carrier
            finesseA = 2*pi/ifo.Optics.ITM.TransmittanceD1
            finesseB = 2*pi/ifo.Optics.ITM.TransmittanceD2 
            pbsA = ifo.Laser.PBSD1
            pbsB = ifo.Laser.PBSD2
            print('Finesse for carrier A:  %7.2f' % finesseA)
            print('Finesse for carrier B:  %7.2f' % finesseB)
            print('Power Recycling Factor: %7.2f' % ifo.PRCgain)
            print('Arm power for carrier A:%7.2f kW' % (finesseA*2/pi*pbsA/2/1000))
            print('Arm power for carrier B:%7.2f kW' % (finesseB*2/pi*pbsB/2/1000))
            print('Power on beam splitter for carrier A: %7.2f W' % pbsA)
            print('Power on beam splitter for carrier B: %7.2f W' % pbsB)
            print('Laser Power for Carrier A:     %7.2f Watt' % ifo.LP1)
            print('Laser Power for Carrier B:     %7.2f Watt' % ifo.LP2)
            print('SRM Detuning for Carrier A:    %7.2f degree' % (ifo.Optics.SRM.TunephaseD1*180/pi))
            print('SRM Detuning for Carrier B:    %7.2f degree' % (ifo.Optics.SRM.TunephaseD2*180/pi))
            print('SRM transmission for Carrier A:%9.4f' % ifo.Optics.SRM.TransmittanceD1)
            print('SRM transmission for Carrier B:%9.4f' % ifo.Optics.SRM.TransmittanceD2)
            print('ITM transmission for Carrier A:%9.4f' % ifo.Optics.ITM.TransmittanceD1)
            print('ITM transmission for Carrier B:%9.4f' % ifo.Optics.ITM.TransmittanceD2)
            print('PRM transmission for both:     %9.4f' % ifo.Optics.PRM.Transmittance)
        else:
            pbs = ifo.Laser.Power * prfactor
            ifo.Laser.ArmPower = finesse*2/pi * pbs/2
            print('Laser Power:            %7.2f Watt' % ifo.Laser.Power)
            print('SRM Detuning:           %7.2f degree' % (ifo.Optics.SRM.Tunephase*180/pi))
            print('SRM transmission:       %9.4f' % ifo.Optics.SRM.Transmittance)
            print('ITM transmission:       %9.4f' % ifo.Optics.ITM.Transmittance)
            print('PRM transmission:       %9.4f' % ifo.Optics.PRM.Transmittance)
            print('Finesse:                %7.2f' % finesse)
            print('Power Recycling Gain:   %7.2f' % prfactor)
            print('Arm Power:              %7.2f kW' % (ifo.Laser.ArmPower/1000))
            print('Power on BS:            %7.2f W' % pbs)

        # coating and substrate thermal load on the ITM
        PowAbsITM = (pbs/2) * \
                    np.hstack([(finesse*2/pi) * ifo.Optics.ITM.CoatingAbsorption,
                               (2 * ifo.Materials.MassThickness) * ifo.Optics.ITM.SubstrateAbsorption])

        # Stefan's Mysterious TCS Section ~~~~ `~~ ` ` 324@#%@#$ !    
        M = np.array([[ifo.TCS.s_cc, ifo.TCS.s_cs], [ifo.TCS.s_cs, ifo.TCS.s_ss]])
        S_uncorr = PowAbsITM.T*M*PowAbsITM
        TCSeff = 1-sqrt(ifo.TCS.SRCloss/S_uncorr)

        print('Thermal load on ITM:    %8.3f W' % sum(PowAbsITM))
        print('Thermal load on BS:     %8.3f W' %
              (ifo.Materials.MassThickness*ifo.Optics.SubstrateAbsorption*pbs))
        #fprintf(['Required TCS efficiency: %8.3f' ...
        #              '(estimate, see IFOModel.m for definition)\n'],    TCSeff);  
        if (ifo.Laser.Power*prfactor != pbs):
            print('Lensing limited input power: %7.2f W' % (pbs/prfactor))

        if source:
            print('BNS Inspiral Range:     ' + str(sss.effr0ns) + ' Mpc/ z = ' + str(sss.zHorizonNS))
            print('BBH Inspiral Range:     ' + str(sss.effr0bh) + ' Mpc/ z = ' + str(sss.zHorizonBH))
            print('Stochastic Omega: %4.1g Universes' % sss.Omega)

        plot.plot_noise(noises)

    return retval
