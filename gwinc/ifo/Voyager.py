from __future__ import division, print_function
from numpy import pi, NaN
from ..util import SpotSizes
import scipy.constants
import scipy.special
from scipy.io import loadmat
from scipy.io.matlab.mio5_params import mat_struct
import os


def IFOModel():
    """IFOMODEL returns a structure describing an IFO for use in
    benchmark programs and noise simulator. Part of the gwinc
    package, which provides science-grounded figures of merit for
    comparing interferometric gravitational wave detector designs."""

    ifo = mat_struct()

    ## Infrastructure----------------------------------------------------------

    ifo.Infrastructure = mat_struct()
    ifo.Infrastructure.Length                     = 3995      # m
    ifo.Infrastructure.ResidualGas = mat_struct()
    ifo.Infrastructure.ResidualGas.pressure       = 4.0e-7    # Pa
    ifo.Infrastructure.ResidualGas.mass           = 3.35e-27  # kg,   Mass of H_2 (ref. 10)
    ifo.Infrastructure.ResidualGas.polarizability = 0.81e-30  # m^3  (H_2, DOI: 10.1116/1.1479360)

    ## Physical and other constantMaterialss; All in SI units------------------

    ifo.Constants = mat_struct()
    #New version of constants
    #ifo.Constants.E0      = 8.8541878176e-12;                 % F / m; Permittivity of Free Space
    #ifo.Constants.hbar    = 1.054572e-34;                     % J - s; (Plancks constant) / (2 * pi)
    #ifo.Constants.kB      = 1.380658e-23;                     % J / K; Boltzman Constant
    #ifo.Constants.h       = ifo.Constants.hbar * 2 * pi;      % J - s; Planks constant
    #ifo.Constants.R       = 8.31447215;                       % J / (K * mol); Gas Constant
    #ifo.Constants.m_e     = 9.10938291e-31;                   % kg; electron mass
    #ifo.Constants.c       = 2.99792458e8;                     % m / s; speed of light in vacuum
    ifo.Constants.Temp    = 295                              # K; Temperature of the Vacuum
    #ifo.Constants.yr      = 365.2422 * 86400;                 % sec; Seconds in a year
    #ifo.Constants.M_earth = 5.972e24;                         % mass of Earth [kg]
    #ifo.Constants.R_earth = 6.3781e6;                         % radius of Earth [m]
    #ifo.Constants.fs      = 16384;                            % Sampling frequency (Hz)
    #ifo.Constants.AU      = 149597870700;                     % m; Astronomical unit, IAU 2012 Resolution B2
    #ifo.Constants.parsec  = ifo.Constants.AU * (648000 / pi); % m, IAU 2015 Resolution B2
    #ifo.Constants.Mpc     = ifo.Constants.parsec * 1e6;       % m, IAU 2015 Resolution B2
    #ifo.Constants.SolarMassParameter = 1.3271244e20;          % m^3 / s^2; G * MSol, IAU 2015 Resolution B3
    #ifo.Constants.G       = 6.67408e-11;                      % m^3 / (kg  s^2); Grav. const
    #                                                          % http://arxiv.org/abs/1507.07956
    #ifo.Constants.MSol    = ifo.Constants.SolarMassParameter / ifo.Constants.G; % kg; Solar mass
    #ifo.Constants.g       = 9.806;                            % m / s^2; grav. acceleration 
    #                                                          % http://physics.nist.gov/cuu/Constants/ 
    #ifo.Constants.H0      = 67110;                            % ms^( - 1); Hubble const.
    #                                                          % http://arxiv.org/pdf/1303.5076v3.pdf
    #ifo.Constants.omegaM  = 0.3175;                           % Mass density parameter 
    #                                                          % http://arxiv.org/pdf/1303.5076v3.pdf
    #ifo.Constants.omegaLambda = 1 - ifo.Constants.omegaM;     % Cosmological constant density parameter
    #                                                          % omegaK = 0 (flat universe) is assumed

  
    #ifo.Constants.fInspiralMin = 3;  % cut-off for inspiral range (Hz, see int73)
  
    ## Parameter describing thermal lensing --------------------------------------
    # The presumably dominant effect of a thermal lens in the ITMs is an increased
    # mode mismatch into the SRC, and thus an increased effective loss of the SRC.
    # This increase is estimated by calculating the round-trip loss S in the SRC as
    # 1-S = |<Psi|exp(i*phi)|Psi>|^2, where
    # |Psi> is the beam hitting the ITM and
    # phi = P_coat*phi_coat + P_subs*phi_subs
    # with phi_coat & phi__subs the specific lensing profiles
    # and P_coat & P_subst the power absorbed in coating and substrate
    #
    # This expression can be expanded to 2nd order and is given by
    # S= s_cc P_coat^2 + 2*s_cs*P_coat*P_subst + s_ss*P_subst^2
    # s_cc, s_cs and s_ss where calculated analytically by Phil Wilems (4/2007)
    ifo.TCS = mat_struct()
    ifo.TCS.s_cc=7.024 # Watt^-2
    ifo.TCS.s_cs=7.321 # Watt^-2
    ifo.TCS.s_ss=7.631 # Watt^-2
  
    # The hardest part to model is how efficient the TCS system is in
    # compensating this loss. Thus as a simple Ansatz we define the
    # TCS efficiency TCSeff as the reduction in effective power that produces
    # a phase distortion. E.g. TCSeff=0.99 means that the compensated distortion
    # of 1 Watt absorbed is eqivalent to the uncompensated distortion of 10mWatt.
    # The above formula thus becomes:
    # S= s_cc P_coat^2 + 2*s_cs*P_coat*P_subst + s_ss*P_subst^2 * (1-TCSeff)^2
    #
    # To avoid iterative calculation we define TCS.SCRloss = S as an input
    # and calculate TCSeff as an output.
    # TCS.SRCloss is incorporated as an additional loss in the SRC
    ifo.TCS.SRCloss = 0.00


    ## Seismic and Gravity Gradient Parameters---------------------------------
    ifo.Seismic = mat_struct()
    ifo.Seismic.Site = 'LHO'                      # LHO or LLO (only used for Newtonian noise)
    ifo.Seismic.KneeFrequency = 10                # Hz; freq where 'flat' noise rolls off
    ifo.Seismic.LowFrequencyLevel = 1e-9          # m/rtHz; seismic noise level below f_knee
    ifo.Seismic.Gamma = .8                        # abruptness of change at f_knee
    ifo.Seismic.Rho = 1.8e3                       # kg/m^3; density of the ground nearby
    ifo.Seismic.Beta = 1                          # quiet times beta = 0.35-0.60
                                                  # noisy times beta = 0.15-1.4
    ifo.Seismic.Omicron = 10                      # Feedforward cancellation factor
    ifo.Seismic.darmSeiSusFile = 'CryogenicLIGO/Sensitivity/GWINC/' + 'seismic.mat'

  
    ## Suspension: SI Units----------------------------------------------------
    ifo.Suspension = mat_struct()
    ifo.Suspension.BreakStress      = 750e6           # Pa; ref. K. Strain
    # Suspension fiber temperatures  [TOP UIM PUM UTM]
    ifo.Suspension.Temp             = [300, 300, 300, 123]
    ifo.Suspension.VHCoupling = mat_struct()
    ifo.Suspension.VHCoupling.theta = 1e-3        # vertical-horizontal x-coupling
  
    # new Silicon parameters added for gwincDev   RA  April, 20, 2010 ~~~~~~~~~~~~~~~~~~~
    # new Silicon parameters added for gwincDev   RA  Feb 25, 2012 ~~~~~~~~~~~~~~~~~~~
    # http://www.ioffe.ru/SVA/NSM/Semicond/Si/index.html
    # all properties should be for T ~ 120 K
    ifo.Suspension.Silicon = mat_struct()
    ifo.Suspension.Silicon.Rho       = 2329          # Kg/m^3   density
    ifo.Suspension.Silicon.C         = 300           # J/kg/K   heat capacity
    ifo.Suspension.Silicon.K         = 700           # W/m/K    thermal conductivity
    ifo.Suspension.Silicon.Alpha     = 1e-10         # 1/K      thermal expansion coeff

    # from Gysin, et. al. PRB (2004)  E(T) = E0 - B*T*exp(-T0/T)
    # E0 = 167.5e9 Pa   T0 = 317 K   B = 15.8e6 Pa/K
    ifo.Suspension.Silicon.dlnEdT    = -2e-5         # (1/K)    dlnE/dT  T = 120K

    ifo.Suspension.Silicon.Phi       = 2e-9          # Nawrodt (2010)      loss angle  1/Q
    ifo.Suspension.Silicon.Y         = 155.8e9       # Pa       Youngs Modulus
    ifo.Suspension.Silicon.Dissdepth = 1.5e-3        # 10x smaller surface loss depth (Nawrodt (2010))
    ifo.Suspension.FiberType         = 1             # 0 = round, 1 = ribbons
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ifo.Suspension.Silica = mat_struct()
    ifo.Suspension.Silica.Rho    = 2200           # Kg/m^3
    ifo.Suspension.Silica.C      = 772            # J/Kg/K
    ifo.Suspension.Silica.K      = 1.38           # W/m/kg
    ifo.Suspension.Silica.Alpha  = 3.9e-7         # 1/K
    ifo.Suspension.Silica.dlnEdT = 1.52e-4        # (1/K), dlnE/dT
    ifo.Suspension.Silica.Phi    = 4.1e-10        # from G Harry e-mail to NAR 27April06
    ifo.Suspension.Silica.Y      = 72e9           # Pa; Youngs Modulus
    ifo.Suspension.Silica.Dissdepth = 1.5e-2      # from G Harry e-mail to NAR 27April06

    ifo.Suspension.C70Steel = mat_struct()
    ifo.Suspension.C70Steel.Rho    =  7800
    ifo.Suspension.C70Steel.C      =  486
    ifo.Suspension.C70Steel.K      =  49
    ifo.Suspension.C70Steel.Alpha  =  12e-6
    ifo.Suspension.C70Steel.dlnEdT = -2.5e-4
    ifo.Suspension.C70Steel.Phi    =  2e-4
    ifo.Suspension.C70Steel.Y      = 212e9        # measured by MB for one set of wires

    ifo.Suspension.MaragingSteel = mat_struct()
    ifo.Suspension.MaragingSteel.Rho = 7800
    ifo.Suspension.MaragingSteel.C   = 460
    ifo.Suspension.MaragingSteel.K   = 20
    ifo.Suspension.MaragingSteel.Alpha  = 11e-6
    ifo.Suspension.MaragingSteel.dlnEdT = 0
    ifo.Suspension.MaragingSteel.Phi  = 1e-4
    ifo.Suspension.MaragingSteel.Y  = 187e9
    # consistent with measured blade spring constants NAR

    ifo.Suspension.Type         = 'BQuad'               # 0 for cylindrical suspension

    # Note stage numbering: mirror is at beginning of stack, not end
    # these mass numbers are from v8 of the Voyager design doc
    ifo.Suspension.Stage1 = mat_struct()
    ifo.Suspension.Stage2 = mat_struct()
    ifo.Suspension.Stage3 = mat_struct()
    ifo.Suspension.Stage4 = mat_struct()

    ###addpath('../../QuadModel/')        # add path of saved file with optimized masses
    ###load(quad_optimized_masses_for_PUM_with_springs) # Load saved file with otpimized mass. Masses are optimized for longitudinal isolation assuming the PUM has springs
    susmat = loadmat('CryogenicLIGO/QuadModel/quad_optimized_masses_for_PUM_with_springs.mat')
    ifo.Suspension.Stage1.Mass = susmat['testmass_mass'][0,0]   # kg; this is redefined below for some reason
    ifo.Suspension.Stage2.Mass = susmat['PUMmass'][0,0]
    ifo.Suspension.Stage3.Mass = susmat['UIMmass'][0,0]
    ifo.Suspension.Stage4.Mass = susmat['topmass_mass'][0,0]

    ifo.Suspension.Stage1.Length = 0.4105        # m
    ifo.Suspension.Stage2.Length = 0.4105        # m
    ifo.Suspension.Stage3.Length = 0.4105        # m
    ifo.Suspension.Stage4.Length = 0.4105        # m

    ifo.Suspension.Stage1.Dilution = NaN
    ifo.Suspension.Stage2.Dilution = 106
    ifo.Suspension.Stage3.Dilution = 80
    ifo.Suspension.Stage4.Dilution = 87

    ifo.Suspension.Stage1.K = NaN
    ifo.Suspension.Stage2.K = 5200              # N/m; vertical spring constant
    ifo.Suspension.Stage3.K = 3900              # N/m; vertical spring constant
    ifo.Suspension.Stage4.K = 3400              # N/m; vertical spring constant

    ifo.Suspension.Stage1.WireRadius = NaN
    ifo.Suspension.Stage2.WireRadius = 310e-6
    ifo.Suspension.Stage3.WireRadius = 350e-6
    ifo.Suspension.Stage4.WireRadius = 520e-6

    # For Ribbon suspension
    ifo.Suspension.Ribbon = mat_struct()
    ifo.Suspension.Fiber = mat_struct()
    ifo.Suspension.Ribbon.Thickness = 115e-6      # m
    ifo.Suspension.Ribbon.Width     = 1150e-6     # m
    ifo.Suspension.Fiber.Radius     = 205e-6      # m

    ifo.Suspension.Stage1.Blade = NaN            # blade thickness
    ifo.Suspension.Stage2.Blade = 4200e-6
    ifo.Suspension.Stage3.Blade = 4600e-6
    ifo.Suspension.Stage4.Blade = 4300e-6

    ifo.Suspension.Stage1.NWires = 4
    ifo.Suspension.Stage2.NWires = 4
    ifo.Suspension.Stage3.NWires = 4
    ifo.Suspension.Stage4.NWires = 2


    ## Amorphous Silicon / Silica coating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  high index material: a-Si    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  https://wiki.ligo.org/OPT/AmorphousSilicon
    ifo.Materials = mat_struct()
    ifo.Materials.Coating = mat_struct()
    ifo.Materials.Coating.Yhighn     = 80e9
    ifo.Materials.Coating.Sigmahighn = 0.22
    ifo.Materials.Coating.CVhighn    = 345.6*2250       # volume-specific heat capacity (J/K/m^3); http://journals.aps.org/prl/pdf/10.1103/PhysRevLett.96.055902
    ifo.Materials.Coating.Alphahighn = 1e-9             # zero crossing at 123 K
    ifo.Materials.Coating.Betahighn  = 1.4e-4           # dn/dT
    ifo.Materials.Coating.ThermalDiffusivityhighn = 1   # W/m/K (this is a misnomer, meant to be thermal conductivity not diffusivity)
    ifo.Materials.Coating.Phihighn   = 3e-5             # just a guess (depends on prep)
    ifo.Materials.Coating.Indexhighn = 3.5

    ## low index material: silica
    #  https://wiki.ligo.org/OPT/SilicaCoatingProp
    ifo.Materials.Coating.Ylown      = 72e9              # Young's modulus (Pa)
    ifo.Materials.Coating.Sigmalown  = 0.17              # Poisson's ratio
    ifo.Materials.Coating.CVlown     = 1.6412e6          # volume-specific heat capacity (J/K/m^3); Crooks et al, Fejer et al
    ifo.Materials.Coating.Alphalown  = 5.1e-7            # Fejer et al
    ifo.Materials.Coating.Betalown   = 8e-6              # dn/dT,  (ref. 14)
    ifo.Materials.Coating.ThermalDiffusivitylown = 1.38  # Fejer et al (this is a misnomer, meant to be thermal conductivity not diffusivity)
    ifo.Materials.Coating.Philown    = 1e-4              # ?

    # calculated for 123 K and 2000 nm following 
    # Ghosh, et al (1994):  http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=317500
    ifo.Materials.Coating.Indexlown  = 1.436             # calculated (RXA)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ## Substrate Material parameters--------------------------------------------
    # Silicon @ 120K (http://www.ioffe.ru/SVA/NSM/Semicond/Si/index.html)

    ifo.Materials.Substrate = mat_struct()
                                                          #  phi_sub = c2 * f^(MechLossExp)
    ifo.Materials.Substrate.c2                = 3e-13     # Coeff of freq dep. term for bulk loss (Lam & Douglass, 1981)
    ifo.Materials.Substrate.MechanicalLossExponent = 1    # Exponent for freq dependence of silicon loss
    ifo.Materials.Substrate.Alphas            = 5.2e-12   # Surface loss limit ???
    ifo.Materials.Substrate.MirrorY           = 155.8e9   # N/m^2; Youngs modulus (ioffe) -- what about anisotropy??
    ifo.Materials.Substrate.MirrorSigma       = 0.27      # kg/m^3; Poisson ratio (ioffe) -- what about anisotropy??
    ifo.Materials.Substrate.MassDensity       = 2329      # kg/m^3; (ioffe)
    ifo.Materials.Substrate.MassAlpha         = 1e-9      # 1/K; CTE = 0 @ 120 K
    ifo.Materials.Substrate.MassCM            = 0.3*1000  # J/kg/K; specific heat (ioffe @ 120K)
    ifo.Materials.Substrate.MassKappa         = 700       # W/(m*K); thermal conductivity (ioffe @ 120)
    ifo.Materials.Substrate.RefractiveIndex   = 3.5       # 3.38 * (1 + 4e-5 * T)   (ioffe)
    ifo.Materials.Substrate.dndT              = 1e-4      # ~123K & 1900 nm : http://arxiv.org/abs/physics/0606168


    ifo.Materials.MassRadius    = 0.450/2             # m
    ifo.Materials.MassThickness = 0.55;

    ifo.Materials.Substrate.Temp = 123            # mirror temperature [K]

    ## Laser-------------------------------------------------------------------
    ifo.Laser = mat_struct()
    ifo.Laser.Wavelength                   = 2000e-9;      # m
    ifo.Laser.Power                        = 150          # W                              % W;

    ## Optics------------------------------------------------------------------
    ifo.Optics = mat_struct()
    ifo.Optics.Type = 'SignalRecycled'

    ifo.Optics.ITM = mat_struct()
    ifo.Optics.ETM = mat_struct()
    ifo.Optics.PRM = mat_struct()
    ifo.Optics.SRM = mat_struct()
    ifo.Optics.SRM.CavityLength         = 55      # m; ITM to SRM distance
    ifo.Optics.PhotoDetectorEfficiency  = 0.95    # photo-detector quantum efficiency
    ifo.Optics.Loss                     = 10e-6   # average per mirror power loss
                                                     # factor of 4 for 1064 -> 2000
    ifo.Optics.BSLoss  = 0.5e-3                   # power loss near beamsplitter
    ifo.Optics.coupling = 1.0                  # mismatch btwn arms & SRC modes; used to
                                               # calculate an effective r_srm
    ifo.Optics.Curvature = mat_struct()
    ifo.Optics.Curvature.ITM = 1800               # RoC of ITM
    ifo.Optics.Curvature.ETM = 2500               # RoC of ETM
    ifo.Optics.SubstrateAbsorption = 0.3e-4       # 1/m; 0.3 ppm/cm for Hereaus
    ifo.Optics.ITM.SubstrateAbsorption = 10e-6 / 0.01    # 1/m; 10 ppm/cm for MCZ Si
    ifo.Optics.pcrit = 10                         # W; tolerable heating power (factor 1 ATC)

    # calculate arm cavity spot sizes
    L = ifo.Infrastructure.Length
    w1,w2,junk = SpotSizes(1 - L / ifo.Optics.Curvature.ITM,
                           1 - L / ifo.Optics.Curvature.ETM,
                           L, ifo.Laser.Wavelength)
    ifo.Optics.ITM.BeamRadius = w1                     # m; 1/e^2 power radius
    ifo.Optics.ETM.BeamRadius = w2                     # m; 1/e^2 power radius


    ifo.Optics.ITM.CoatingAbsorption = 1e-6            # absorption of ITM
    ifo.Optics.ITM.Transmittance  = 0.008                # Transmittance of ITM
    ifo.Optics.ETM.Transmittance  = 5e-6                 # Transmittance of ETM
    ifo.Optics.SRM.Transmittance  = 0.16                 # Transmittance of SRM
    ifo.Optics.PRM.Transmittance  = 0.03

    # coating layer optical thicknesses - mevans June 2008
    ifo.Optics.ITM.CoatingThicknessLown = 0.308
    ifo.Optics.ITM.CoatingThicknessCap  = 0.5

    ifo.Optics.ETM.CoatingThicknessLown = 0.27
    ifo.Optics.ETM.CoatingThicknessCap  = 0.5

    #ifo.Optics.SRM.Tunephase = 0.23;           % SRM tuning, 795 Hz narrowband
    ifo.Optics.SRM.Tunephase = 0.0             # SRM tuning [radians]
    ifo.Optics.Quadrature = mat_struct()
    ifo.Optics.Quadrature.dc = pi/2            # homoDyne phase [radians]

    ## Squeezer Parameters------------------------------------------------------

    # Define the squeezing you want:
    #   None = ignore the squeezer settings
    #   Freq Independent = nothing special (no filter cavties)
    #   Freq Dependent = applies the specified filter cavites
    #   Optimal = find the best squeeze angle, assuming no output filtering
    #   OptimalOptimal = optimal squeeze angle, assuming optimal readout phase
    ifo.Squeezer = mat_struct()
    ifo.Squeezer.Type = 'Freq Dependent'
    ifo.Squeezer.AmplitudedB = 10                  # SQZ amplitude [dB]
    ifo.Squeezer.InjectionLoss = 0.05              # power loss to sqz
    ifo.Squeezer.SQZAngle = 0                      # SQZ phase [radians]

    # Parameters for frequency dependent squeezing
    ifo.Squeezer.FilterCavity = mat_struct()
    ifo.Squeezer.FilterCavity.fdetune = -31.5      # detuning [Hz]
    ifo.Squeezer.FilterCavity.L   = 300            # cavity length [m]
    ifo.Squeezer.FilterCavity.Ti  = 800e-6         # input mirror trasmission [Power]
    ifo.Squeezer.FilterCavity.Te  = 0e-6           # end mirror trasmission
    ifo.Squeezer.FilterCavity.Lrt = 10e-6          # round-trip loss in the cavity
    ifo.Squeezer.FilterCavity.Rot = 0 * pi/180     # phase rotation after cavity

    ## Variational Output Parameters--------------------------------------------
    # Define the output filter cavity chain
    #   None = ignore the output filter settings
    #   Chain = apply filter cavity chain
    #   Optimal = find the best readout phase
    ifo.OutputFilter = mat_struct()
    ifo.OutputFilter.Type = 'None'

    ifo.OutputFilter.FilterCavity = mat_struct()
    ifo.OutputFilter.FilterCavity.fdetune = -30   # detuning [Hz]
    ifo.OutputFilter.FilterCavity.L = 4000        # cavity length
    ifo.OutputFilter.FilterCavity_Ti = 10e-3      # input mirror trasmission [Power]
    ifo.OutputFilter.FilterCavity.Te = 0          # end mirror trasmission
    ifo.OutputFilter.FilterCavity.Lrt = 100e-6    # round-trip loss in the cavity
    ifo.OutputFilter.FilterCavity.Rot = 0         # phase rotation after cavity

    ## parameters for semiconductor optics
    ifo.Materials.Substrate.isSemiConductor = True      # we are doing semiconductor optics
    ifo.Materials.Substrate.CarrierDensity = 1e13 * 1e6 # 1/m^3; carrier density for phosphorous-doped silicon
    ifo.Materials.Substrate.ElectronDiffusion = 97 * 1e-4 # m^2/s; electron diffusion coefficient for silicon at 120 K
    ifo.Materials.Substrate.HoleDiffusion = 35 * 1e-4 # m**2/s; hole diffusion coefficient for silicon at 120 K

    ifo.Materials.Substrate.ElectronEffMass = 1.07 * scipy.constants.m_e # kg; effective mass of each electron
    ifo.Materials.Substrate.HoleEffMass = 0.88 * scipy.constants.m_e # kg; effective mass of each hole
    ifo.Materials.Substrate.ElectronIndexGamma = -8.8e-22 * 1e-6 # m**3; dependence of index of refraction on electron carrier density
    ifo.Materials.Substrate.HoleIndexGamma = -10.2e-22 * 1e-6 # m**3; dependence of index of refraction on hole carrier density

    ## Incorporate PSO results--------------------------------------------------
    # quantum - load latest results
    qopt_mat = sorted(os.listdir('CryogenicLIGO/Sensitivity/GWINC/optRuns'))[-1]
    zz = loadmat('CryogenicLIGO/Sensitivity/GWINC/optRuns/' + qopt_mat)

    ifo.Laser.Power                    = zz['x'][0][0]
    ifo.Squeezer.FilterCavity.fdetune  = zz['x'][0][1]
    ifo.Squeezer.FilterCavity.Ti       = zz['x'][0][2]
    ifo.Optics.ITM.Transmittance       = zz['x'][0][3]
    ifo.Optics.SRM.Transmittance       = zz['x'][0][4]
    #ifo.Optics.SRM.Tunephase           = x(6)*0;
    ifo.Optics.Quadrature.dc           = zz['x'][0][5]

    # coating
    itm = loadmat('CryogenicLIGO/Sensitivity/coating/aSi/Data/ITM_layers_151221_2237.mat')
    etm = loadmat('CryogenicLIGO/Sensitivity/coating/aSi/Data/ETM_layers_151221_2150.mat')
    ifo.Optics.ITM.CoatLayerOpticalThickness = itm['TNout']['L'][0][0].T
    ifo.Optics.ETM.CoatLayerOpticalThickness = etm['TNout']['L'][0][0].T

    return ifo


# parameters for quad pendulum suspension updated 3rd May 2006, NAR
# References:
# LIGO-T000012-00-D
# 	* Differentiate between silica and sapphire substrate absorption
# 	* Change ribbon suspension aspect ratio
# 	* Change pendulum frequency
# * References:
# * 1. Electro-Optic Handbook, Waynant & Ediger (McGraw-Hill: 1993)
# * 2. LIGO/GEO data/experience
# * 3. Suspension reference design, LIGO-T000012-00
# * 4. Quartz Glass for Optics Data and Properties, Heraeus data sheet,
# *    numbers for suprasil
# * 5. Y.S. Touloukian (ed), Thermophysical Properties of Matter 
# *    (IFI/Plenum,1970)
# * 6. Marvin J. Weber (ed) CRC Handbook of laser science and technology, 
# *    Vol 4, Pt 2
# * 7. R.S. Krishnan et al.,Thermal Expansion of Crystals, Pergamon Press
# * 8. P. Klocek, Handbook of infrared and optical materials, Marcel Decker, 
# *    1991
# * 9. Rai Weiss, electronic log from 5/10/2006
# * 10. Wikipedia online encyclopedia, 2006
# * 11. D.K. Davies, The Generation and Dissipation of Static Charge on
# * dielectrics in a Vacuum, page 29
# * 12. Gretarsson & Harry, Gretarsson thesis
# * 13. Fejer
# * 14. Braginsky
