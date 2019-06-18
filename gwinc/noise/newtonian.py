from __future__ import division
import scipy.constants
from numpy import pi, sqrt, exp


def gravg(f, ifo):
    """Return estimate of newtonian noise contribribution to |h(f)|^2

    N = GRAVG(F, IFO) returns the gravity gradient (Newtonian) noise
    contribution in strain^2/Hz for four mirrors.

    References:

     Saulson 1984,           http://dx.doi.org/10.1103/PhysRevD.30.732
     Hughes and Thorne 1998, http://dx.doi.org/10.1103/PhysRevD.58.122002

     Driggers and Harms 2011, ``Results of Phase 1 Newtonian Noise
     Measurements at the LIGO Sites,'' February-March 2011.  T1100237.
     https://dcc.ligo.org/cgi-bin/private/DocDB/ShowDocument?docid=60064

    Written by Enrico Camagna (?)

    added to Bench by Gregg Harry 8/27/03
    seismic spectrum modified by Jan Harms 05/11/2010
    Calculates gravity gradient noise for four mirrors

    """

    fk = ifo.Seismic.KneeFrequency
    a = ifo.Seismic.LowFrequencyLevel
    L = ifo.Infrastructure.Length
    gamma = ifo.Seismic.Gamma
    ggcst = scipy.constants.G
    rho = ifo.Seismic.Rho
    # factor to account for correlation between masses
    # and the height of the mirror above the ground
    beta = ifo.Seismic.Beta
    h = ifo.Seismic.TestMassHeight
    c_rayleigh = ifo.Seismic.RayleighWaveSpeed

    if 'Omicron' in ifo.Seismic:
        omicron = ifo.Seismic.Omicron
    else:
        omicron = 1

    # a sort of theta function (Fermi distr.)
    coeff = 3**(-gamma*f)/(3**(-gamma*f) + 3**(-gamma*fk))

    # modelization of seismic noise (vertical)
    ground = a*coeff + a*(1-coeff)*(fk/f)**2
    if 'Site' in ifo.Seismic and ifo.Seismic.Site == 'LLO':
        ground = a*coeff*(fk/f) + a*(1-coeff)*(fk/f)**2

    # effective GG spring frequency, with G gravitational
    fgg = sqrt(ggcst * rho) / (2*pi)

    # fixed numerical factors, 5/9/06, PF
    n = (beta*4*pi/L*(fgg**2/f**2)*ground)**2

    # The following two terms are corrections due to Jan Harms
    # https://git.ligo.org/rana-adhikari/CryogenicLIGO/issues/45
    # (1) projection of NN force onto the direction of the arm
    n = n * 1/2
    # (2) exponential cutoff at frequency (seismic speed)/(test mass height)
    n = n * exp(-4*pi*f*h/c_rayleigh)

    # Feedforward cancellation
    n /= (omicron**2)

    return n * ifo.gwinc.sinc_sqr

def atmois(f, ifo):
    import scipy.special as scisp

    a_if = ifo.Atmospheric.InfrasoundLevel1Hz
    e_if = ifo.Atmospheric.InfrasoundExponent
    p_air = ifo.Atmospheric.AirPressure
    rho_air = ifo.Atmospheric.AirDensity
    ai_air = ifo.Atmospheric.AdiabaticIndex
    c_sound = ifo.Atmospheric.SoundSpeed

    L = ifo.Infrastructure.Length
    ggcst = scipy.constants.G
    h = ifo.Seismic.TestMassHeight

    w = 2 * pi * f
    k = w / c_sound

    # Pressure spectrum
    psd_if = (a_if * f**e_if)**2

    # Harms LRR (2015), eq. 172
    # https://doi.org/10.1007/lrr-2015-3
    # With an extra factor 2 for two arms
    # And with the Bessel terms ignored... for 4 km this amounts to a 10%
    # correction at 10 Hz and a 30% correction at 1 Hz
    coupling_if = 4./3 * (4 * pi / (k * L * w**2) * ggcst * rho_air / (ai_air * p_air))**2

    n_if = coupling_if * psd_if

    return n_if * ifo.gwinc.sinc_sqr
