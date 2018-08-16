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
