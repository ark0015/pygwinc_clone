from __future__ import division, print_function
import scipy.constants
from numpy import pi, sqrt

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
    
    See also: ground.m
    
    Written by Enrico Camagna (?)
    
    added to Bench by Gregg Harry 8/27/03
    seismic spectrum modified by Jan Harms 05/11/2010
    Calculates gravity gradient noise for four mirrors
    ** --- Add reference here -- **"""

    L     = ifo.Infrastructure.Length
    G     = scipy.constants.G
    rho   = ifo.Seismic.Rho
    beta  = ifo.Seismic.Beta        # factor to account for correlation between masses
                                    # and the height of the mirror above the ground

    if 'Spectrum' in ifo.Seismic:
        seismic = ifo.Seismic.Spectrum
    else:
        seismic = ground(ifo.Seismic, f)

    # effective GG spring frequency, with G gravitational
    fgg = sqrt(G * rho) / (2*pi)

    # fixed numerical factors, 5/9/06, PF
    n = (beta*4*pi/L*(fgg**2/f**2)*seismic)**2

    # Feedforward cancellation
    n = n / (ifo.Seismic.Omicron)**2
    return n


def ground(Seismic, f):
    """Returns estimate of seismic displacement spectrum
    N = GROUND(SEISMIC, F) returns the estimated ground (seismic)
    displacement spectrum in meters / rtHz at frequencies F using
    parameters SEISMIC (usually a field of the IFOMODEL).
    
    Example usage:
    
        ifo = IFOModel();
        f = logspace(log10(1), log10(60), 101);
        n = ground(ifo.Seismic, f);
        loglog(f, n)
        ylabel('meters / rtHz');
        xlabel('Hz');
    
      This function is currently only used by gravg.m (Newtonian Noise).
    
      Jan says: The seismic NN is 90th percentile. But this was just to
      lift the spectra above the aLIGO GWINC sensitivity. 90th, 95th,
      whatever percentile does not matter since you want to guarantee that
      NN is not limiting at any time. So in the spectral plot the only
      important information is that NN can be high enough to be seen in
      the detector.
    
      References:
        Waldman and Fritschel, 'Reference Seismic Data for LLO', T0900312
        https://dcc.ligo.org/cgi-bin/private/DocDB/ShowDocument?docid=3315"""

    fk = Seismic.KneeFrequency
    a1 = Seismic.LowFrequencyLevel
    a2 = a1*100
    gamma = Seismic.Gamma

    # a sort of theta function (Fermi distr.)
    coeff = 1/(1 + 3**(gamma*(f-fk)))

    # modelization of seismic noise (velocity)
    if 'Site' not in Seismic.__dict__:
        print('defaulting to Livingston site')
        Seismic.Site = 'LLO'

    if Seismic.Site == 'LHO':
        n = (2*pi*f)**(4/3)*(a1*coeff + a1*(1-coeff)*(fk/f)**(9/3))
    elif Seismic.Site == 'LLO':
        n = a2*coeff + a2*(1-coeff)*(fk/f)**2
    else:
        raise 'Unknown seismic.site'

    # convert into displacement
    n = n/(2*pi*f)
    return n
