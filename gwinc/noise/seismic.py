from __future__ import division, print_function
import numpy as np
from numpy import log10
from scipy.interpolate import interp1d


def seismic(f, ifo):
    """seismic noise psd at frequencies f for given ifo.
      n = seismic(f,ifo)
      [nh,nv] = seismic(f,ifo)
      [n,nh,nv] = seismic(f,ifo)
    
    Modified to include realistic SEI + SUS models (Rana, 11/2005)"""

    # Interpolate the log10 onto the ifo frequency array
    #   n = interp1(ifo.Seismic.darmseis_f, ...
    #     log10(ifo.Seismic.darmseis_x), f, 'cubic', -30);

    ##########################################################
    # Suspension TFs
    ##########################################################
    hTable = ifo.Suspension.hTable
    vTable = ifo.Suspension.vTable

    # and vertical to beamline coupling angle
    theta = ifo.Suspension.VHCoupling.theta

    ##########################################################

    # noise input, horizontal and vertical
    nx, np_ = seisBSC(f)

    # horizontal noise total
    nh = (abs(hTable)**2) * nx**2

    # vertical noise total
    nv = (abs(theta * vTable)**2) * nx**2

    # new total noise
    n = nv + nh

    # Convert into Strain PSD (4 TMs)
    nh *= 4 / ifo.Infrastructure.Length**2
    nv *= 4 / ifo.Infrastructure.Length**2
    n  *= 4 / ifo.Infrastructure.Length**2

    return n


def seisBSC(f):
    """get rough ISI noise source spectra
    nx - ISI translational DOFs
    np - ISI rotational DOFs"""


    # translational DOFs (from Rana's bsc_seismic.m)
    SEI_F = np.array([0.01, 0.03, 0.1, 0.2, 0.5, 1, 10, 30, 300])
    SEI_X = np.array([3e-6, 1e-6, 2e-7, 2e-7, 8e-10, 1e-11, 3e-13, 3e-14, 3e-14])
    nx = 10**(interp1d(SEI_F,log10(SEI_X), bounds_error=False, fill_value=-14)(f))
  
    # rotational DOFs
    SEI_P = np.array([1e-8, 3e-8, 2e-8, 1e-8, 4e-10, 1e-11, 3e-13, 3e-14, 3e-14])
    np_ = 10**(interp1d(SEI_F,log10(SEI_P), bounds_error=False, fill_value=-14)(f))
  
    return nx, np_
