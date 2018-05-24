from __future__ import division
from numpy import pi, imag, zeros
import numpy as np
import scipy.constants


def susptherm(f, ifo):
    """Suspention thermal noise.

    Assumes suspension transfer functions and V-H coupling have been
    pre-calculated and populated into the relevant `ifo` struct
    fields.

    """
    # Assign Physical Constants
    kB   = scipy.constants.k
    Temp = ifo.Suspension.Temp

    # and vertical to beamline coupling angle
    theta = ifo.Suspension.VHCoupling.theta

    noise = zeros((1, f.size))

    # if the temperature is uniform along the suspension
    if np.isscalar(Temp) or len(Temp) == 1:
        ##########################################################
        # Suspension TFs
        ##########################################################
        hForce = ifo.Suspension.hForce
        vForce = ifo.Suspension.vForce

        ##########################################################
        # Thermal Noise Calculation
        ##########################################################

        # convert to beam line motion
        #  theta is squared because we rotate by theta into the suspension
        #  basis, and by theta to rotate back to the beam line basis
        dxdF = hForce + theta**2 * vForce

        # thermal noise (m^2/Hz) for one suspension
        w = 2*pi*f
        noise = 4 * kB * Temp * abs(imag(dxdF)) / w

    # if the temperature is set for each suspension stage
    else:
        ##########################################################
        # Suspension TFs
        ##########################################################
        hForce = ifo.Suspension.hForce_singlylossy[:,:]
        vForce = ifo.Suspension.vForce_singlylossy[:,:]

        ##########################################################
        # Thermal Noise Calculation
        ##########################################################

        dxdF = zeros(hForce.shape, dtype=complex)
        for ii in range(len(Temp)):
            # add up the contribution from each stage

            # convert to beam line motion
            #  theta is squared because we rotate by theta into the suspension
            #  basis, and by theta to rotate back to the beam line basis
            dxdF[ii,:] = hForce[ii,:] + theta**2 * vForce[ii,:]

            # thermal noise (m^2/Hz) for one suspension
            w = 2*pi*f
            noise += 4 * kB * Temp[ii] * abs(imag(dxdF[ii,:])) / w

    # 4 masses, turn into gravitational wave strain
    noise *= 4 * ifo.gwinc.dhdl_sqr

    return np.squeeze(noise)
