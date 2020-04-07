'''Functions to calculate suspension thermal noise

'''
from __future__ import division
import numpy as np
from numpy import pi, imag

from .. import const


def suspension_thermal(f, sus):
    """Suspention thermal noise for a single suspended test mass

    :f: frequency array in Hz
    :sus: gwinc suspension structure

    :returns: displacement noise power spectrum at :f:

    :Temp: must either be set for each stage individually, or globally
    in :sus.Temp:.  The latter will be preferred if specified, so if
    you wish to use per-stage tempurature you must remove :sus.Temp:.

    Assumes suspension transfer functions and V-H coupling have been
    pre-calculated and populated into the relevant `sus` fields.

    """
    # Assign Physical Constants
    kB = const.kB

    # and vertical to beamline coupling angle
    theta = sus.VHCoupling.theta

    noise = np.zeros((1, f.size))

    # if the temperature is uniform along the suspension
    if 'Temp' in sus:
        ##########################################################
        # Suspension TFs
        ##########################################################

        hForce = sus.hForce
        vForce = sus.vForce

        ##########################################################
        # Thermal Noise Calculation
        ##########################################################

        # convert to beam line motion
        #  theta is squared because we rotate by theta into the suspension
        #  basis, and by theta to rotate back to the beam line basis
        dxdF = hForce + theta**2 * vForce

        # thermal noise (m^2/Hz) for one suspension
        w = 2*pi*f
        noise = 4 * kB * sus.Temp * abs(imag(dxdF)) / w

    # if the temperature is set for each suspension stage
    else:
        ##########################################################
        # Suspension TFs
        ##########################################################

        hForce = sus.hForce_singlylossy[:, :]
        vForce = sus.vForce_singlylossy[:, :]

        ##########################################################
        # Thermal Noise Calculation
        ##########################################################

        dxdF = np.zeros(hForce.shape, dtype=complex)
        for n, stage in enumerate(reversed(sus.Stage)):
            # add up the contribution from each stage

            # convert to beam line motion.  theta is squared because
            # we rotate by theta into the suspension basis, and by
            # theta to rotate back to the beam line basis
            dxdF[n, :] = hForce[n, :] + theta**2 * vForce[n, :]

            # thermal noise (m^2/Hz) for one suspension
            w = 2*pi*f
            noise += 4 * kB * stage.Temp * abs(imag(dxdF[n, :])) / w

    return np.squeeze(noise)
