# Ketju-integrator Python bindings
# Copyright (C) 2020-2022 Matias Mannerkoski and contributors.
# Licensed under the GPLv3 license.
# See LICENSE.md or <https://www.gnu.org/licenses/> for details.

"""
Bindings to the BH merger remnant calculation for testing.
"""
import os
import ctypes

import numpy as np

__all__ = ["normalized_merger_properties"]

class BHMergerInput(ctypes.Structure):
    _fields_ = [
    ('G', ctypes.c_double),
    ('c', ctypes.c_double),
    ('m1', ctypes.c_double),
    ('m2', ctypes.c_double),
    ('x1', ctypes.c_double*3),
    ('x2', ctypes.c_double*3),
    ('v1', ctypes.c_double*3),
    ('v2', ctypes.c_double*3),
    ('S1', ctypes.c_double*3),
    ('S2', ctypes.c_double*3),
    ]

    def __init__(self, m1, m2,  x1, x2, v1, v2, S1, S2, G=1., c=1.):
        self.m1 = m1
        self.m2 = m2
        self.G = G
        self.c = c
        np.ctypeslib.as_array(self.x1)[:] = np.asarray(x1)
        np.ctypeslib.as_array(self.x2)[:] = np.asarray(x2)
        np.ctypeslib.as_array(self.v1)[:] = np.asarray(v1)
        np.ctypeslib.as_array(self.v2)[:] = np.asarray(v2)
        np.ctypeslib.as_array(self.S1)[:] = np.asarray(S1)
        np.ctypeslib.as_array(self.S2)[:] = np.asarray(S2)


class BHMergerRemnantProperties(ctypes.Structure):
    _fields_ = [
    ('m', ctypes.c_double),
    ('vkick', ctypes.c_double*3),
    ('S', ctypes.c_double*3),
    ]

def _load_shim():
    _this_dir = os.path.dirname(os.path.realpath(__file__))
    shim_file = os.path.join(_this_dir, "libketju_integrator_python_shim.so")
    shim = ctypes.cdll.LoadLibrary(shim_file)

    shim.ketju_calculate_bh_merger_remnant_properties.argtypes = [
                ctypes.POINTER(BHMergerInput),
                ctypes.POINTER(BHMergerRemnantProperties)]
    return shim
shim = _load_shim()


def normalized_merger_properties(q, R, V, chi1_vec, chi2_vec):
    """
    Wrapper for the merger remnant calculation function used by the integrator.

    Uses units with G=c=1, with the binary total mass set to 1.

    Arguments
    ---------
    q: mass ratio of the binary.
    R: separation of the binary along the x-axis.
    V: relative velocity of the binary along the y-axis.
    chi1_vec, chi2_vec: dimensionless spin vectors of the BHs (chi = (spin angular momentum)/m**2).

    Returns
    -------
    A dict containing the merger remnant properties, with keys:
        mass: mass of the remnant.
        vkick: recoil kick velocity.
        S: spin angular momentum vector.
        chi: dimensionless spin magnitude.
    """

    m1 = 1/(1+q)
    m2 = 1 - m1
    S1 = m1**2 * np.asarray(chi1_vec)
    S2 = m2**2 * np.asarray(chi2_vec)
    bh_in = BHMergerInput(m1, m2, [m2*R, 0, 0], [-m1*R, 0, 0], [0, m2*V, 0],
                          [0, -m1*V, 0], S1, S2)
    bh_out = BHMergerRemnantProperties()
    shim.ketju_calculate_bh_merger_remnant_properties(bh_in, bh_out)

    return dict(mass=bh_out.m,
                vkick=np.array(bh_out.vkick),
                S=np.array(bh_out.S),
                chi=np.linalg.norm(bh_out.S)/bh_out.m**2)

