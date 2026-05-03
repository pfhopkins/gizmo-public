# Ketju-integrator Python bindings
# Copyright (C) 2020-2023 Matias Mannerkoski and contributors.
# Licensed under the GPLv3 license.
# See LICENSE.md or <https://www.gnu.org/licenses/> for details.

"""
An example of using the bindings to the merger remnant properties model used
in the integrator.
Produces a distribution of remnant properties from random spin directions.
Uses the ketju_utils package for unit conversions.
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import uniform_direction

from ketju_integrator import normalized_merger_properties
import ketju_utils

# Binary parameters, the binary total mass is 1
# mass ratio
q = .5
# spin magnitudes
chi1 = .8
chi2 = .1

# Use Schwarzschild ISCO values for the orbit for simplicity, this affects the
# results a bit compared to more accurate values
R = 6
V = 1/np.sqrt(3)


# Generate random spin directions
dist = uniform_direction(3)
N=50000
chi1_dirs = dist.rvs(N)
chi2_dirs = dist.rvs(N)
kick_velocity = []
mass_loss = []
chi = []

for d1, d2 in zip(chi1_dirs, chi2_dirs):
    props = normalized_merger_properties(q, R, V, chi1*d1, chi2*d2)
    mass_loss.append(1-props['mass'])
    kick_velocity.append(np.linalg.norm(props['vkick'])/ketju_utils.units.km_per_s)
    chi.append(props['chi'])

fig, axes = plt.subplots(1,3)

axes[0].hist(kick_velocity, bins=100, histtype='step')
axes[0].set_xlabel(r"Kick velocity $v_\mathrm{kick}/\mathrm{km\,s^{-1}}$")

axes[1].hist(mass_loss, bins=100, histtype='step')
axes[1].set_xlabel("Relative mass loss")

axes[2].hist(chi, bins=100, histtype='step')
axes[2].set_xlabel(r"Dimensionless spin magnitude $\chi$")

plt.show()
