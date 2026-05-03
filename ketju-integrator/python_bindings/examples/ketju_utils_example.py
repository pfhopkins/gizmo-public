# Ketju-integrator Python bindings
# Copyright (C) 2020-2023 Matias Mannerkoski and contributors.
# Licensed under the GPLv3 license.
# See LICENSE.md or <https://www.gnu.org/licenses/> for details.

"""
An example of setting up the system and analysing the integration results
using the ketju_utils package (provided in a separate repository, see README).
ketju_utils uses an unit system with G=c=1 and an unit mass of one solar mass internally.
"""

import numpy as np
import matplotlib.pyplot as plt

import ketju_integrator
import ketju_utils

def keplerian_orbit(m1, m2, a, e):
    # Keplerian orbit at apocentre, G==1 units
    M = m1+m2
    n = np.sqrt(M/a**3)
    phi = np.pi
    r = a*(1+e)

    phidot = np.sqrt((1-e)/(1+e)) * n/(1+e)
    x = np.array([-r, 0., 0.])
    v = np.array([0., r*phidot, 0.])
    x1 = -m2/M*x
    x2 = m1/M*x
    v1 = -m2/M*v
    v2 = m1/M*v

    return ketju_utils.Particle(0,m1,x1,v1), ketju_utils.Particle(0,m2,x2,v2)

# masses in solar masses
m1 = 15
m2 = 10
Rs = 2*(m1+m2)

# (approximate) initial orbital parameters
a0 = 25 # in units of combined Schwarzschild radius Rs
e0 = 0.1

# spins
chi1 = .3 
chi1_dir = np.array([1.,0,0])

chi2 = .7
chi2_dir = np.array([0.,0,1])

chi1_dir /= np.linalg.norm(chi1_dir) 
chi2_dir /= np.linalg.norm(chi2_dir)


bh1, bh2 = keplerian_orbit(m1, m2, a0*Rs, e0)

# set spins
bh1.spin[0] = m1**2 * chi1 * chi1_dir
bh2.spin[0] = m2**2 * chi2 * chi2_dir

init_state = ketju_utils.particles_to_integrator_state([bh1,bh2])

# Set integration time based on initial approximate orbital period T
T = 2*np.pi*(Rs*a0)**1.5/(m1+m2)**.5
prop_t = 1e3 * T


# Setting the PN corrections used can be done like this
PNF = ketju_integrator.PNFlags
#pnflags = PNF.PN_1_0_ACC + PNF.PN_2_5_ACC + PNF.PN_3_5_ACC
pnflags = PNF.PN_ALL

# Integration output times
ts = np.linspace(0, prop_t, 10000)

# Run integrator in subprocess mode, see the direct mpi example for how to
# change this if it doesn't work on your system.
ret = ketju_integrator.run_integrator(
                     init_state, ts, nproc=2, gbs_relative_tolerance=1e-8,
                     PN_flags=pnflags, G=1, c=1,
# BH merger distance set to a lower value than the default to illustrate the
# problematic small separation behavior
                     bh_merger_distance_in_Rs=6,
                     )

if ret['merger_infos']:
    print("Had merger")
    print(ret['merger_infos'][0])

# Transform to ketju_utils.Particle objects and crop the particles to the same
# length in case of merger during the integration.
bh1, bh2 = ketju_utils.integrator_states_to_particles(ret['states'])
bh1 = bh1[:len(bh2)]

#plot orbit
fig, (ax1,ax2) = plt.subplots(1,2)
ax1.plot(bh1.x[:,0]/Rs, bh1.x[:,1]/Rs, alpha=.6)
ax1.plot(bh2.x[:,0]/Rs, bh2.x[:,1]/Rs, alpha=.6)
ax2.plot(bh1.x[:,0]/Rs, bh1.x[:,2]/Rs, alpha=.6)
ax2.plot(bh2.x[:,0]/Rs, bh2.x[:,2]/Rs, alpha=.6)
ax1.set_aspect('equal')
ax1.set_xlabel('$x/R_s$')
ax1.set_ylabel('$y/R_s$')

ax2.set_aspect('equal')
ax2.set_xlabel('$x/R_s$')
ax2.set_ylabel('$z/R_s$')
if ax2.get_ylim()[1] - ax2.get_ylim()[0] < 1:
    ax2.set_ylim(-.5*a0, .5*a0)


# PN accurate orbital parameters, not accounting for spin effects.
# Non-zero spins will cause oscillations in the calculated parameters.
# PN_level should be set to match the used integer level PN corrections:
# PN_1_0 -> 2, PN_2_0 -> 4, PN_3_0 -> 6
pars = ketju_utils.orbital_parameters(bh1, bh2, PN_level=6)

fig, (ax0,ax1,ax2) = plt.subplots(3,1, sharex='col')

# Convert to seconds from G=c=1 units
t_in_s = pars['t']/ketju_utils.units.s 
ax0.plot(t_in_s, pars['R']/Rs)
ax0.set_ylabel('$R/R_s$')

ax1.plot(t_in_s, pars['a_R']/Rs)
ax1.set_ylabel('$a/R_s$')

ax2.plot(t_in_s, pars['e_t'])
ax2.set_ylabel('$e_t$')
ax2.set_xlabel('$t/\mathrm{s}$')


# CoM velocity and binary separation, useful for checking that the PN level
# errors in the CoM velocity aren't too large for you merger separation setting.
fig, (ax0, ax1) = plt.subplots(2,1, sharex='col')

ax0.plot(t_in_s, pars['R']/Rs)
ax0.set_ylabel('$R/R_s$')

ax1.plot(t_in_s, ketju_utils.binary_CoM_vel(bh1,bh2,PN_level=6)/ketju_utils.units.km_per_s)
ax1.set_ylabel(r'$v_\mathrm{CoM}/\mathrm{km\,s^{-1}}$')
ax1.set_xlabel('$t/$s')


plt.show()
