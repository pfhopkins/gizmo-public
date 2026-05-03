# Ketju-integrator Python bindings
# Copyright (C) 2020-2022 Matias Mannerkoski and contributors.
# Licensed under the GPLv3 license.
# See LICENSE.md or <https://www.gnu.org/licenses/> for details.

"""
A simple example of using the bindings with two BHs
using subprocesses to run the calculation.
Run like any other python script.
Run either with the module installed or from this directory or with this directory on the python 
search path for the below import to work.
This subprocess mode may not work in all environments,
see the mpi mode example in that case.
"""
from ketju_integrator import *

import numpy as np
import matplotlib.pyplot as plt

# Initial state of the system
init_state = State(0., np.array([1., 1.]),
                   np.array([[40.,0,0], [-40.,0,0]]),
                   np.array([[0.,.085, 0], [0.,-0.085,0]]),
                   np.array([[1.,0,0], [0.,1.,0.]])
                   )
# Output times
ts = np.linspace(0,1e5, 1000)

# run the integrator
ret = run_integrator(init_state, ts, nproc=2, PN_flags=PNFlags.PN_ALL)

print("Perf counters:")
for k, v in ret['perf'].items():
    print(k, v)

ret_states = ret['states']


# Plot orbits
pos1 = np.array([st.pos[0] for st in ret_states])
pos2 = np.array([st.pos[1] for st in ret_states])
plt.plot(*pos1.T[:2])
plt.plot(*pos2.T[:2])
plt.show()

