# Ketju-integrator Python bindings
# Copyright (C) 2020-2023 Matias Mannerkoski and contributors.
# Licensed under the GPLv3 license.
# See LICENSE.md or <https://www.gnu.org/licenses/> for details.

"""
A simple example of using the bindings directly with mpirun for a system with two BHs.
Run using e.g. "mpirun -np 2 python3 ketju_integrator_example_mpi.py". 
Run either with the module installed or from this directory or with this directory on the python 
search path for the below import to work.
"""
from ketju_integrator import *

import numpy as np
import matplotlib.pyplot as plt

# To avoid depending on additional mpi libraries, for this example the mpi
# setup is handled using some functions from the bindings shim.
# For more serious use, consider using something like mpi4py.
# With mpi4py you only need to do
# 
# from mpi4py import MPI
# mpi_rank = MPI.COMM_WORLD.Get_rank()
# 
# and can remove the init and finalize calls that use the shim

from ketju_integrator.bindings import shim
mpi_rank = shim.init_mpi()

# Initial state of the system
# initial state is only needed on root task

if mpi_rank == 0:
    init_state = State(0., np.array([1., 1.]),
                       np.array([[40.,0,0], [-40.,0,0]]),
                       np.array([[0.,.085, 0], [0.,-0.085,0]]),
                       np.array([[1.,0,0], [0.,1.,0.]])
                       )
# Output times
    ts = np.linspace(0,1e5, 1000)
else:
    init_state = ts = None

# run the integrator
ret = run_integrator(init_state, ts, mpi_rank=mpi_rank, PN_flags=PNFlags.PN_ALL)

if mpi_rank == 0:
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


shim.finalize_mpi()

