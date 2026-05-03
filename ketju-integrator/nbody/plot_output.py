# Ketju-integrator sample N-body driver program.
# Copyright (C) 2020-2022 Matias Mannerkoski and contributors.
# Licensed under the GPLv3 license.
# See LICENSE.md or <https://www.gnu.org/licenses/> for details.

import sys

import numpy as np
import matplotlib.pyplot as plt

# Reads the file specified as input, produces a simple plot
datafile = sys.argv[1]

ts = []
mass = []
pos = []
vel = []
spin = []

def grouped_lines(f):
    return zip(*([iter(f)]*5)) 

with open(datafile, 'r') as f:
    for head, ml, pl, vl, sl in grouped_lines(f):
        head_ss = head.split()
        ts.append(float(head_ss[0]))
        N_part = int(head_ss[1])
        N_PN = int(head_ss[2])
        mass.append(np.fromstring(ml, sep=' ', dtype=float,count=N_part))
        pos.append(np.fromstring(pl, sep=' ', dtype=float,count=3*N_part).reshape(N_part,3))
        vel.append(np.fromstring(vl, sep=' ', dtype=float,count=3*N_part).reshape(N_part,3))
        if N_PN > 0:
            spin.append(np.fromstring(sl, sep=' ', dtype=float,count=3*N_PN).reshape(N_PN,3))
        else:
            spin.append([])

# array lenghts change when BHs merge, need to handle each section with different counts separately
merger_indices = np.append(np.nonzero(np.diff(list(map(len,spin))))[0], [-1])

prev_ind = 0
for ind in merger_indices:
    pos_arr = np.array(pos[prev_ind:ind+1 if ind != -1 else None])
    N_BH = len(spin[ind])
# plot only positions in xy-plane for now
    for i in range(N_BH):
        plt.plot(*pos_arr[:,i,:2].T, color='k')
    for i in range(N_BH, len(pos_arr[0])):
        plt.plot(*pos_arr[:,i,:2].T, color='gray')

plt.plot(*pos_arr[-1,:N_BH,:2].T, '.', color='k')
plt.plot(*pos_arr[-1,N_BH:,:2].T, '*', color='gray')

plt.show()

