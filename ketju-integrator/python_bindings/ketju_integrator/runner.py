# Ketju-integrator Python bindings
# Copyright (C) 2020-2023 Matias Mannerkoski and contributors.
# Licensed under the GPLv3 license.
# See LICENSE.md or <https://www.gnu.org/licenses/> for details.

"""
A helper module to run the integrator using subprocesses.
Launched from the run_integrator() and load_state_dump() functions,
and communicates with the main process using pickle and piped stdout.
"""

import sys
import pickle
import ctypes

from .bindings import shim, _integrate

def _integrate_as_subprocess(this_task):
    if this_task == 0:
# Read the input data
        indata = pickle.loads(sys.stdin.buffer.read())
        instate = indata['instate']
        ts = indata['ts']
        c = indata['c']
        G = indata['G']
        integrator_options = indata['integrator_options']
    else:
        instate = ts = integrator_options = c = G = None

    retval = _integrate(this_task, instate, ts, c, G, integrator_options)

# Dump the results to stdout
    if this_task == 0:
        sys.stdout.buffer.write(pickle.dumps(retval, protocol=pickle.HIGHEST_PROTOCOL))
        sys.stdout.flush()
    

def _load_dump(fname, this_task):
    system = CSystem()
    fnb = fname.encode('utf-8')
    shim.load_state_dump(system, fnb)
    state = State.from_c_system(system)
    options = struct_to_dict(system.options[0])
    constants = struct_to_dict(system.constants[0])
    retval = dict(state=state, options=options, constants=constants)
    shim.free_system(system)
    if this_task == 0:
        sys.stdout.buffer.write(pickle.dumps(retval, protocol=pickle.HIGHEST_PROTOCOL))


def _run_main(this_task):
    if sys.argv[1] == 'integrate':
        _integrate_as_subprocess(this_task)
    elif sys.argv[1] == 'load_dump':
        _load_dump(sys.argv[2], this_task)

if __name__ == '__main__':
    this_task = shim.init_mpi()
    _run_main(this_task)
    shim.finalize_mpi()
