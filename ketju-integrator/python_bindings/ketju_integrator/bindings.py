# Ketju-integrator Python bindings
# Copyright (C) 2020-2023 Matias Mannerkoski and contributors.
# Licensed under the GPLv3 license.
# See LICENSE.md or <https://www.gnu.org/licenses/> for details.

"""
Implementation of the bindings.
"""
import subprocess
import io
import os
import sys
import enum
import ctypes
import pickle

import numpy as np

__all__ = ['State', 'PNFlags', 'run_integrator', 'load_state_dump']

# Struct definitions
class CSystemPhysicalState(ctypes.Structure):
    _fields_ = [
        ('time', ctypes.c_double),
        ('mass', ctypes.POINTER(ctypes.c_double)),
        ('pos' , ctypes.POINTER(ctypes.c_double)),
        ('vel' , ctypes.POINTER(ctypes.c_double)),
        ('spin', ctypes.POINTER(ctypes.c_double)),
    ]

class CIntegratorOptions(ctypes.Structure):
    _fields_ = [
        ('gbs_kmax', ctypes.c_int),
        ('gbs_relative_tolerance', ctypes.c_double),
        ('gbs_absolute_tolerance', ctypes.c_double),
        ('output_time_relative_tolerance', ctypes.c_double),
        ('output_time_absolute_tolerance', ctypes.c_double),
        ('PN_flags', ctypes.c_uint),
        ('steps_between_mst_reconstruction', ctypes.c_int),
        ('max_tree_distance', ctypes.c_int),
        ('max_step_count', ctypes.c_int),
        ('mst_algorithm', ctypes.c_int),
        ('star_star_softening', ctypes.c_double),
        ('no_star_star_gravity', ctypes.c_bool),
        ('enable_bh_mergers', ctypes.c_bool),
        ('bh_merger_distance_in_Rs', ctypes.c_double),
        ('enable_bh_merger_kicks', ctypes.c_bool),
    ]

class CNaturalConstants(ctypes.Structure):
    _fields_ = [
        ('G', ctypes.c_double),
        ('c', ctypes.c_double),
    ]

class CPerformanceCounters(ctypes.Structure):
    _fields_ = [
        ("time_force", ctypes.c_double),
        ("time_force_comm", ctypes.c_double),
        ("time_coord", ctypes.c_double),
        ("time_mst", ctypes.c_double),
        ("time_gbscomm", ctypes.c_double),
        ("time_gbsextr", ctypes.c_double),
        ("time_total", ctypes.c_double),
        
        ("time_tprim1", ctypes.c_double),
        ("time_tprim2", ctypes.c_double),
        ("time_tprim3", ctypes.c_double),
        ("time_tprim4", ctypes.c_double),
        ("time_tprim5", ctypes.c_double),

        ("successful_steps", ctypes.c_int),
        ("failed_steps", ctypes.c_int),
        ("total_work", ctypes.c_double),
        ("relative_energy_error", ctypes.c_double),
    ]

class CBHMergerInfo(ctypes.Structure):
    _fields_ = [
        ("t_merger", ctypes.c_double),
        ("m1", ctypes.c_double),
        ("m2", ctypes.c_double),
        ("m_remnant", ctypes.c_double),
        ("chi1", ctypes.c_double),
        ("chi2", ctypes.c_double),
        ("chi_remnant", ctypes.c_double),
        ("vkick", ctypes.c_double),
        ("extra_index1", ctypes.c_int),
        ("extra_index2", ctypes.c_int),
    ]

class CSystem(ctypes.Structure):
    _fields_ = [
        ('num_particles', ctypes.c_int),
        ('num_pn_particles', ctypes.c_int),
        ('physical_state', ctypes.POINTER(CSystemPhysicalState)),
        ('internal_state', ctypes.c_void_p),
        ('perf', ctypes.POINTER(CPerformanceCounters)),
        ('options', ctypes.POINTER(CIntegratorOptions)),
        ('constants', ctypes.POINTER(CNaturalConstants)),

        ('particle_extra_data', ctypes.c_void_p),
        ('particle_extra_data_elem_size', ctypes.c_size_t),

        ('num_mergers', ctypes.c_int),
        ('merger_infos', ctypes.POINTER(CBHMergerInfo)),
    ]

def struct_to_dict(struct):
    return {field : getattr(struct, field) for field, _ in struct._fields_}

class PNFlags(enum.IntFlag):
    PN_NONE = 0
    PN_1_0_ACC = 1 << 0
    PN_THREEBODY = 1 << 1
    PN_2_0_ACC = 1 << 2
    PN_2_5_ACC = 1 << 3
    PN_3_0_ACC = 1 << 4
    PN_3_5_ACC = 1 << 5
    PN_CONSERVATIVE = PN_1_0_ACC | PN_2_0_ACC | PN_3_0_ACC
    PN_DYNAMIC_ALL = 0xFF 

    PN_1_5_SPIN_ACC = 1 << 8
    PN_2_0_SPIN_ACC = 1 << 9
    PN_SPIN_EVOL = 1 << 10
    PN_SPIN_ALL = 0x700

    PN_ALL = PN_DYNAMIC_ALL | PN_SPIN_ALL

class State:
    """
    State of the system at a given time t.
    Used for input and output of results.
    
    Attributes:
    time: current physical time
    mass: numpy array (shape (num_particles,)) containing masses of the system
    pos: numpy array (shape (num_particles, 3)) containing particle coordinates
    vel: numpy array (shape (num_particles, 3)) containing particle velocities
    spin: numpy array (shape (num_pn_particles, 3)) containing particle spins
    ids: numpy array (shape (num_particles,), dtype=int) containing particle ids.
         These can be used to track particles through potential mergers during
         the integration.
    """
    def __init__(self, time, mass, pos, vel, spin=None, ids=None):
        """
        Create a state from given arraylike values (see above).
        Spin field is optional and its length determines the number of PN particles if given.
        ids field is optional and will be autogenerated if not given.
        """
        self.time = time
        
        self.mass = np.asarray(mass)
        if len(self.mass.shape) != 1:
            raise ValueError("mass must be 1d array-like")

        self.pos = np.atleast_2d(pos)
        if len(self.pos.shape) != 2 or self.pos.shape[0] != len(mass) or self.pos.shape[1] != 3:
            raise ValueError("pos must have shape (len(mass),3)")

        self.vel = np.atleast_2d(vel)
        if self.vel.shape != self.pos.shape:
            raise ValueError("vel must have shape (len(mass), 3)")


        if spin is not None:
            self.spin = np.atleast_2d(spin)
            if self.spin.shape[0] > self.pos.shape[0] or self.spin.shape[1] != 3:
                raise ValueError("spin must have shape (N, 3), N must be at most the total number of particles")
        else:
            self.spin = None
        if ids is None:
            self.ids = np.arange(0, len(self.mass), dtype=int)
        else:
            self.ids = np.asarray(ids, dtype=int)
            if self.ids.shape != self.mass.shape:
                raise ValueError("ids must have shape (len(mass),)")

    @classmethod
    def from_c_system(cls, system):
        st = system.physical_state[0]
        N = int(system.num_particles)
        time = float(st.time)
        mass = np.copy(np.ctypeslib.as_array(st.mass, shape=(N,)))
        pos = np.copy(np.ctypeslib.as_array(st.pos, shape=(N,3)))
        vel = np.copy(np.ctypeslib.as_array(st.vel, shape=(N,3)))
        if system.num_pn_particles > 0:
            spin = np.copy(np.ctypeslib.as_array(st.spin, shape=(system.num_pn_particles,3)))
        else:
            spin = None
        ids_ctype = np.ctypeslib.as_ctypes_type(np.dtype(int))
        if system.particle_extra_data is not None:
            ids = np.copy(np.ctypeslib.as_array(
                            ctypes.cast(system.particle_extra_data,ctypes.POINTER(ids_ctype)),
                            shape=(N,),))
        else:
            ids = None
        return cls(time, mass, pos, vel, spin, ids)

    def to_c_system(self, system):
        st = system.physical_state[0]
        N = int(system.num_particles)
        st.time = self.time
        np.ctypeslib.as_array(st.mass, shape=(N,))[:] = self.mass[:]
        np.ctypeslib.as_array(st.pos, shape=(N,3))[:] = self.pos[:]
        np.ctypeslib.as_array(st.vel, shape=(N,3))[:] = self.vel[:]
        if self.spin is not None:
            np.ctypeslib.as_array(st.spin, shape=(system.num_pn_particles,3))[:] = self.spin[:]

# XXX: the data for the ids is stored in this object, returned from this function
# The caller must keep it alive for the duration of the integration.
        id_storage = np.copy(self.ids)
        system.particle_extra_data = id_storage.ctypes.data
        ids_ctype = np.ctypeslib.as_ctypes_type(np.dtype(int))
        system.particle_extra_data_elem_size = ctypes.sizeof(ids_ctype)
        return id_storage



def run_integrator(instate, ts, c=1., G=1., nproc=1, mpi_rank=None, **integrator_options):
    """
    Run the integrator starting from the State in instate storing states at the
    times given in ts.

    Can be run in two modes:
        subprocess mode (default): launches nproc processes using mpirun to perform the integration,
            communicating with them through their stdout. Not very robust in cluster environments,
            but allows easy running on typical desktop machines e.g. from an interactive session.
        MPI mode: specified by setting mpi_rank to an integer value specifying the current rank.
            In this case this function is called from python processes launched with e.g. mpirun,
            and the MPI context has been initialized and is managed using e.g. mpi4py.
            All tasks in COMM_WORLD must call this function, but on tasks with mpi_rank != 0
            all the other arguments can be set to None and are ignored.

    Arguments
    ---------
    instate: A state object describing the system at the initial time.
    ts: an array of times t for output states from the integration, does not need to contain the initial time.
    c,G: speed of light and gravitational constant to use for the integration, defines the unit system used.
    nproc: how many processes to use for the integration, only relevant for subprocess mode.
    mpi_rank: the rank of the current task in MPI_COMM_WORLD when run in MPI mode.
    **integrator_options: any options that can be set in the ketju_integrator_options struct in C.

    Returns
    -------
    A dict with keys:
        'states': a list of State objects at the requested times
        'perf': a dict of the performance counters from the integrator
        'merger_infos': a list of dicts of describing the BH mergers during the integration.

    Only returned on the root task in when run in MPI mode.
    """
    if mpi_rank is None:
        data = dict(ts=ts, instate=instate, integrator_options=integrator_options, c=c, G=G)
# Run the runner module using multiple tasks with mpirun.
#TODO mpirun might not be the correct one to call on all systems, maybe make specifiable.
# Pass the input data through stdin, get the output data through stdout as pickled objects
#XXX: killing/interrupting the parent task can sometimes leave the workers behind,
# even though mpirun should handle killing them in that case. How to fix that?
        ret = subprocess.check_output(['mpirun', '-np', str(nproc), sys.executable,
                                       '-m', 'ketju_integrator.runner',
                                       'integrate'],
                        input=pickle.dumps(data, protocol=pickle.HIGHEST_PROTOCOL))
        return pickle.loads(ret)

    else:
        return _integrate(mpi_rank, instate, ts, c, G, integrator_options)

def load_state_dump(fname):
    return _unpickle_subprocess_output(
            subprocess.check_output([sys.executable, '-m', 'ketju_integrator.runner',
                                     'load_dump', fname]))


def _load_shim():
# Load the shim library
    _this_dir = os.path.dirname(os.path.realpath(__file__))
    shim_file = os.path.join(_this_dir, "libketju_integrator_python_shim.so")
    shim = ctypes.cdll.LoadLibrary(shim_file)
    shim.init_system.argtypes = [ctypes.POINTER(CSystem), ctypes.c_int, ctypes.c_int]
    shim.comm_initial_data.argtypes = [ctypes.POINTER(CSystem)]
    shim.comm_initial_data.restype = ctypes.c_int
    shim.init_mpi.restype = ctypes.c_int
    shim.free_system.argtypes = [ctypes.POINTER(CSystem)]
    shim.run_integrator.argtypes = [ctypes.POINTER(CSystem), ctypes.c_double]
    shim.load_state_dump.argtypes = [ctypes.POINTER(CSystem), ctypes.c_char_p]
    return shim
shim = _load_shim()

def _integrate(this_task, instate, ts, c, G, integrator_options):
# Since I didn't want to require an additional dependency like mpi4py for the
# minor amount of comms needed for setup, they are instead handled by the shim
# library and some values here are actually meaningful only on the root task
    num_ts = 0
    system = CSystem()
    n_pn_part = n_other_part = -1
    if this_task == 0:
        num_ts = len(ts)
        n_pn_part = len(instate.spin) if instate.spin is not None else 0
        n_other_part = len(instate.pos) - n_pn_part
# init the system (needs to be done at the same time on all tasks),
# the particle counts are communicated in the shim
    shim.init_system(system, n_pn_part, n_other_part)

    if this_task == 0:
# set the initial values
        id_storage = instate.to_c_system(system)
        system.constants[0].c = c
        system.constants[0].G = G
        for k,v in integrator_options.items():
            setattr(system.options[0], k, v)
# The ids are only handled on the root task.
# This is fine since all the extra data management stuff is handled separately
# on each task in the integrator

# comm the initial data and how many steps we need to take
    num_ts = shim.comm_initial_data(system, num_ts)

# run the integration
    if this_task == 0:
        outstates = []
        for t in ts:
            shim.run_integrator(system, t - system.physical_state[0].time)
            outstates.append(State.from_c_system(system))
    else:
        for _ in range(num_ts):
# The integration time is communicated in the shim.
            shim.run_integrator(system, 0.)

    if this_task == 0:
        perf = struct_to_dict(system.perf[0])
        merger_infos = [struct_to_dict(system.merger_infos[i]) for i in range(system.num_mergers)]
        retval = dict(states=outstates, perf=perf, merger_infos=merger_infos)
    else:
        retval = None

    shim.free_system(system)
    return retval
