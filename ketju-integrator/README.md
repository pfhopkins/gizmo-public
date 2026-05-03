# Ketju-integrator (MSTAR)

A highly parallel regularized N-body integrator library.

This code forms the core of the KETJU project, which adds detailed small-scale
modelling of supermassive black hole dynamics to galactic simulation codes such as
GADGET-4.
See the [KETJU project website] for more details on KETJU.
In addition to being used in our GADGET-based KETJU codes,
this integrator code is designed to be easily used in other codes as well.
To help with this, the repository includes both a simple example `nbody` program written in C,
as well as bindings that allow the integrator to be called from Python.

The MSTAR integrator is based on algorithmic regularization together with a
minimum spanning tree -based relative coordinate system to allow accurate integration
of even collisional orbits while reducing the impact of rounding errors.
The code also includes post-Newtonian corrections for modelling black hole dynamics and their mergers.
Good scaling behaviour to systems of several thousands of particles on hundreds of CPU cores
is achieved through MPI parallelization of both the calculation of gravitational interactions
and the Gragg–Bulirsch–Stoer extrapolation algorithm.
For more detail on the technical aspects of the integrator, see the [code paper] as well
as comments within the code.

When publishing work utilizing this code, please cite the [code paper].

If you encounter issues when using the code, you can report them using the [issue tracker].

[code paper]: https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.4131R
[KETJU project website]: https://www.mv.helsinki.fi/home/phjohans/ketju
[GADGET4-KETJU repository]: https://bitbucket.org/helsinkiextragalacticgroup/gadget4-ketju
[issue tracker]: https://bitbucket.org/helsinkiextragalacticgroup/ketju-integrator/issues 

## License

This software is licensed under the GPLv3 license:
you are free to modify and distribute it under certain conditions,
and it comes without any warranty.
See [LICENSE.md](LICENSE.md) for details.

The code included in the `thirdparty` directory is distributed under the licenses specified there.


## Building
The code requires an MPI library/environment (e.g. OpenMPI) to be installed before building.

To build the integrator library, first configure the CMake build system by running  

```
mkdir build
cd build
cmake .. [CMAKE_OPTIONS]
```

Configuration options such as specific compiler optimization flags can then be
changed by editing the `CMakeCache.txt` file either directly or using e.g. `ccmake`,
but doing this is not typically required.

After configuration, the library and other included components can be
built by running `make` in the build directory.
The `nbody` executable is then found under `bin/`.
See the [`nbody` README](nbody/README.md) for details on using it.

## Usage in other codes

The integrator is easiest to use from programs written in C or C++,
or through the included Python wrapper.
For use from other languages,
you may need to write a small wrapper for the library using C first.
You also need to include the MPI library in your code and perform the necessary
initialization for it before calling the integrator.

Examples of using the integrator in C and C++ can be found in the included `nbody`
code as well in the [GADGET4-KETJU repository].

### Use in C or C++

The public interface of the library is included in the file 
`include/ketju_integrator/ketju_integrator.h`,
so in your C or C++ code add  
`#include <ketju_integrator/ketju_integrator.h>`  
to use the definitions.
Details of the interface are also explained in the comments included in that header file.

The general outline of using the integrator looks something like this:

```C++
// Set up your MPI environment etc. with mpi_communicator being a communicator
// over the tasks you want to use for the integration (can be MPI_COMM_WORLD).
// All the tasks on this communicator must call the following functions and
// modify the integrator setting and initial state in the same way.

// Create the integrator data structures
// num_pn_particles is the number of particles (usually black holes)
// with PN terms enabled (depending on options).
// num_other_particles is the number of particles (usually stars)
// using only Newtonian gravity, possibly with softening.
struct ketju_system R;
ketju_create_system(&R, num_pn_particles, num_other_particles, mpi_communicator);

// Set the configuration in R.options and R.constants if you want to change the defaults.
// Initialize values in R.physical_state:
//   Set mass, pos, vel for R.num_particles (= num_pn_particles + num_other_particles) entries.
//   Set spins for the first num_pn_particles.
...

// Optional:
// If the initial data you just set is far from the center-of-mass frame of the system,
// you may benefit from shifting the system to its CoM frame using
double CoM_pos[3], CoM_vel[3];
ketju_into_CoM_frame(&R, CoM_pos, CoM_vel);
// to maintain the best numerical resolution for the system coordinates during the integration.
// CoM_pos and CoM_vel give the original coordinates of the CoM for later processing.

// Repeat as desired:
// Run the integration
ketju_run_integrator(&R, integration_time);
// Read the resulting state from R.physical_state and process it
...

// Finally free the allocated memory and communicators
ketju_free_system(&R);
```

You can have multiple independent integrator systems within your program at the same time,
as all of the integrator state is contained within the `ketju_system` struct.

#### Error handling

At the moment the integrator implements only simple error handling by calling an
internal `terminate()` function to print an error message and 
kill the program (by calling `MPI_Abort`)
in case the code ends up in an unrecoverable situation, likely due to an error in
the inputs or a bug in the code, or due to exceeding the number of steps given
by `max_step_count`.
A dump of the integrator state is also produced for some of these errors to help
with debugging.
These state dumps can be loaded for analysis using the `ketju_debug_read_integrator_state` function,
e.g. inside `gdb`, or through the `load_state_dump` function in the python bindings.

Otherwise, if the `ketju_run_integrator()` call returns, the integration has completed without errors.

### Building your code and linking the library

Finally, to build your code that uses the integrator, you need to make some changes to
your build system.
If you're building your code with `cmake`, you can use the following commands
in your `CMakeLists.txt` to add the necessary include and linker flags:  

```cmake
set(ketju-integrator_DIR "/path/to/ketju-integrator/") 
find_package(ketju-integrator 1.0 REQUIRED)  
add_executable(nbody "main.c") # using the nbody program as an example  
target_link_libraries(nbody ketju-integrator)
```

Otherwise, simply add (assuming gcc/clang)   
`-I/path/to/ketju-integrator/include/`  
to your compiler flags and    
`-L/path/to/ketju-integrator/lib -lketju-integrator`  
to your linker flags.
naturally you will also need to include the necessary MPI flags or use the 
MPI compiler wrapper to compile your code to link the necessary MPI library.



### Python bindings

Simple Python 3 bindings depending only on the standard library and NumPy are
included as the `ketju_integrator` module under `python_bindings`.

They function via a small C library shim and allow parallel integration by
running the integrator in subprocesses with MPI, or by running python directly with MPI.
Examples of using the bindings are given in
`python_bindings/examples/`.
After building the library as described above,
you can install the python bindings with `pip` by running 

```
python3 -m pip install --user ./python_bindings/
```

(you can omit `--user` if in a venv), you may need to first upgrade `pip` with
`python3 -m pip install --upgrade pip`.

The [`ketju_utils` package](https://bitbucket.org/helsinkiextragalacticgroup/ketju-utils)
provides some utilities for analysing BH orbits with post-Newtonian corrections,
and an example utilizing it is included.

## Test suite

If you are making changes to the integrator code itself, it is useful to also
build and run the included test suite.
This can be useful also for checking that the code functions correctly on your system.
To build the test suite, run `git submodule update --init`
to pull in the unit test libraries required for building the tests.
Then set `UNIT_TESTS=ON` in the CMake configuration and run the build as normal.

The tests can be run using `make run-unit-tests` or by running the test binary
found under `bin/` directly.
Note that some of the are currently somewhat sensitive to rounding errors caused
by different compilers and optimizations, so in case the tests fail inspect the
output to see if the errors in the results are only slightly over the assigned tolerances.


