Ketju-integrator demo N-body program
====================================

This is a simple demonstration program that can use the integrator to simulate
N-body systems.
After building by following the instructions in the main README, the `nbody`
executable is found in the `bin` directory.

This program is intended mainly as a simple example which can serve as a basis
for building more complete programs.
The features are therefore quite limited, and for science applications it is 
likely best to implement your own program to drive the integrator.

Initial conditions
------------------

The initial conditions file is a text file containing a header row with two integers  
`num_bhs num_stars`  
followed by `num_bhs` rows of BH particle data  
`mass x y z vx vy vz Sx Sy Sz`  
where x, y, z are the coordinates, vx, vy, vz the velocity components and
Sx, Sy, Sz the spin angular momentum components.
These are then followed by `num_stars` rows of star particle data
`mass x y z vx vy vz`  

The units are by default based on GADGET default units, i.e. 1e10 Msol for mass,
kpc for coordinates and km/s for velocities.  
To change the unit system, you can edit the code to set the values of G and c
in `R.constants`


There is very little error checking in the input code to keep it simple, 
so mistakes in the IC file can lead to unclear failures.

A sample IC file is included as `nbody_sample_ic.txt`.

Running
-------
Run by optionally calling through the MPI launcher on your system (mpirun etc.),
and giving as arguments the IC file, the integration time, number of output times
spread over the integration time, and the integration tolerance to use.

For example, assuming current working directory is the repository main directory:   
`mpirun -np 2 bin/nbody nbody/nbody_sample_ic.txt 1e-3 10 1e-8`  
to run the code using two MPI tasks.


Output
------

The code writes output in a `nbody_output.txt` file in the current working directory.
The file format consists of rows containing the different state variables of the system,
see the `write_output_file` function in `main.c` for details.

The output can be plotted using the included `plot_output.py` script, by running  
`python3 nbody/plot_output.py nbody_output.txt`  
This provides a simple plot of the trajectories of the particles.
