The CHIMES flag is the 'top-level' switch. When this is switched on, the following changes occur: 
  - A 'global variables' structure for the CHIMES parameters is created. 
  - An array of 'gas variables' structures is created, one for each SPH particle. The array is 
    'MaxPartSph' long on each processor (same as SphP). However, memory for the abundance 
    arrays is only allocated for the actual particles. For example, if an Sph particle is 
    moved to another processor, the memory for the abundance array of that particle is freed on 
    the local processor. This reduces the amount of memory needed (rather than always having 
    memory allocated for abundance arrays of 'MaxPartSph' particles on each processor). 
  - Parameters specific to CHIMES are read in from the parameter file. 
  - In the 'DoCooling()' routine, the internal energy is updated by calling the CHIMES chemistry 
    solver for that particle, which evolves the abundances and temperature together. 
  - The CHIMES abundance arrays are read in from ICs, written to snapshots, and are included in 
    restart files. Currently, the entire array is written out, which can be quite large. In the future, 
    we might want to look at how to reduce this - we could define a new 'class' of small snapshot that 
    includes only a handful of species and is written out more frequently (e.g. for creating movies), 
    or we could specify that the full abundance array is only written out every Nth snapshot. 
  - To control which elements are included in CHIMES, you can use the 'IncludeXXX' parameters in the 
    parameter file. Excluding elements will make the network run faster and will reduce the memory 
    needed for the abundance arrays, but the cooling from that element will then be excluded 
    altogether, so this should only be done for elements that you are sure will not contribute 
    to the cooling. Note that the included elements will also depend on other Config flags: 
    - If METALS is NOT defined, He is set according to XH (from allvars.h) and all element abundances of 
      Carbon and above are set to zero. Note that if the 'IncludeXXX' parameters are set to 1, they will still 
      allocate memory to the metal ions, they will just be set to zero. 
    - To include metals, you need to set both METALS and COOL_METAL_LINES_BY_SPECIES. Then each element 
      abundance is set from the mass fractions given in the P[i].Metallicity array, and we calculate XH 
      from X_He and Z, i.e. XH = 1.0 - X_He - Z. 
      NOTE: If you set METALS and CHIMES but not COOL_METAL_LINES_BY_SPECIES, the code will exit. To 
      include metal cooling in CHIMES, it needs to know the abundances of each metal element separately, 
      but these are not tracked separately in GIZMO without COOL_METAL_LINES_BY_SPECIES, hence both 
      are required. 
    WARNING: Currently, the ion/molecule abundances are set in the ICs file. However, there is no check 
    to ensure that the particle metallicities (and hence element_abundances, which are computed from 
    the particle metallicities) are consistent with the initial ion/molecule abundances. If they are 
    inconsistent, the chemistry solver will force the ion/molecule abundances to be consistent with 
    the metallicities, BUT the result may not be what you expected (for example, it may no longer be 
    in an initial equilibrium). 
  - In the domain decomposition, extra weight is given to high-density SPH particles, to try to better 
    balance the work load of the chemistry. NB: Currently this is commented out - we need to find what 
    works best. 
  - When CHIMES is switched on, cpu.txt contains an additional line for 'sfrcoolimbal'. This is CPU time 
    spent by tasks waiting for the other tasks to finish looping through their active particles and 
    doing the cooling (and star formation). 

There are then a few options for how to deal with self-shielding from the UV field. If none of the below 
options are set, then the gas is treated as optically thin to the applied UV field. To include self-shielding, 
you need to include one (and only one) of the following options, and also set the CellSelfShielding_On 
parameter to '1' or '2' in the parameter file (see CHIMES manual for details). 
  - CHIMES_SOBOLEV_SHIELDING - uses hydrogen column density N_H = n_H * L_sob, where: 
    L_sob = Shielding_length_factor * (rho / abs(grad(rho))), 
    and the Shielding_length_factor is set in the parameter file. 
Also, the H2 self-shielding will depend on Doppler broadening of the H2 lines. The thermal Doppler broadening 
is calculated from the gas temperature, and is automatically included within the chemistry solver, but the 
non-thermal broadening (due to turbulent motions) does need to be specified. Currently, this is just set 
to b = 7.1 km/s in the code (in cooling.c), but in the future we could add an option to compute this 
explicitly from the gas motions in the simulation. 

Additional CHIMES Config flags: 
  - CHIMES_REDUCED_OUTPUT - if this option is enabled, the full CHIMES abundance array will only be output
    every N_chimes_full_output_freq snapshots. For the remaining snapshots, only the abundances of 
    electrons, HI and H2 will be output. 
  - CHIMES_HYDROGEN_ONLY - sets the HYDROGEN_MASS_FRACTION to 1. This option is ignored if METALS are also set. 
  - CHIMES_NH_OUTPUT - writes out the column densities used for local self-shielding in CHIMES. 

Changes unrelated to CHIMES:
  - SOLAR_ABUNDANCES_WIERSMA09 - this gives an option to set the solar abundances to the values from Table 1 
    of Wiersma et al. 2009, MNRAS, 393, 33 (as used by e.g. OWLS, EAGLE). If this flag is not used, it will 
    revert to GIZMO's default solar abundances. 
  - Previously, the 'InitMetallicity' parameter was only included when certain FB flags were switched on.
    I have now changed this (in begrun.c, init.c and allvars.h) so that 'InitMetallicity' is included whenever 
    METALS are switched on.
  - I have also made changes to gravity/analytic_gravity.h, so use a Hernquist static potential, but these 
    changes SHOULD NOT be added to the repo! I have saved the 'original' version as 
    gravity/analytic_gravity_ORIGINAL.h. Make sure you revert back to this version!! 
  - GALSF_FB_NOENRICHMENT - this disables the injection of metals from SNe and winds, so that metallicity can 
    be held constant. Mass, momentum and energy are still injected. This change has been implemented in 
    galaxy_sf/mechanical_fb.c, lines 920-923. 


There are some Config flags that CHIMES does not currently work with, or its behavior may be unexpected: 
  - GALSF_EFFECTIVE_EQS - not currently compatible with CHIMES. 
  - GALSF_TURNOFF_COOLING_WINDS - currently, if cooling has been disabled for a wind particle, it also 
    will not evolve the chemical abundances either, i.e. the abundances become 'frozen' for the 
    delay time. 
  - GALSF_FB_FIRE_RT_HIIHEATING - if this option is used with CHIMES, then, for gas particles that are flagged as 
    being 'ionised' by a local star, the minimum temperature is set to HIIRegion_Temp. In other words, 
    particles below this are instantly heated to this temperature, and are not allowed to cool below it 
    for a time of DelayTimeHII. However, the ions are left to evolve in non-equilibrium, i.e. we do not 
    instantly ionise the particle (just heat it). Also, if CHIMES_LOCAL_UV is enabled, we explicitly include 
    the UV radiation from local stars, both for photoionisation and for photoheating, so this should explicitly 
    model HII regions, and then GALSF_FB_FIRE_RT_HIIHEATING is not required. 
  - GALSF_FB_FIRE_RT_UVHEATING - currently, switching this on would have no effect on CHIMES. We need to add a 
    UV field in CHIMES and use the SphP[i].Rad_Flux_UV of each gas particle. Also, note that, in accel.c, the 
    Rad_Flux_UV of each SPH particle is attenuated by local self-shielding by calling selfshield_local_incident_uv_flux(). 
    We do NOT want to do this when using it with CHIMES, because the chemistry solver will apply local self-shielding 
    itself. Therefore, you need to add an exception for CHIMES in accel.c. ALSO, what is the difference between 
    SphP[i].Rad_Flux_UV and SphP[i].Rad_Flux_EUV?
  - BLACK_HOLES - currently, CHIMES is not compatible with the inclusion of black holes, because when a gas particle 
    is accreted, we will need to delete the corresponding gasVariables structure (this is not currently done). 
  - RT_CHEM_PHOTOION - this is currently incompatible with CHIMES, as it uses its own routines to compute chemistry 
    abundances for the radiative transport. 
  - OUTPUTCOOLRATE - not currently compatible with CHIMES. 
  - COSMIC_RAY_FLUID - not currently compatible with CHIMES. 
  - Hydro methods: CHIMES is currently only compatible with Lagrangian methods (i.e. SPH or MFM), because the 
    chemical species are 'fixed' to the particles. If you run with e.g. MFV or moving mesh, the species would 
    need to be advected between cells, which is not currently implemented. Also, CHIMES is incompatible with 
    TURB_DIFF_METALS, because the chemical  species would need to follow the diffusion of metals, which is 
    not currently implemented. 
