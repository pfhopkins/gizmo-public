#ifndef ALLVARS_H
#include "allvars.h"
#endif
#include "gravity/forcetree.h"
#include "domain.h"
#ifdef COOLING
#include "cooling/cooling.h"
#endif
#ifdef BLACK_HOLES
#include "./galaxy_sf/blackholes/blackhole.h"
#endif

/* declarations of functions throughout the code */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * in part (adding/removing routines as necessary)
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifdef HAVE_HDF5
void write_header_attributes_in_hdf5(hid_t handle);
void read_header_attributes_in_hdf5(char *fname);
void write_parameters_attributes_in_hdf5(hid_t handle);
void write_units_attributes_in_hdf5(hid_t handle);
void write_constants_attributes_in_hdf5(hid_t handle);
#endif
void output_compile_time_options(void);

void report_pinning(void);
void pin_to_core_set(void);
void get_core_set(void);

void init_turb(void);
void set_turb_ampl(void);
void add_turb_accel(void);
void log_turb_temp(void);

void mpi_report_comittable_memory(long long BaseMem);
long long report_comittable_memory(long long *MemTotal,
                                   long long *Committed_AS,
                                   long long *SwapTotal,
                                   long long *SwapFree);

void merge_and_split_particles(void);
int does_particle_need_to_be_merged(int i);
int does_particle_need_to_be_split(int i);
double target_mass_renormalization_factor_for_mergesplit(int i);
void merge_particles_ij(int i, int j);
//void split_particle_i(int i, int n_particles_split, int i_nearest, double r2_nearest);
void split_particle_i(int i, int n_particles_split, int i_nearest);
double gamma_eos(int i);
void do_first_halfstep_kick(void);
void do_second_halfstep_kick(void);
double matrix_invert_ndims(double T[3][3], double Tinv[3][3]);
#ifdef HERMITE_INTEGRATION
int eligible_for_hermite(int i);
void do_hermite_prediction(void);
void do_hermite_correction(void);
#endif
#ifdef ADAPTIVE_TREEFORCE_UPDATE
int needs_new_treeforce(int i);
#endif
void find_timesteps(void);
#ifdef GALSF
void compute_stellar_feedback(void);
#endif
void compute_hydro_densities_and_forces(void);
void compute_grav_accelerations(void);
void calc_memory_checksum(void *base, size_t bytes);

void get_disk_forces(double RR, double zz, double *f_R, double *f_z);
double get_disk_mass(double time);
void growing_disk_init(void);

double get_turb_pot(double x, double y, double z);

void   sub_turb_move_perturbers(double t0, double t1);
void   sub_turb_add_forces(void);
void   sub_turb_read_table(void);
void   sub_turb_parent_halo_accel(double dx, double dy, double dz, double *acc);
double sub_turb_enclosed_mass(double r, double msub, double vmax, double radvmax, double c);

void interpolate_fluxes_opacities_gasgrains(void);
#if defined(RT_OPACITY_FROM_EXPLICIT_GRAINS)
double return_grain_extinction_efficiency_Q(int i, int k_freq);
#endif
int powerspec_turb_find_nearest_evaluate(int target, int mode, int *nexport, int *nsend_local);
void powerspec_turb_calc_dispersion(void);
double powerspec_turb_obtain_fields(void);
void powerspec_turb_save(char *fname, double *disp);
void powerspec_turb_collect(void);
void powerspec_turb(int filenr);
void compute_additional_forces_for_all_particles(void);


void set_cosmo_factors_for_current_time(void);
void drift_sph_extra_physics(int i, integertime tstart, integertime tend, double dt_entr);

void set_non_standard_physics_for_current_time(void);
void calculate_non_standard_physics(void);
void compute_statistics(void);
void execute_resubmit_command(void);
void make_list_of_active_particles(void);
void output_extra_log_messages(void);


static inline double WRAP_POSITION_UNIFORM_BOX(double x)
{
    while(x >= All.BoxSize) {x -= All.BoxSize;}
    while(x < 0) {x += All.BoxSize;}
    return x;
}

static inline double DMAX(double a, double b) { return (a > b) ? a : b; }
static inline double DMIN(double a, double b) { return (a < b) ? a : b; }
static inline int IMAX(int a, int b) { return (a > b) ? a : b; }
static inline int IMIN(int a, int b) { return (a < b) ? a : b; }
static inline double MINMOD(double a, double b) {return (a>0) ? ((b<0) ? 0 : DMIN(a,b)) : ((b>=0) ? 0 : DMAX(a,b));}
/* special version of MINMOD below: a is always the "preferred" choice, b the stability-required one. here we allow overshoot, just not opposite signage */
static inline double MINMOD_G(double a, double b) {return a;}
static inline double sigmoid_sqrt(double x) {return 0.5*(1 + x/sqrt(1+x*x));} /* Sigmoid ("turn-on") function (1 + x/(1+x^2))/2, interpolates between 0 as x->-infty and 1 as x->infty. Useful for cheaply doing smooth fits of e.g. EOS where different thermo processes turn on at certain temps */


#ifdef BOX_SHEARING
void calc_shearing_box_pos_offset(void);
#endif


int ngb_treefind_variable_threads_targeted(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
                                           int mode, int *exportflag, int *exportnodecount, int *exportindex,
                                           int *ngblist, int TARGET_BITMASK);

int ngb_treefind_pairs_threads_targeted(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
                                           int mode, int *exportflag, int *exportnodecount, int *exportindex,
                                           int *ngblist, int TARGET_BITMASK);



void do_distortion_tensor_kick(int i, double dt_gravkick);
void set_predicted_sph_quantities_for_extra_physics(int i);
void do_sph_kick_for_extra_physics(int i, integertime tstart, integertime tend, double dt_entr);
#if (SINGLE_STAR_TIMESTEPPING > 0)
void do_fewbody_kick(int i, double fewbody_kick_dv[3], double dt);
#endif

void check_particle_for_temperature_minimum(int i);

double get_pressure(int i);
double return_user_desired_target_density(int i);
double return_user_desired_target_pressure(int i);
#ifdef EOS_TILLOTSON
double calculate_eos_tillotson(int i);
void tillotson_eos_init(void);
#endif

void read_fof(int num);
int fof_compare_ID_list_ID(const void *a, const void *b);

void myfree_msg(void *p, char *msg);
void kspace_neutrinos_init(void);

#ifdef OUTPUT_TWOPOINT_ENABLED
void twopoint(void);
void twopoint_save(void);
int twopoint_ngb_treefind_variable(MyDouble searchcenter[3], MyFloat rsearch, int target, int *startnode, int mode, int *nexport, int *nsend_local);
int twopoint_count_local(int target, int mode, int *nexport, int *nsend_local);
#endif

void powerspec(int flag, int *typeflag);
double PowerSpec_Efstathiou(double k);
void powerspec_save(void);
void foldonitself(int *typelist);
void dump_potential(void);

int snIaheating_evaluate(int target, int mode, int *nexport, int *nSend_local);
void snIa_heating(void);
void voronoi_setup_exchange(void);

double get_neutrino_powerspec(double k, double ascale);
double get_powerspec(double k, double ascale);
void init_transfer_functions(void);


int MPI_Check_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                       int dest, int sendtag, void *recvbufreal, int recvcount,
                       MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status * status);

int MPI_Sizelimited_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
			     int dest, int sendtag, void *recvbuf, int recvcount,
			     MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status * status);

int mpi_calculate_offsets(int *send_count, int *send_offset, int *recv_count, int *recv_offset, int send_identical);
void sort_based_on_field(void *data, int field_offset, int n_items, int item_size, void **data2ptr);
void mpi_distribute_items_to_tasks(void *data, int task_offset, int *n_items, int *max_n, int item_size);

void parallel_sort_special_P_GrNr_ID(void);
void calculate_power_spectra(int num, long long *ntot_type_all);


int pmforce_is_particle_high_res(int type, MyDouble *pos);

void compare_partitions(void);
void assign_unique_ids(void);
int permut_data_compare(const void *a, const void *b);
void  generate_permutation_in_active_list(void);

void conduction(void);
void conduction_matrix_multiply(double *in, double *out);
double conduction_vector_multiply(double *a, double *b);
int conduction_evaluate(int target, int mode, double *in, double *out, double *sum,
			int *nexport, int *nsend_local);


void fof_get_group_center(double *cm, int gr);
void fof_get_group_velocity(double *cmvel, int gr);
int fof_find_dmparticles_evaluate(int target, int mode, int *nexport, int *nsend_local);
void fof_compute_group_properties(int gr, int start, int len);

#ifdef TURB_DIFF_DYNAMIC
double INLINE_FUNC Get_Particle_Size_for_turb(int i);
#endif
void parallel_sort(void *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *));
void parallel_sort_comm(void *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *), MPI_Comm comm);
int compare_IDs(const void *a, const void *b);
void test_id_uniqueness(void);
int compare_densities_for_sort(const void *a, const void *b);
int io_compare_P_ID(const void *a, const void *b);
int io_compare_P_GrNr_SubNr(const void *a, const void *b);
void drift_particle(int i, integertime time1);
void put_symbol(double t0, double t1, char c);
void write_cpu_log(void);
int get_timestep_bin(integertime ti_step);
const char* svn_version(void);
void find_particles_and_save_them(int num);
void lineofsight_output(void);
void sum_over_processors_and_normalize(void);
void absorb_along_lines_of_sight(void);
void output_lines_of_sight(int num);
integertime find_next_lineofsighttime(integertime time0);
integertime find_next_gridoutputtime(integertime ti_curr);
void add_along_lines_of_sight(void);
void do_the_kick(int i, integertime tstart, integertime tend, integertime tcurrent, int mode);
void x86_fix(void) ;

void *mymalloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, int linenr);
void *mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file, int line);

void *myrealloc_fullinfo(void *p, size_t n, const char *func, const char *file, int line);
void *myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line);

void myfree_fullinfo(void *p, const char *func, const char *file, int line);
void myfree_movable_fullinfo(void *p, const char *func, const char *file, int line);

void mymalloc_init(void);
void dump_memory_table(void);
void report_detailed_memory_usage_of_largest_task(size_t *OldHighMarkBytes, const char *label, const char *func, const char *file, int line);

double get_shear_viscosity(int i);

void kinetic_feedback_mhm(void);
int kin_compare_key(const void *a, const void *b);
void kinetic_evaluate(int target, int mode);

void bubble(void);
void multi_bubbles(void);
void bh_bubble(double bh_dmass, MyFloat center[3], MyIDType BH_id);
void find_CM_of_biggest_group(void);
int compare_length_values(const void *a, const void *b);
double rho_dot(double z, void *params);
double bhgrowth(double z1, double z2);

int fof_find_dmparticles_evaluate(int target, int mode, int *nexport, int *nsend_local);

double INLINE_FUNC Get_Particle_Size(int i);
double INLINE_FUNC Get_Gas_density_for_energy_i(int i);
double INLINE_FUNC Get_Particle_Expected_Area(double h);
double get_cell_Bfield_in_microGauss(int i);
double Get_Gas_Ionized_Fraction(int i);
double CR_calculate_adiabatic_gasCR_exchange_term(int i, double dt_entr, double gamma_minus_eCR_tmp, int mode);
double INLINE_FUNC Get_CosmicRayEnergyDensity_cgs(int i);
double CR_gas_heating(int target, double n_elec, double nH0, double nHcgs);
double Get_CosmicRayIonizationRate_cgs(int i);
#ifdef COSMIC_RAY_FLUID
void CalculateAndAssign_CosmicRay_DiffusionAndStreamingCoefficients(int i);
double INLINE_FUNC Get_Gas_CosmicRayPressure(int i, int k_CRegy);
double Get_CosmicRayGradientLength(int i, int k_CRegy);
double Get_CosmicRayStreamingVelocity(int i, int k_CRegy);
double CosmicRay_Update_DriftKick(int i, double dt_entr, int mode);
double CR_cooling_and_gas_heating(int target, double n_elec, double nH_cgs, double dtime_cgs, int mode);
double CR_energy_spectrum_injection_fraction(int k_CRegy, int source_type, double shock_vel, int return_index_in_bin, int target);
double return_cosmic_ray_anisotropic_closure_function_threechi(int target, int k_CRegy);
void inject_cosmic_rays(double CR_energy_to_inject, double injection_velocity, int source_type, int target, double *dir);
double return_CRbin_M1speed(int k_CRegy);
double evaluate_cr_transport_reductionfactor(int target, int k_CRegy, int mode);
double Get_AlfvenMachNumber_Local(int i, double vA_idealMHD_codeunits, int use_shear_corrected_vturb_flag);
double diffusion_coefficient_constant(int target, int k_CRegy);
double diffusion_coefficient_extrinsic_turbulence(int mode, int target, int k_CRegy, double M_A, double L_scale, double b_muG, double vA_noion, double rho_cgs, double temperature, double cs_thermal, double nh0, double nHe0, double f_ion);
double diffusion_coefficient_self_confinement(int mode, int target, int k_CRegy, double M_A, double L_scale, double b_muG, double vA_noion, double rho_cgs, double temperature, double cs_thermal, double nh0, double nHe0, double f_ion);
double return_CRbin_numberdensity_in_cgs(int target, int k_CRegy);
double return_CRbin_CR_energies_in_GeV(int target, int k_CRegy);
double return_CRbin_CR_charge_in_e(int target, int k_CRegy);
int return_CRbin_CR_species_ID(int k_CRegy);
double return_CRbin_kinetic_energy_in_GeV(int target, int k_CRegy);
double return_CRbin_gamma_factor(int target, int k_CRegy);
double gamma_eos_of_crs_in_bin(int k_CRegy);
double return_CRbin_beta_factor(int target, int k_CRegy);
double get_cell_Urad_in_eVcm3(int i);
void CR_cooling_and_losses(int target, double n_elec, double nHcgs, double dtime_cgs);
double return_CRbin_CRmass_in_mp(int target, int k_CRegy);
double return_CRbin_CR_rigidity_in_GV(int target, int k_CRegy);
double CR_get_streaming_loss_rate_coefficient(int target, int k_CRegy);
double Get_Gas_ion_Alfven_speed_i(int i);
double return_CRbin_nuplusminus_asymmetry(int i, int k_CRegy);
#endif
#ifdef EOS_ELASTIC
void elastic_body_update_driftkick(int i, double dt_entr, int mode);
#endif
#if defined(EOS_ELASTIC) || defined(EOS_TILLOTSON)
double get_negative_pressure_tensilecorrfac(double r, double h_i, double h_j);
#endif
double INLINE_FUNC convert_internalenergy_soundspeed2(int i, double u);
double INLINE_FUNC Get_Gas_effective_soundspeed_i(int i);
double INLINE_FUNC Get_Gas_thermal_soundspeed_i(int i);
double Get_Gas_Alfven_speed_i(int i);
double Get_Gas_Fast_MHD_wavespeed_i(int i);
double Get_Gas_Mean_Molecular_Weight_mu(double T_guess, double rho, double *xH0, double *ne_guess, double urad_from_uvb_in_G0, int target);
void update_explicit_molecular_fraction(int i, double dtime_cgs);
double return_dust_to_metals_ratio_vs_solar(int i);
double yhelium(int target);
double Get_Gas_Molecular_Mass_Fraction(int i, double temperature, double neutral_fraction, double free_electron_ratio, double urad_from_uvb_in_G0);
double INLINE_FUNC Get_Gas_BField(int i_particle_id, int k_vector_component);
#ifdef MAGNETIC
double Get_DtB_FaceArea_Limiter(int i);
#ifdef DIVBCLEANING_DEDNER
double INLINE_FUNC Get_Gas_PhiField(int i_particle_id);
double INLINE_FUNC Get_Gas_PhiField_DampingTimeInv(int i_particle_id);
#endif
#endif
#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
double INLINE_FUNC Get_Particle_Size_AGS(int i);
double get_particle_volume_ags(int j);
#endif

double INLINE_FUNC hubble_function(double a);
#ifdef GR_TABULATED_COSMOLOGY
double DarkEnergy_a(double);
double DarkEnergy_t(double);
#ifdef GR_TABULATED_COSMOLOGY_W
void fwa_init(void);
double INLINE_FUNC fwa(double);
double INLINE_FUNC get_wa(double);
#ifdef GR_TABULATED_COSMOLOGY_G
double INLINE_FUNC dHfak(double a);
double INLINE_FUNC dGfak(double a);
#endif
#endif
#endif

#ifdef GR_TABULATED_COSMOLOGY_H
double INLINE_FUNC hubble_function_external(double a);
#endif

void blackhole_accretion(void);
#ifdef BH_WIND_SPAWN
void get_random_orthonormal_basis(int seed, double *nx, double *ny, double *nz);
void get_wind_spawn_direction(int i, int num_spawned_this_call, int mode, double *ny, double *nz, double *veldir, double *dpdir);
#ifdef MAGNETIC
void get_wind_spawn_magnetic_field(int j, int mode, double *ny, double *nz,  double *dpdir, double d_r);
#endif
int blackhole_spawn_particle_wind_shell( int i, int dummy_sph_i_to_clone, int num_already_spawned );
void spawn_bh_wind_feedback(void);
#endif
int blackhole_evaluate(int target, int mode, int *nexport, int *nsend_local);
int blackhole_evaluate_swallow(int target, int mode, int *nexport, int *nsend_local);

int  blackhole_compare_key(const void *a, const void *b);


void fof_fof(int num);
void fof_import_ghosts(void);
void fof_course_binning(void);
void fof_find_groups(void);
void fof_check_cell(int p, int i, int j, int k);
void fof_find_minids(void);
int fof_link_accross(void);
void fof_exchange_id_lists(void);
int fof_grid_compare(const void *a, const void *b);
void fof_compile_catalogue(void);
void fof_save_groups(int num);
void fof_save_local_catalogue(int num);
void fof_find_nearest_dmparticle(void);
int fof_find_nearest_dmparticle_evaluate(int target, int mode, int *nexport, int *nsend_local);

int fof_compare_key(const void *a, const void *b);
void fof_link_special(void);
void fof_link_specialpair(int p, int s);
void fof_make_black_holes(void);

int io_compare_P_GrNr_ID(const void *a, const void *b);

void write_file(char *fname, int readTask, int lastTask);

void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *primary_taskID, int *last);

int get_values_per_blockelement(enum iofields blocknr);

int get_datatype_in_block(enum iofields blocknr);
void get_dataset_name(enum iofields blocknr, char *buf);


int blockpresent(enum iofields blocknr);
void fill_write_buffer(enum iofields blocknr, int *pindex, int pc, int type);
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type);

long get_particles_in_block(enum iofields blocknr, int *typelist);

int get_bytes_per_blockelement(enum iofields blocknr, int mode);

void read_file(char *fname, int readTask, int lastTask);

void get_Tab_IO_Label(enum iofields blocknr, char *label);


void long_range_init_regionsize(void);

int find_files(char *fname);

int metals_compare_key(const void *a, const void *b);
void enrichment_evaluate(int target, int mode);

void pm_init_nonperiodic_allocate(void);

void  pm_init_nonperiodic_free(void);

double get_random_number(MyIDType id);
void set_random_numbers(void);

int grav_tree_compare_key(const void *a, const void *b);
int dens_compare_key(const void *a, const void *b);
int hydro_compare_key(const void *a, const void *b);

int data_index_compare(const void *a, const void *b);
int peano_compare_key(const void *a, const void *b);

void mysort_dataindex(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_domain(void *b, size_t n, size_t s);
void mysort_idlist(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_pmnonperiodic(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_peano(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));

void check_wind_creation(void);
void treat_outflowing_particles(void);
void set_injection_accel(void);


int density_isactive(int n);
int GasGrad_isactive(int i);

#ifdef HYDRO_VOLUME_CORRECTIONS
void cellcorrections_calc(void);
void cellcorrections_final_operations_and_cleanup(void);
#endif

size_t sizemax(size_t a, size_t b);

void reconstruct_timebins(void);
void init_peano_map(void);
peanokey peano_hilbert_key(int x, int y, int z, int bits);
peanokey peano_and_morton_key(int x, int y, int z, int bits, peanokey *morton);
peanokey morton_key(int x, int y, int z, int bits);

void catch_abort(int sig);
void catch_fatal(int sig);
void terminate_processes(void);
void enable_core_dumps_and_fpu_exceptions(void);
void write_pid_file(void);
#ifdef PAUSE_RUN_TO_ATTACH_DEBUGGER
void pause_run_to_attach_debugger();
#endif
void pm_init_periodic_allocate(void);
void pm_init_periodic_free(void);
void move_particles(integertime time1);
void find_next_sync_point_and_drift(void);
void find_dt_displacement_constraint(double hfac);
#ifdef WAKEUP
void process_wake_ups(void);
#endif
void set_units_sfr(void);
void gravity_forcetest(void);
void allocate_commbuffers(void);
void allocate_memory(void);
void begrun(void);
void check_omega(void);
void compute_accelerations(void);
void compute_global_quantities_of_system(void);
void compute_potential(void);
void construct_timetree(void);
void star_formation_parent_routine(void);
#if defined(TURB_DRIVING)
void do_turb_driving_step_first_half(void);
void do_turb_driving_step_second_half(void);
double st_return_dt_between_updates(void);
double st_return_mode_correlation_time(void);
double st_return_rms_acceleration(void);
double st_turbdrive_get_gaussian_random_variable(void);
void st_turbdrive_init_ouseq(void);
void st_turbdrive_calc_phases(void);
double st_return_driving_scale(void);
#endif
double evaluate_NH_from_GradRho(MyFloat gradrho[3], double hsml, double rho, double numngb_ndim, double include_h, int target);


#ifdef GALSF
int is_particle_single_star_eligible(long i);
double evaluate_stellar_age_Gyr(double stellar_tform);
double evaluate_light_to_mass_ratio(double stellar_age_in_gyr, int i);
double calculate_relative_light_to_mass_ratio_from_imf(double stellar_age_in_gyr, int i);
double calculate_individual_stellar_luminosity(double mdot, double mass, long i);
double return_probability_of_this_forming_bh_from_seed_model(int i);

// this structure needs to be defined here, because routines for feedback event rates, etc, are shared among files //
struct addFB_evaluate_data_in_
{
    MyDouble Pos[3], Vel[3], Msne, unit_mom_SNe;
    MyFloat Hsml, V_i, SNe_v_ejecta;
#ifdef GALSF_FB_MECHANICAL
    MyFloat Area_weighted_sum[AREA_WEIGHTED_SUM_ELEMENTS];
#endif
#ifdef METALS
    MyDouble yields[NUM_METAL_SPECIES];
#endif
    int NodeList[NODELISTLENGTH];
}
*addFB_evaluate_DataIn_, *addFB_evaluate_DataGet_;

void particle2in_addFB_fromstars(struct addFB_evaluate_data_in_ *in, int i, int fb_loop_iteration);
double mechanical_fb_calculate_eventrates(int i, double dt);
double Z_for_stellar_evol(int i);
#ifdef METALS
void get_jet_yields(double *yields, int i);
#endif
#endif

#ifdef SINGLE_STAR_FB_JETS
double single_star_jet_velocity(int n);
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
double single_star_feedback_velocity_fortimestep(int n);
#endif

#ifdef GRAIN_FLUID
void apply_grain_dragforce(void);
#endif

#ifdef RT_INFRARED
double get_min_allowed_dustIRrad_temperature(void);
double get_rt_ir_lambdadust_effective(double T, double rho, double *nH0_guess, double *ne_guess, int target, int update_Tdust);
#endif

#if defined(FLAG_NOT_IN_PUBLIC_CODE) || (defined(RT_CHEM_PHOTOION) && defined(GALSF))
double particle_ionizing_luminosity_in_cgs(long i);
#endif



#ifdef GALSF_FB_MECHANICAL
void determine_where_SNe_occur(void);
void mechanical_fb_calc(int fb_loop_iteration);
void mechanical_fb_calc_toplevel(void);
void verify_and_assign_local_mechfb_integrals(void);
#endif

#ifdef GALSF_FB_THERMAL
void determine_where_addthermalFB_events_occur(void);
void thermal_fb_calc(void);
#endif

#ifdef COOL_METAL_LINES_BY_SPECIES
/*double GetMetalLambda(double, double);*/
double getSpCoolTableVal(long i,long j,long k,long tblK);
#ifndef CHIMES
double GetCoolingRateWSpecies(double nHcgs, double logT, double *Z);
double GetLambdaSpecies(long k_index, long index_x0y0, long index_x0y1, long index_x1y0, long index_x1y1, double dx, double dy, double dz, double mdz);
void LoadMultiSpeciesTables(void);
void ReadMultiSpeciesTables(int iT);
char *GetMultiSpeciesFilename(int i, int hk);
#endif
#endif

double bh_angleweight(double bh_lum_input, MyFloat bh_angle[3], double dx, double dy, double dz);
double bh_angleweight_localcoupling(int j, double cos_theta, double r, double H_bh);

#if defined(GALSF_SUBGRID_WINDS)
void assign_wind_kick_from_sf_routine(int i, double sm, double dtime, double* pvtau_return);
#endif


#if defined(BLACK_HOLES)
int blackhole_evaluate_PREPASS(int target, int mode, int *nexport, int *nSend_local);
#endif

#ifdef GALSF_SUBGRID_WINDS
#if (GALSF_SUBGRID_WIND_SCALING==2)
void disp_setup_smoothinglengths(void);
void disp_density(void);
#endif
#endif




#ifdef CHIMES
double chimes_convert_u_to_temp(double u, double rho, int target);
void chimes_update_gas_vars(int target);
void chimes_gizmo_exit(void);
#ifdef COOL_METAL_LINES_BY_SPECIES
void chimes_update_element_abundances(int i);
#endif
#ifdef CHIMES_TURB_DIFF_IONS
void chimes_update_turbulent_abundances(int i, int mode);
#endif
#ifdef CHIMES_METAL_DEPLETION
void chimes_init_depletion_data(void);
double chimes_jenkins_linear_fit(double nH, double T, double Ax, double Bx, double zx);
void chimes_compute_depletions(double nH, double T, int thread_id);
#endif
#else
#endif
void cooling_parent_routine(void);
void count_hot_phase(void);
void delete_node(int i);
void density(void);
void density_decouple(void);
void determine_interior(void);
int dissolvegas(void);
void do_box_wrapping(void);
double enclosed_mass(double R);
void energy_statistics(void);
void ensure_neighbours(void);

void output_log_messages(void);
void ewald_corr(double dx, double dy, double dz, double *fper);
void ewald_force(int ii, int jj, int kk, double x[3], double force[3]);
void ewald_force_ni(int iii, int jjj, int kkk, double x[3], double force[3]);
void ewald_init(void);
double ewald_psi(double x[3]);
double ewald_pot_corr(double dx, double dy, double dz);
int find_ancestor(int i);
integertime find_next_outputtime(integertime time);
void find_next_time(void);
integertime find_next_time_walk(int node);
void free_memory(void);
void advance_and_find_timesteps(void);
integertime get_timestep(int p, double *a, int flag);

void determine_PMinterior(void);
void gravity_tree(void);
void hydro_force(void);
void init(void);
void do_the_cooling_for_particle(int i);
double get_equilibrium_dust_temperature_estimate(int i, double shielding_factor_for_exgalbg);
double return_electron_fraction_from_heavy_ions(int target, double temperature, double density_cgs, double n_elec_HHe);
void apply_pm_hires_region_clipping_selection(int i);
double get_starformation_rate(int i, int mode);
void update_internalenergy_for_galsf_effective_eos(int i, double tcool, double tsfr, double cloudmass_fraction, double rateOfSF);
void init_clouds(void);
void integrate_sfr(void);
void insert_node(int i);
int mark_targets(void);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
void mpi_printf(const char *fmt, ...);
void open_outputfiles(void);
void write_outputfiles_header(void);
void peano_hilbert_order(void);
double pot_integrand(double xx);
void predict(double time);
void predict_collisionless_only(double time);
void predict_sph_particles(double time);
void prepare_decouple(void);
void read_ic(char *fname);
int read_outputlist(char *fname);
void read_parameter_file(char *fname);
void rearrange_particle_sequence(void);
void swap_treewalk_pointers(int i, int j);
void remove_particle_from_tree(int i);
void reorder_gas(void);
void reorder_particles(void);
void restart(int modus);
void run(void);
void savepositions(int num);
void savepositions_ioformat1(int num);
double my_second(void);
void set_softenings(void);
void set_sph_kernel(void);
void set_units(void);
void setup_smoothinglengths(void);
void apply_special_boundary_conditions(int i, double mass_for_dp, int mode);

void minimum_large_ints(int n, long long *src, long long *res);
void sumup_large_ints(int n, int *src, long long *res);
void sumup_longs(int n, long long *src, long long *res);

void statistics(void);
double timediff(double t0, double t1);
void veldisp(void);
void veldisp_ensure_neighbours(int mode);
int binarySearch(const double * arr, const double x, const int l, const int r, const int total);

double get_gravkick_factor(integertime time0, integertime time1);
double drift_integ(double a, void *param);
double gravkick_integ(double a, void *param);
double growthfactor_integ(double a, void *param);
double hydrokick_integ(double a, void *param);
void init_drift_table(void);
double get_drift_factor(integertime time0, integertime time1);
double measure_time(void);
double report_time(void);

/* on some DEC Alphas, the correct prototype for pow() is missing,
   even when math.h is included ! */

double pow(double, double);


void long_range_init(void);
void long_range_force(void);
void pm_init_periodic(void);
void pmforce_periodic(int mode, int *typelist);
void pm_init_regionsize(void);
void pm_init_nonperiodic(void);
int pmforce_nonperiodic(int grnr);

int pmpotential_nonperiodic(int grnr);
void pmpotential_periodic(void);

void readjust_timebase(double TimeMax_old, double TimeMax_new);

double enclosed_mass(double R);
void pm_setup_nonperiodic_kernel(void);


#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE)
#ifdef CHIMES_STELLAR_FLUXES
double chimes_G0_luminosity(double stellar_age, double stellar_mass);
double chimes_ion_luminosity(double stellar_age, double stellar_mass);
int rt_get_source_luminosity_chimes(int i, int mode, double *lum, double *chimes_lum_G0, double *chimes_lum_ion);
#endif
int rt_get_source_luminosity(int i, int mode, double *lum);
int rt_get_donation_target_bin(int bin);
int rt_get_lum_band_stellarpopulation(int i, int mode, double *lum);
int rt_get_lum_band_agn(int i, int mode, double *lum);
int rt_get_lum_band_singlestar(int i, int mode, double *lum);
void eddington_tensor_dot_vector(double ET[6], double vec_in[3], double vec_out[3]);
double return_flux_limiter(int target, int k_freq);
double rt_kappa(int j, int k_freq);
int check_if_absorbed_photons_can_be_reemitted_into_same_band(int kfreq);
double rt_absorb_frac_albedo(int j, int k_freq);
double rt_absorption_rate(int i, int k_freq);
double rt_diffusion_coefficient(int i, int k_freq);
void rt_eddington_update_calculation(int j);
void rt_update_driftkick(int i, double dt_entr, int mode);
#endif
#ifdef RT_SOURCE_INJECTION
void rt_source_injection(void);
#endif

#ifdef RADTRANSFER
void rt_set_simple_inits(int RestartFlag);
#if defined(RT_EVOLVE_INTENSITIES)
void rt_init_intensity_directions(void);
#endif
void rt_get_lum_gas(int target, double *je);
#ifdef RT_ISRF_BACKGROUND
void rt_apply_boundary_conditions(int i);
void get_background_isrf_urad(int i, double *urad);
#endif
double slab_averaging_function(double x);
double blackbody_lum_frac(double E_lower, double E_upper, double T_eff);
double stellar_lum_in_band(int i, double E_lower, double E_upper);

#ifdef RT_DIFFUSION_CG
void rt_diffusion_cg_solve(void);
#endif

#ifdef RT_CHEM_PHOTOION
double rt_return_photon_number_density(int i, int k);
double rt_photoion_chem_return_temperature(int i, double internal_energy);
void rt_update_chemistry(void);
void rt_get_sigma(void);
double rt_GetCoolingTime(int i, double u, double rho);
double rt_cooling_photoheating(int i, double dt);
double rt_DoCooling(int i, double dt);
double rt_DoHeating(int i, double dt);
double rt_get_cooling_rate(int i, double internal_energy);
void rt_write_chemistry_stats(void);
#endif

#endif


void find_block(char *label,FILE *fd);

#ifdef GDE_DISTORTIONTENSOR
void get_half_kick_distortion(int pindex, MyBigFloat half_kick_add[6][6]);
void analyse_phase_space(int pindex,
                         MyBigFloat *s_1, MyBigFloat *s_2, MyBigFloat *s_3,
                         MyBigFloat *smear_x, MyBigFloat *smear_y, MyBigFloat *smear_z,
                         MyBigFloat *second_deriv, MyBigFloat *sigma);
MyBigFloat get_analytic_annihilation(MyBigFloat s_1, MyBigFloat s_2, MyBigFloat s_3,
                                   MyBigFloat s_1_prime, MyBigFloat s_2_prime, MyBigFloat s_3_prime,
                                   MyBigFloat second_deriv, MyBigFloat second_deriv_prime,
                                   MyBigFloat dt, MyBigFloat sigma);
MyBigFloat get_max_caustic_density(MyBigFloat s_1, MyBigFloat s_2,
                                 MyBigFloat s_1_prime, MyBigFloat s_2_prime,
                                 MyBigFloat second_deriv, MyBigFloat second_deriv_prime,
                                 MyBigFloat sigma,
                                 int pindex);
void get_current_ps_info(int pindex, MyBigFloat * flde, MyBigFloat * psde);
void do_phase_space_drift(int i, double dt_drift);
void do_the_phase_space_kick(int i, double dt_gravkick);
void do_long_range_phase_space_kick(int i, double dt_gravkick);
/* some math functions we need from phasespace_math.c */
void ludcmp(MyBigFloat ** a, int n, int *indx, MyBigFloat * d);
MyBigFloat **matrix(long nrl, long nrh, long ncl, long nch);
MyBigFloat *vector(long nl, long nh);
void free_matrix(MyBigFloat ** m, long nrl, long nrh, long ncl, long nch);
void free_vector(MyBigFloat * v, long nl, long nh);
void mult_matrix(MyBigFloat ** matrix_a, MyBigFloat ** matrix_b, int dimension, MyBigFloat ** matrix_result);
void mult_matrix_transpose_A(MyBigFloat ** matrix_a, MyBigFloat ** matrix_b, int dimension, MyBigFloat ** matrix_result);
void luinvert(MyBigFloat ** input_matrix, int n, MyBigFloat ** inverse_matrix);
void eigsrt(MyBigFloat d[], MyBigFloat ** v, int n);
void jacobi(MyBigFloat ** a, int n, MyBigFloat d[], MyBigFloat ** v, int *nrot);
#endif


#if defined(COMPUTE_TIDAL_TENSOR_IN_GRAVTREE)
#ifdef PMGRID
void pmtidaltensor_periodic_diff(void);
int pmtidaltensor_nonperiodic_diff(int grnr);
void pmtidaltensor_periodic_fourier(int component);
int pmtidaltensor_nonperiodic_fourier(int component, int grnr);
#endif
#endif

int ags_gravity_kernel_shared_BITFLAG(short int particle_type_primary);
#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
void ags_setup_smoothinglengths(void);
void ags_density(void);
int ags_density_isactive(int i);
double ags_return_maxsoft(int i);
double ags_return_minsoft(int i);
void AGSForce_calc(void);
#endif

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
void advect_mesh_point(int i, double dt);
double calculate_face_area_for_cartesian_mesh(double *dp, double rinv, double l_side, double *Face_Area_Vec);
#endif

#if (SINGLE_STAR_TIMESTEPPING > 0)
void subtract_companion_gravity(int i);
#endif


void hydro_gradient_calc(void);
int GasGrad_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int gradient_iteration);
void construct_gradient(double *grad, int i);
void local_slopelimiter(double *grad, double valmax, double valmin, double alim, double h, double shoot_tol, int pos_preserve, double d_max, double val_cen);

#ifdef TURB_DIFF_DYNAMIC
void dynamic_diff_vel_calc(void);
void dynamic_diff_calc(void);
int DynamicDiff_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int dynamic_iteration);
void *DynamicDiff_evaluate_primary(void *p, int dynamic_iteration);
void *DynamicDiff_evaluate_secondary(void *p, int dynamic_iteration);
#endif

#ifdef PARTICLE_EXCISION
void apply_excision();
#endif

#ifdef DM_SIDM
double prob_of_interaction(double mass, double r, double h_si, double dV[3], double dt);
double g_geo(double r);
void calculate_interact_kick(double dV[3], double kick[3], double m);
void init_geofactor_table(void);
double geofactor_integ(double x, void * params);
double geofactor_angle_integ(double u, void * params);
void init_self_interactions();
#ifdef GRAIN_COLLISIONS
double return_grain_cross_section_per_unit_mass(int i);
double prob_of_grain_interaction(double cx_per_unitmass, double mass, double r, double h_si, double dV[3], double dt, int j_ngb);
#endif
#endif



#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
double do_cbe_nvt_inversion_for_faces(int i);
#endif

#ifdef DM_FUZZY
void do_dm_fuzzy_initialization(void);
void do_dm_fuzzy_drift_kick(int pindex, double dt_entr, int mode);
void DMGrad_gradient_calc(void);
void do_dm_fuzzy_flux_computation(double HLLwt, double dt, double prev_a, double dv[3],
                                  double GradRho_L[3], double GradRho_R[3],
                                  double GradRho2_L[3][3], double GradRho2_R[3][3],
                                  double rho_L, double rho_R, double dv_Right_minus_Left,
                                  double Area[3], double fluxes[3], double AGS_Numerical_QuantumPotential_L, double AGS_Numerical_QuantumPotential_R, double *dt_egy_Numerical_QuantumPotential);
void do_dm_fuzzy_flux_computation_old(double HLLwt, double dt, double m0, double prev_a, double dp[3], double dv[3],
                                  double GradRho_L[3], double GradRho_R[3],
                                  double GradRho2_L[3][3], double GradRho2_R[3][3],
                                  double rho_L, double rho_R, double dv_Right_minus_Left,
                                  double Area[3], double fluxes[3], double AGS_Numerical_QuantumPotential, double *dt_egy_Numerical_QuantumPotential);
void dm_fuzzy_reconstruct_and_slopelimit_sub(double *u_R_f, double *u_L_f, double q_R, double dq_R_0[3], double q_L, double dq_L_0[3], double dx[3]);
void dm_fuzzy_reconstruct_and_slopelimit(double *u_R, double du_R[3], double *u_L, double du_L[3],
                                         double q_R, double dq_R[3], double d2q_R[3][3],
                                         double q_L, double dq_L[3], double d2q_L[3][3],
                                         double dx[3]);
#endif


#ifdef SINGLE_STAR_TIMESTEPPING
void kepler_timestep(int i, double dt, double kick_dv[3], double drift_dx[3], int mode);
void odeint_super_timestep(int i, double dt_super, double kick_dv[3], double drift_dx[3], int mode);
double gravfac(double r, double mass);
double gravfac2(double r, double mass);
void grav_accel_jerk(double mass, double dx[3], double dv[3], double accel[3], double jerk[3]);
double eccentric_anomaly(double mean_anomaly, double ecc);
#endif
