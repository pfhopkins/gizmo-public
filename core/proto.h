#ifndef ALLVARS_H
#include "../declarations/allvars.h"
#endif
#include "../gravity/forcetree.h"
#include "../domain/domain.h"
#ifdef COOLING
#include "../cooling/cooling.h"
#endif
#ifdef SINK_PARTICLES
#include "../sinks/sink.h"
#endif

/*! declarations of functions throughout the code */
/*!
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * (adding/removing routines as necessary)
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


void write_header_attributes_in_hdf5(hid_t handle);
void read_header_attributes_in_hdf5(char *fname);

void output_compile_time_options(void);

void report_pinning(void);
void pin_to_core_set(void);
void get_core_set(void);

void init_turb(void);
void set_turb_ampl(void);
int new_turbforce_needed_this_timestep(void);
void add_turb_accel(void);
void log_turb_temp(void);
double get_random_number(MyIDType id);
double mpi_report_comittable_memory(long long BaseMem, int verbose);
long long report_comittable_memory(long long *MemTotal,
                                   long long *Committed_AS,
                                   long long *SwapTotal,
                                   long long *SwapFree);

void merge_and_split_particles(void);
int does_particle_need_to_be_merged(int i);
int does_particle_need_to_be_split(int i);
double target_mass_renormalization_factor_for_mergesplit(int i, int split_key);
#if defined(FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT) || defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
int check_if_sufficient_mergesplit_time_has_passed(int i);
int is_particle_a_special_zoom_target(int i);
#endif
int merge_particles_ij(int i, int j);
int split_particle_i(int i, int n_particles_split, int i_nearest);
void do_first_halfstep_kick(void);
void do_second_halfstep_kick(void);
double matrix_invert_ndims(Mat3<double>& T, Mat3<double>& Tinv);
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

void calculate_and_assign_nonideal_mhd_coefficients(int i, struct particle_data *pp, struct gas_cell_data *cell);
void calculate_and_assign_conduction_and_viscosity_coefficients(int i, struct particle_data *pp, struct gas_cell_data *cell);
#ifdef TURB_DIFFUSION
void calculate_and_assign_turbulent_diffusion_coefficients(int i, struct particle_data *pp = P, struct gas_cell_data *cell = CellP);
#endif

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
void drift_extra_physics(int i, integertime tstart, integertime tend, double dt_entr);

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

/* DMAX/DMIN/IMAX/IMIN are now macros in macros.h (GPU-safe, no linkage issues). */
GIZMO_GPU_FUNCTION static inline double MINMOD(double a, double b) {return (a>0) ? ((b<0) ? 0 : DMIN(a,b)) : ((b>=0) ? 0 : DMAX(a,b));}
/* special version of MINMOD below: a is always the "preferred" choice, b the stability-required one. here we allow overshoot, just not opposite signage */
GIZMO_GPU_FUNCTION static inline double MINMOD_G(double a, double b) {return a;}

static inline double c_light_code_reduced(int k_freq, struct particle_data *pp, struct gas_cell_data *cell) {
#if defined(RT_FLUXLIMITER)
    return C_LIGHT_CODE;
#else
    return RT_SPEEDOFLIGHT_REDUCTION * C_LIGHT_CODE;
#endif
}
static inline double rsol_correction_factor_for_velocity_terms(int k_freq, struct particle_data *pp, struct gas_cell_data *cell) {return RT_SPEEDOFLIGHT_REDUCTION;}

double ForceSoftening_KernelRadius(int p);
double compute_force_softening_kernel_radius(int p); /* internal: source-of-truth softening computation */
void   compute_all_force_softening(int mode);        /* mode=0 active-only; mode=1 all NumPart (init/restart) */
GIZMO_GPU_FUNCTION inline double sigmoid_sqrt(double x) {return 0.5*(1 + x/sqrt(1+x*x));} /* inline for GPU single-TU; definition also in global.cc for non-GPU TUs */
/* velocity_gradient_norm is now a member function of gas_cell_data — use cell[i].velocity_gradient_norm() */

#ifdef BOX_SHEARING
void calc_shearing_box_pos_offset(void);
#endif


/* ngb_treefind_*_targeted decls retired with mesh/ngb.cc in Step 5 Phase D5. */



void set_predicted_quantities_for_extra_physics(int i);
void do_kick_for_extra_physics(int i, integertime tstart, integertime tend, double dt_entr);
#if (SINGLE_STAR_TIMESTEPPING > 0)
#endif

void set_eos_pressure(int i, struct particle_data *pp = P, struct gas_cell_data *cell = CellP);
double return_user_desired_target_density(int i);
double return_user_desired_target_pressure(int i);
#ifdef EOS_TILLOTSON
void tillotson_eos_init(void);
#endif
#ifdef EOS_ANEOS
#include "../eos/aneos.h"
#endif

#ifdef SPECIAL_POINT_WEIGHTED_MOTION
double weight_function_for_weighted_motion_smoothing(double r, int mode);
#endif

#ifdef EOS_SUBSTELLAR_ISM
void hydrogen_molecule_partitionfunc(double temp, double result[3]);
void hydrogen_molecule_zrot_mixture(double temp, double result[3]);
void hydrogen_molecule_zvib(double temp, double result[3]);
double hydrogen_molecule_energy(double temp);
double hydrogen_molecule_gamma(double temp);
#endif

void read_fof(int num);
int fof_compare_ID_list_ID(const void *a, const void *b);

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


int MPI_Check_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                       int dest, int sendtag, void *recvbufreal, int recvcount,
                       MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status * status);

int MPI_Sizelimited_Sendrecv(void *sendbuf0, size_t sendcount, MPI_Datatype sendtype,
                             int dest, int sendtag, void *recvbuf0, size_t recvcount,
                             MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
                             MPI_Status *status);

int mpi_calculate_offsets(int *send_count, int *send_offset, int *recv_count, int *recv_offset, int send_identical);
void sort_based_on_field(void *data, int field_offset, int n_items, int item_size, void **data2ptr);
void mpi_distribute_items_to_tasks(void *data, int task_offset, int *n_items, int *max_n, int item_size);
int getNodeCount(void);

void parallel_sort_special_P_GrNr_ID(void);
void calculate_power_spectra(int num, long long *ntot_type_all);


int pmforce_is_particle_high_res(int type, Vec3<double>& Pos);

void assign_unique_ids(void);

void conduction(void);
int conduction_evaluate(int target, int mode, double *in, double *out, double *sum,
			int *nexport, int *nsend_local);


void fof_compute_group_properties(int gr, int start, int len);

#ifdef TURB_DIFF_DYNAMIC
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
void find_particles_and_save_them(int num);
void lineofsight_output(void);
void sum_over_processors_and_normalize(void);
void absorb_along_lines_of_sight(void);
void output_lines_of_sight(int num);
integertime find_next_lineofsighttime(integertime time0);
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

/* Get_Particle_Size is now a member function of particle_data — use P[i].Get_Particle_Size() or pp[i].Get_Particle_Size() */
double INLINE_FUNC Get_Particle_Expected_Area(double h);
double Get_Gas_Ionized_Fraction(int i, struct particle_data *pp = P, struct gas_cell_data *cell = CellP);
double CR_calculate_adiabatic_gasCR_exchange_term(int i, double dt_entr, double gamma_minus_eCR_tmp, int mode, struct particle_data *pp, struct gas_cell_data *cell);
double INLINE_FUNC Get_CosmicRayEnergyDensity_cgs(int i, struct particle_data *pp, struct gas_cell_data *cell);
double CR_gas_heating(int target, double n_elec, double nH0, double nHcgs, struct particle_data *pp, struct gas_cell_data *cell);
double Get_CosmicRayIonizationRate_cgs(int i, struct particle_data *pp, struct gas_cell_data *cell);
#ifdef COSMIC_RAY_FLUID
/* CR inline functions (Get_Gas_CosmicRayPressure, cosmicrayfluid_rsol_corrfac, etc.) are in
   cosmic_ray_functions.h. Included by gradient_functions.h, hydro_functions.h, and
   cooling_functions.h — NOT here, because cosmic_ray_functions.h calls functions
   declared later in this file (CR_return_mean_rigidity_in_bin_in_GV, etc.). */
double return_CRbin_M1speed(int k_CRegy);
double INLINE_FUNC cosmicrayfluid_rsol_corrfac(int k);
double INLINE_FUNC Get_Gas_CosmicRayPressure(int i, int k_CRegy, struct gas_cell_data *cell);
double gamma_eos_of_crs_in_bin(int k_CRegy);
void CalculateAndAssign_CosmicRay_DiffusionAndStreamingCoefficients(int i, struct particle_data *pp, struct gas_cell_data *cell);
double Get_CosmicRayGradientLength(int i, int k_CRegy, struct particle_data *pp, struct gas_cell_data *cell);
double CosmicRay_Update_DriftKick(int i, double dt_entr, int mode, struct particle_data *pp, struct gas_cell_data *cell);
double CR_energy_spectrum_injection_fraction(int k_CRegy, int source_type, double shock_vel, int return_index_in_bin, int target, struct particle_data *pp, struct gas_cell_data *cell);
double return_cosmic_ray_anisotropic_closure_function_threechi(int target, int k_CRegy, struct gas_cell_data *cell);
void inject_cosmic_rays(double CR_energy_to_inject, double injection_velocity, int source_type, int target, double *dir, struct gas_cell_data *cell);
double evaluate_cr_transport_reductionfactor(int target, int k_CRegy, int mode, struct gas_cell_data *cell);
double Get_AlfvenMachNumber_Local(int i, double vA_idealMHD_codeunits, int use_shear_corrected_vturb_flag, struct gas_cell_data *cell);
double diffusion_coefficient_constant(int target, int k_CRegy, struct gas_cell_data *cell);
double diffusion_coefficient_extrinsic_turbulence(int mode, int target, int k_CRegy, double M_A, double L_scale, double b_muG, double vA_noion, double rho_cgs, double temperature, double cs_thermal, double nh0, double nHe0, double f_ion);
double diffusion_coefficient_self_confinement(int mode, int target, int k_CRegy, double M_A, double L_scale, double b_muG, double vA_noion, double rho_cgs, double temperature, double cs_thermal, double nh0, double nHe0, double f_ion);
double return_CRbin_CR_charge_in_e(int target, int k_CRegy);
int return_CRbin_CR_species_ID(int k_CRegy);
void CR_cooling_and_losses(int target, double n_elec, double nHcgs, double dtime_cgs, struct particle_data *pp, struct gas_cell_data *cell);
double return_CRbin_CRmass_in_mp(int target, int k_CRegy);
double CR_get_streaming_loss_rate_coefficient(int target, int k_CRegy, struct particle_data *pp, struct gas_cell_data *cell);
double Get_Gas_ion_Alfven_speed_i(int i, struct particle_data *pp, struct gas_cell_data *cell);
double return_CRbin_nuplusminus_asymmetry(int i, int k_CRegy, struct gas_cell_data *cell);
#if defined(CRFLUID_EVOLVE_SPECTRUM)
void CR_spectrum_define_bins(void);
void CR_initialize_multibin_quantities(void);
void CR_cooling_and_losses_multibin(int target, double n_elec, double nHcgs, double dtime_cgs, int mode_driftkick, struct particle_data *pp, struct gas_cell_data *cell);
double CR_return_slope_from_number_and_energy_in_bin(double energy_in_code_units, double number_effective_in_code_units, double bin_centered_energy_in_GeV, int k_bin);
double CR_return_new_bin_edge_from_rate(double rate_dt_dimless, double x_m_bin, double x_p_bin, int loss_mode, int NR_key, double additional_variable_dummy);
double CR_coulomb_energy_integrand(double x, double tau, double slope);
double CR_reaccel_energy_integrand(double x, double tau, double slope, double delta_slope, int NR_key);
double CR_compton_energy_integrand(double x, double tau, double slope);
int CR_check_if_bin_is_nonrelativistic(int k_bin);
double CR_return_true_number_in_bin(int target, int k_bin, struct gas_cell_data *cell);
double CR_return_effective_number_in_bin_in_codeunits(int target, int k_bin, struct gas_cell_data *cell);
double CR_return_spectral_slope_target(int target, int k_bin, struct gas_cell_data *cell);
double CR_get_number_in_bin_from_slope(int target, int k_bin, double energy, double slope);
double CR_return_mean_energy_in_bin_in_GeV(int target, int k_bin, struct gas_cell_data *cell);
double CR_return_mean_rigidity_in_bin_in_GV(int target, int k_bin, struct gas_cell_data *cell);
int compare_CR_rigidity_for_sort(const void *a, const void *b);
double return_CRbin_kinetic_energy_in_GeV_binvalsNRR(int k_CRegy);
#endif
#endif
#ifdef EOS_ELASTIC
void elastic_body_update_driftkick(int i, double dt_entr, int mode);
#endif
#if defined(EOS_ELASTIC) || defined(EOS_TILLOTSON) || defined(EOS_ANEOS)
double get_negative_pressure_tensilecorrfac(double r, double h_i, double h_j);
#endif
/* Get_Gas_effective_soundspeed_i is now cell[i].effective_soundspeed() */
/* Get_Gas_Fast_MHD_wavespeed_i is now cell[i].fast_MHD_wavespeed() */
double Get_Gas_Mean_Molecular_Weight_mu(double T_guess, double rho, double *xH0, double *ne_guess, double urad_from_uvb_in_G0, int target, struct particle_data *pp, struct gas_cell_data *cell);
void update_explicit_molecular_fraction(int i, double dtime_cgs, struct particle_data *pp, struct gas_cell_data *cell);
double molecfrac_rootfind_function(double fH2, double x00, double x01, double x_b_0, double x_c, double y_a, double G_LW_dt_unshielded);
double return_dust_to_metals_ratio_vs_solar(int i, double T_dust_manual_override, struct particle_data *pp, struct gas_cell_data *cell);
double INLINE_FUNC yhelium(int target, struct particle_data *pp);
double Get_Gas_Molecular_Mass_Fraction(int i, double temperature, double neutral_fraction, double free_electron_ratio, double urad_from_uvb_in_G0, struct particle_data *pp, struct gas_cell_data *cell);
#ifdef MAGNETIC
double Get_DtB_FaceArea_Limiter(int i, struct particle_data *pp, struct gas_cell_data *cell);
#ifdef DIVBCLEANING_DEDNER
double INLINE_FUNC Get_Gas_PhiField(int i_particle_id);
double INLINE_FUNC Get_Gas_PhiField_DampingTimeInv(int i_particle_id);
#endif
#endif
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
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

void sink_accretion(void);
#ifdef SINK_WIND_SPAWN
void get_random_orthonormal_basis(int seed, Vec3<double>& nx, Vec3<double>& ny, Vec3<double>& nz);
void get_wind_spawn_direction(int i, int num_spawned_this_call, int mode, Vec3<double>& ny, Vec3<double>& nz, Vec3<double>& veldir, Vec3<double>& dpdir);
double get_spawned_cell_launch_speed(int i, struct particle_data *pp);
#ifdef MAGNETIC
void get_wind_spawn_magnetic_field(int j, int mode, Vec3<double>& ny, Vec3<double>& nz, Vec3<double>& dpdir, double d_r);
#endif
int sink_spawn_particle_wind_shell( int i, int dummy_cell_i_to_clone, int num_already_spawned );
void spawn_sink_wind_feedback(void);
#endif
int sink_evaluate(int target, int mode, int *nexport, int *nsend_local);
int sink_evaluate_swallow(int target, int mode, int *nexport, int *nsend_local);


void fof_fof(int num);
void fof_compile_catalogue(void);
void fof_save_groups(int num);
void fof_save_local_catalogue(int num);

void fof_make_sink_particles(void);

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

#ifdef GALSF_MERGER_STARCLUSTER_PARTICLES
double evaluate_starstar_merger_for_starcluster_particle_pair(int i, int j);
int evaluate_starstar_merger_for_starcluster_eligibility(int i);
#endif


void long_range_init_regionsize(void);
int find_files(char *fname);
void pm_init_nonperiodic_allocate(void);
void  pm_init_nonperiodic_free(void);
int data_index_compare(const void *a, const void *b);
int peano_compare_key(const void *a, const void *b);
void mysort_dataindex(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_domain(void *b, size_t n, size_t s);
void mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_pmnonperiodic(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_peano(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));

int density_isactive(int n);
int GasGrad_isactive(int i, struct particle_data *pp, struct gas_cell_data *cell);

#ifdef HYDRO_VOLUME_CORRECTIONS
void cellcorrections_calc(void);
void cellcorrections_final_operations_and_cleanup(void);
#endif

size_t sizemax(size_t a, size_t b);

void reconstruct_timebins(void);
peano1D domain_double_to_int(double d);
peanokey peano_hilbert_key(peano1D x, peano1D y, peano1D z, int bits);
peanokey peano_and_morton_key(peano1D x, peano1D y, peano1D z, int bits, peanokey *morton);
peanokey morton_key(peano1D x, peano1D y, peano1D z, int bits);

void pm_init_periodic_allocate(void);
void pm_init_periodic_free(void);
void move_particles(integertime time1);
void gizmo_full_drift_to(integertime time1);
#ifdef __cplusplus
extern "C" {
#endif
void gizmo_full_drift_invalidate(void);
#ifdef __cplusplus
}
#endif
void ghost_exchange(double safety_factor);
void ghost_exchange_hydro(double safety_factor);
void ghost_exchange_hydro_oneway(double safety_factor);
struct ghost_exchange_spec_t;
extern "C" void ghost_exchange_run(const struct ghost_exchange_spec_t *spec);
void ghost_exchange_cleanup(void);
int ghost_exchange_needs_redo(void);
/* Canonical out-of-line host accessors for `All.*` — defined in
 * declarations/allvars.cc (the TU that owns `All`). Useful where
 * explicit host-extern intent helps readability. With the central
 * AllDevice mirror's device-pass-gated redirect, host wrapper code in
 * GPU TUs may also read `All.*` directly. */
#ifdef __cplusplus
extern "C" {
#endif
struct global_data_all_processes *gizmo_host_all_ptr(void);
integertime gizmo_host_ti_current(void);
#ifdef __cplusplus
}
#endif
#ifdef __cplusplus
extern "C" {
#endif
void ghost_exchange_local_tree_invalidate_drift(void);
void ghost_exchange_local_tree_invalidate_full(void);
void ghost_exchange_local_tree_mark_h_dirty_indices(const int *indices, int n);
void ghost_exchange_local_tree_mark_h_dirty_range(int start, int end);
#ifdef __cplusplus
}
#endif
int ghost_get_previous_count(void);
void find_next_sync_point_and_drift(void);
void find_dt_displacement_constraint(double hfac);
void process_wake_ups(void);
void set_units_sfr(void);
void allocate_memory(void);
void begrun(void);
void check_omega(void);
void compute_global_quantities_of_system(void);
void star_formation_parent_routine(void);
#if defined(GRAIN_FLUID) && defined(GRAIN_FLUID_PROMOTION)
void grain_promotion_parent_routine(void);
#endif
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
GIZMO_GPU_FUNCTION double evaluate_NH_from_GradRho(MyFloat gradrho[3], double rkern, double rho, double numngb_ndim, double include_h, int target, struct particle_data *pp = P);
GIZMO_GPU_FUNCTION inline double evaluate_NH_from_GradRho(const Vec3<MyFloat>& gradrho, double rkern, double rho, double numngb_ndim, double include_h, int target, struct particle_data *pp = P) { return evaluate_NH_from_GradRho(const_cast<MyFloat*>(gradrho.data), rkern, rho, numngb_ndim, include_h, target, pp); }
double evaluate_time_since_t_initial_in_Gyr(double t_initial);

#ifdef GALSF
int is_particle_single_star_eligible(long i);
double evaluate_stellar_age_Gyr(long i);
double evaluate_light_to_mass_ratio(double stellar_age_in_gyr, int i);
double calculate_relative_light_to_mass_ratio_from_imf(double stellar_age_in_gyr, int i, int mode);
double calculate_individual_stellar_luminosity(double mdot, double mass, long i);
double return_probability_of_this_forming_sink_from_seed_model(int i);

void particle2in_addFB_fromstars(struct addFB_evaluate_data_in_ *in, int i, int fb_loop_iteration);
double mechanical_fb_calculate_eventrates(int i, double dt);
double Z_for_stellar_evol(int i);
#ifdef METALS
void get_jet_yields(double *yields, int i);
#endif
#if defined(GALSF_FB_MECHANICAL) && defined(GALSF_FB_FIRE_STELLAREVOLUTION)
double mechanical_fb_calculate_eventrates_SNe(int i, double dt);
void mechanical_fb_calculate_eventrates_Winds(int i, double dt);
void mechanical_fb_calculate_eventrates_Rprocess(int i, double dt);
void mechanical_fb_calculate_eventrates_Agetracers(int i, double dt);
void particle2in_addFB_SNe(struct addFB_evaluate_data_in_ *in, int i);
void particle2in_addFB_winds(struct addFB_evaluate_data_in_ *in, int i);
void particle2in_addFB_Rprocess(struct addFB_evaluate_data_in_ *in, int i);
void particle2in_addFB_ageTracer(struct addFB_evaluate_data_in_ *in, int i);
#ifdef METALS
void get_wind_yields(double *yields, int i);
void get_SNe_yields(double *yields,int i,double t_gyr,int SNeIaFlag, double *Msne);
#endif
#ifdef GALSF_FB_FIRE_AGE_TRACERS
#ifdef GALSF_FB_FIRE_AGE_TRACERS_CUSTOM
int read_agetracerlist(char *fname);
#endif
int get_age_tracer_bin(double age);
double get_age_tracer_bin_start_time(int k);
#endif
#endif
#endif


#if defined(GALSF_ISMDUSTCHEM_MODEL)
void Initialize_ISMDustChem_Global_Variables();
void Initialize_ISMDustChem_Particle_Variables(int i, struct particle_data *pp, struct gas_cell_data *cell);
void update_dust_processes(int i, double dtime_gyr, struct particle_data *pp, struct gas_cell_data *cell);
void ISMDustChem_get_SNe_dust_yields(double *yields, int i, double t_gyr, int SNeIaFlag, double Msne, struct particle_data *pp, struct gas_cell_data *cell);
void ISMDustChem_get_wind_dust_yields(double *yields, int i, struct gas_cell_data *cell);
double specific_Z_AGB_dust(int spec_indx, double star_age, int z_bound);
double cumulative_AGB_dust_returns(int dust_type, double star_age, double z);
void update_dense_molecular_fields(int i, double temp, double rho, double nh0, double ne, struct particle_data *pp, struct gas_cell_data *cell);
void update_dust_accretion(int i, double dtime_gyr, double temp, double rho, struct particle_data *pp, struct gas_cell_data *cell);
void update_dust_sputtering(int i, double dtime_gyr, double temp, double rho, struct particle_data *pp, struct gas_cell_data *cell);
double Lambda_Dust_HighTemperature_Gas_ISM(int target, double T, double n_elec, struct particle_data *pp, struct gas_cell_data *cell);
double return_ismdustchem_species_of_interest_for_diffusion_and_yields(int i, int k, double mass, struct gas_cell_data *cell);
double ISMDustChem_Return_Mass_Where_Dust_Shocked(double rho_cell_in_code_units, double Esne51_into_cell, double mass_preshock_in_code_units, double Z_cell);
void update_ISMDustChem_after_mechanical_injection(int j, double mass_shocked, double m0, double mf, double *Z_injected);
void ISMDustChem_update_iron_inclusions(int i, struct particle_data *pp, struct gas_cell_data *cell);
void ISMDustChem_get_elem_yields_from_species_yields(double *dust_yields, double *species_yields);
void ISMDustChem_get_species_key_elem(int spec_indx, double *dust_metallicity, int *key_elem, double *key_num_atoms, double *key_mass);
void ISMDustChem_get_species_properties(int spec_indx, double *dust_atomic_weight, double *bulk_dens);
void ISMDustChemEvo_renormalize_dust_fields(int i, struct particle_data *pp, struct gas_cell_data *cell);
void check_dust_fields(int i, int update_process);
#if defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
void Initialize_ISMDustChemEvo_Particle_Variables(int i, struct particle_data *pp, struct gas_cell_data *cell);
double get_ISMDustChemEvo_bin_mass(int i, int j, int k, struct gas_cell_data *cell);
void update_ISMDustChemEvo_bin_number_and_slope(int i, int j, int k, double number_in_bin, double mass_in_bin, struct gas_cell_data *cell);
void check_for_slope_limiting(int k, double bulk_dens, double *number_in_bin, double *slope_in_bin, double mass_in_bin);
void ISMDustChemEvo_get_SNe_dust_grain_size_yields(double *yields, int i, int SNeIaFlag, double Msne, struct particle_data *pp, struct gas_cell_data *cell);
void ISMDustChemEvo_get_wind_dust_grain_size_yields(double *yields, double Msne);
void ISMDustChemEvo_update_bins_given_grain_size_change(int i, int j, double *bin_da, double mass_limit, struct gas_cell_data *cell);
void update_dust_shattering_and_coagulation(int i, double dtime_gyr, double temp, double rho, struct particle_data *pp, struct gas_cell_data *cell);
void update_dust_photodestruction(int i, double dtime_gyr, struct gas_cell_data *cell);
double shattering_coagulation_polynomial(int i, int spec_indx, int bin_i, int bin_j, struct gas_cell_data *cell);
void ISMDustChemEvo_update_bins_given_mass_change(int i, int j, double *bin_dM, double bulk_dens, struct gas_cell_data *cell);
void ISMDustChemEvo_get_new_bin_N_and_slope_given_mass_change(double *bin_dM, double *bin_M, double *bin_N, double *bin_slope, double *new_bin_N, double *new_bin_slope, double bulk_dens);
void ISMDustChem_SNe_sputtering_step(int spec_indx, double *init_bin_N, double *init_bin_slope, double *init_bin_M, double *final_bin_N, double *final_bin_slope, double *final_bin_M, double bulk_dens);
void ISMDustChem_SNe_shattering_step(int spec_indx, double *init_bin_N, double *init_bin_slope, double *init_bin_M, double *final_bin_N, double *final_bin_slope, double *final_bin_M, double bulk_dens);
// Below functions only for debugging
void ISMDustChemEvo_check_bins_after_update(int i, int update_process, double mass); 
void ISMDustChemEvo_check_yields_before_update(double *bin_nums, double *bin_slopes, double *bin_masses, int yields_process, int species_num, double total_mass);
#endif
#endif

#if defined(GALSF_SFR_IMF_SAMPLING_DISTRIBUTE_SF)
void update_stellarnumber_and_timedistribofstarformation(void);
#endif

#ifdef SINGLE_STAR_FB_JETS
double single_star_jet_velocity(int n);
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
double single_star_feedback_velocity_fortimestep(int n);
#endif
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION
void singlestar_subgrid_protostellar_evolution_update_track(int n, double dm, double dt);
#if (SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION == 2)
double ps_polytropic_index(int stage, double mdot, double mass);
double ps_rhoc(double m, double n_ad, double r);
double ps_Pc(double m, double n_ad, double r);
double ps_Tc(double rhoc, double Pc);
double ps_beta(double m, double n_ad, double rhoc, double Pc);
double ps_betac(double rhoc, double Pc, double Tc);
double ps_dlogbeta_dlogm(double m, double r, double n_ad, double beta, double rhoc, double Pc);
double ps_dlogbetaperbetac_dlogm(double m, double r, double n_ad, double beta, double rhoc, double Pc, double Tc);
double ps_lum_I(double mdot);
double ps_lum_MS(double m);
double ps_radius_MS_in_solar(double m);
double ps_lum_Hayashi_BB(double m, double r);
#endif
double stellar_lifetime_in_Gyr(int n);
double single_star_relic_SN_mass(int n);
#if defined(SINGLE_STAR_FB_WINDS)
double single_star_wind_mdot(int n, int set_mode);
double single_star_wind_velocity(int n);
double single_star_escape_speed(int n);
double singlestar_WR_lifetime_Gyr(int n);
#endif
#if defined(SINGLE_STAR_FB_SNE)
double single_star_SN_velocity(int n);
void single_star_SN_init_directions(void);
#endif
#ifdef SINGLE_STAR_MERGE_AWAY_CLOSE_BINARIES
int is_star_eligible_for_binary_merge_away(int j);
#endif
#endif

#ifdef GRAIN_FLUID
void apply_grain_dragforce(void);
#endif

#ifdef RT_INFRARED
double get_min_allowed_dustIRrad_temperature(void);
double dust_dE_cooling(int i, double Tgas, double Tdust, double *Tdust_fixedpoint_1, double *Tdust_fixedpoint_2, struct particle_data *pp, struct gas_cell_data *cell);
double rt_ir_lambdadust(int i, double Tgas, struct particle_data *pp, struct gas_cell_data *cell);
#endif

#if defined(GALSF_FB_FIRE_RT_HIIHEATING) || (defined(RT_CHEM_PHOTOION) && defined(GALSF))
double particle_ionizing_luminosity_in_cgs(long i);
#endif

#ifdef GALSF_FB_FIRE_RT_HIIHEATING
void HII_heating_singledomain(void);
int do_the_local_ionization(int target, double dt, int source);
#endif

#ifdef GALSF_FB_FIRE_RT_LONGRANGE
void selfshield_local_incident_uv_flux(void);
#endif

#ifdef GALSF_FB_MECHANICAL
void determine_where_SNe_occur(void);
void mechanical_fb_calc_toplevel(void);
void verify_and_assign_local_mechfb_integrals(void);
#endif

#ifdef GALSF_FB_THERMAL
void determine_where_addthermalFB_events_occur(void);
void thermal_fb_calc(void);
#endif

#ifdef COOL_METAL_LINES_BY_SPECIES
/*double GetMetalLambda(double, double);*/
#ifndef CHIMES
double GetCoolingRateWSpecies(double nHcgs, double logT, double *Z);
double GetLambdaSpecies(long k_index, long index_x0y0, long index_x0y1, long index_x1y0, long index_x1y1, double dx, double dy, double dz, double mdz);
void LoadMultiSpeciesTables(void);
void ReadMultiSpeciesTables(int iT);
char *GetMultiSpeciesFilename(int i, int hk);
#endif
#endif

double sink_fb_angleweight(double sink_lum_input, Vec3<double>& sink_angle, double dx, double dy, double dz);
double sink_fb_angleweight_localcoupling(int j, double cos_theta, double r, double H_kernel);

#if defined(GALSF_SUBGRID_WINDS)
void assign_wind_kick_from_sf_routine(int i, double sm, double dtime, double* pvtau_return);
#endif

#if defined(GALSF_FB_FIRE_RT_LOCALRP)
void radiation_pressure_winds_consolidated(void);
#endif

#if defined(SINK_PARTICLES)
#endif

#ifdef GALSF_SUBGRID_WINDS
#if (GALSF_SUBGRID_WIND_SCALING==2)
void disp_setup_smoothinglengths(void);
void disp_density(void);
#endif
#endif


#ifdef COSMIC_RAY_SUBGRID_LEBRON
double cr_get_source_injection_rate(int i, struct particle_data *pp, struct gas_cell_data *cell);
double cr_get_source_shieldfac(int i, struct particle_data *pp, struct gas_cell_data *cell);
#endif


#ifdef CHIMES
double chimes_convert_u_to_temp(double u, double rho, int target, struct particle_data *pp, struct gas_cell_data *cell);
void chimes_update_gas_vars(int target, struct particle_data *pp, struct gas_cell_data *cell);
void chimes_gizmo_exit(void);
#ifdef COOL_METAL_LINES_BY_SPECIES
void chimes_update_element_abundances(int i, struct particle_data *pp, struct gas_cell_data *cell);
#endif
#ifdef CHIMES_TURB_DIFF_IONS
void chimes_update_turbulent_abundances(int i, int mode, struct particle_data *pp, struct gas_cell_data *cell);
#endif
#ifdef CHIMES_METAL_DEPLETION
void chimes_init_depletion_data(void);
void chimes_compute_depletions(double nH, double T, int thread_id);
#endif
#else
#endif
void cooling_parent_routine(void);
#ifdef NUCLEAR_NETWORK
void nuclear_parent_routine(void);
void InitNuclearNetwork(void);
#endif
void density(void);
void do_box_wrapping(void);
void energy_statistics(void);
void ensure_neighbours(void);

void output_log_messages(void);
void ewald_corr(double dx, double dy, double dz, double *fper);
void ewald_force(int ii, int jj, int kk, double x[3], double force[3]);
void ewald_init(void);
double ewald_psi(double x[3]);
double ewald_pot_corr(double dx, double dy, double dz);
integertime find_next_outputtime(integertime time);
integertime get_timestep(int p, double *a, int flag);
GIZMO_GPU_FUNCTION double return_timestep_dilation_factor(int i, int mode, struct particle_data *pp = P);
GIZMO_GPU_FUNCTION double timestep_dilation_factor(int i, int mode, struct particle_data *pp = P);
void refresh_timestep_dilation_factors_for_gpu(void);
GIZMO_GPU_FUNCTION double unit_integertime_in_physical(int i, struct particle_data *pp = P);
GIZMO_GPU_FUNCTION double get_physical_timestep_from_timebin(int bin, int i, struct particle_data *pp = P);
GIZMO_GPU_FUNCTION double get_particle_timestep_in_physical(int i, struct particle_data *pp = P);
GIZMO_GPU_FUNCTION double get_particle_feedback_timestep_in_physical(int i, struct particle_data *pp = P);

void gravity_tree(void);
void hydro_force(void);
void init(void);
GIZMO_GPU_FUNCTION void do_the_cooling_for_particle(int i, struct particle_data *pp, struct gas_cell_data *cell);
double GetCoolingTime(double u_old, double rho, double ne_guess, double *ne_eval, int target, struct particle_data *pp, struct gas_cell_data *cell);
double get_equilibrium_dust_temperature_estimate(int i, double shielding_factor_for_exgalbg, double T, struct particle_data *pp, struct gas_cell_data *cell);
double gas_dust_heating_coeff(int i, double T, double Tdust, struct particle_data *pp, struct gas_cell_data *cell);
double rt_eqm_dust_temp(int i, double T, double dust_absorption_rate, struct particle_data *pp, struct gas_cell_data *cell);
double dust_dEdt(int i, double T, double Tdust, double dust_absorption_rate, double fdustmet_init, struct particle_data *pp, struct gas_cell_data *cell);
double return_electron_fraction_from_heavy_ions(int target, double temperature, double density_cgs, double n_elec_HHe, struct particle_data *pp, struct gas_cell_data *cell);
MyFloat return_electron_fraction_from_Cplus(int target, MyFloat temp, MyFloat x_elec, MyFloat shieldfac, struct particle_data *pp, struct gas_cell_data *cell);
MyFloat return_electron_fraction_from_Oplus(int target, MyFloat nHp, struct particle_data *pp, struct gas_cell_data *cell);
MyFloat return_electron_fraction_from_molecular_ions(int target, MyFloat temp, struct particle_data *pp, struct gas_cell_data *cell);
MyFloat return_electron_fraction_from_alkali(int i, MyFloat temp, struct particle_data *pp, struct gas_cell_data *cell);
MyFloat get_FUV_G0(int i, MyFloat shieldfac, int mode, struct particle_data *pp, struct gas_cell_data *cell);
MyFloat f_Cplus(int i, MyFloat temp, MyFloat x_elec, MyFloat shieldfac, struct particle_data *pp, struct gas_cell_data *cell); 
MyFloat f_Oplus(MyFloat nHp);
MyFloat f_CO(int i, MyFloat temp, MyFloat x_elec, MyFloat shieldfac, MyFloat nHp, struct particle_data *pp, struct gas_cell_data *cell);
MyFloat alpha_recomb_grain(int i, MyFloat temp, MyFloat x_slec, MyFloat shieldfac, const char *ion_name, struct particle_data *pp, struct gas_cell_data *cell);
MyFloat grain_charge_psi(int i, MyFloat temp, MyFloat x_elec, MyFloat shieldfac, struct particle_data *pp, struct gas_cell_data *cell);
MyFloat total_ionization_rate_C(int i, MyFloat shieldfac, struct particle_data *pp, struct gas_cell_data *cell);
MyFloat cosmic_ray_ionization_rate_C(int i, struct particle_data *pp, struct gas_cell_data *cell);
MyFloat photoionization_rate_C(int i, MyFloat shieldfac, struct particle_data *pp, struct gas_cell_data *cell);
int ion_name_to_index(const char *ion_name, struct gas_cell_data *cell);
double get_starformation_rate(int i, int mode);
void update_internalenergy_for_galsf_effective_eos(int i, double tcool, double tsfr, double cloudmass_fraction, double rateOfSF);
void init_clouds(void);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
void mpi_printf(const char *fmt, ...);
void open_outputfiles(void);
void peano_hilbert_order(void);
void predict(double time);
void read_ic(char *fname);
void gizmo_register_hdf5_deflate_filter(void);  /* file_io/hdf5_deflate_filter.cc */
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
double my_second(void);
void set_softenings(void);
void set_units(void);
void setup_smoothinglengths(void);
void apply_special_boundary_conditions(int i, double mass_for_dp, int mode);

void minimum_large_ints(int n, long long *src, long long *res);
void sumup_large_ints(int n, int *src, long long *res);
void sumup_longs(int n, long long *src, long long *res);

void statistics(void);
double timediff(double t0, double t1);
void veldisp(void);
int binarySearch(const double * arr, const double x, const int l, const int r, const int total);

double get_gravkick_factor(integertime time0, integertime time1, int i, int mode);
double drift_integ(double a, void *param);
double gravkick_integ(double a, void *param);
double growthfactor_integ(double a, void *param);
void init_drift_table(void);
double get_drift_factor(integertime time0, integertime time1, int i, int mode);
double get_drift_factor_omp_safe(integertime time0, integertime time1, int i, int mode);
double measure_time(void);
double report_time(void);

/* on some DEC Alphas, the correct prototype for pow() is missing,
   even when math.h is included ! */
#ifndef __NVCOMPILER
double pow(double, double);
#endif


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

void pm_setup_nonperiodic_kernel(void);


#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE)
#ifdef CHIMES_STELLAR_FLUXES
double chimes_G0_luminosity(double stellar_age, double stellar_mass);
double chimes_ion_luminosity(double stellar_age, double stellar_mass);
int rt_get_source_luminosity_chimes(int i, int mode, double *lum, double *chimes_lum_G0, double *chimes_lum_ion, struct particle_data *pp, struct gas_cell_data *cell);
#endif
int rt_get_source_luminosity(int i, int mode, double *lum, struct particle_data *pp, struct gas_cell_data *cell);
int rt_get_donation_target_bin(int bin);
int rt_get_lum_band_stellarpopulation(int i, int mode, double *lum, struct particle_data *pp, struct gas_cell_data *cell);
int rt_get_lum_band_agn(int i, int mode, double *lum, struct particle_data *pp, struct gas_cell_data *cell);
int rt_get_lum_band_singlestar(int i, int mode, double *lum, struct particle_data *pp, struct gas_cell_data *cell);
void rt_define_effective_frequencies_in_bands(void);
double rt_kappa(int j, int k_freq, struct particle_data *pp, struct gas_cell_data *cell);
int check_if_absorbed_photons_can_be_reemitted_into_same_band(int kfreq);
double rt_absorb_frac_albedo(int j, int k_freq, struct particle_data *pp, struct gas_cell_data *cell);
double rt_absorption_rate(int i, int k_freq, struct particle_data *pp, struct gas_cell_data *cell);
#ifdef RT_SOLVER_EXPLICIT
GIZMO_GPU_FUNCTION inline double rt_diffusion_coefficient(int i, int k_freq, struct gas_cell_data *cell) {
    return cell[i].flux_limiter(k_freq) * C_LIGHT_CODE_REDUCED / (1.e-45 + cell[i].Rad_Kappa[k_freq] * cell[i].Density*All.cf_a3inv);
}
#endif
void rt_eddington_update_calculation(int j, struct gas_cell_data *cell);
void rt_update_driftkick(int i, double dt_entr, int mode, struct particle_data *pp, struct gas_cell_data *cell);
#endif
#ifdef RT_SOURCE_INJECTION
void rt_source_injection(void);
#endif
MyFloat dust_planck_mean_opacity(MyFloat Trad, MyFloat Tdust);

#ifdef TRANSPORT_SUBCYCLE
void transport_subcycle_exchange_fluxes(void);
void transport_subcycle_kick(void);
#endif

#ifdef RADTRANSFER
void rt_set_simple_inits(int RestartFlag);
#if defined(RT_EVOLVE_INTENSITIES)
void rt_init_intensity_directions(void);
#endif
void rt_get_lum_gas(int target, double *je, struct particle_data *pp, struct gas_cell_data *cell);
#ifdef RT_ISRF_BACKGROUND
void rt_apply_boundary_conditions(int i, struct particle_data *pp, struct gas_cell_data *cell);
void get_background_isrf_urad(int i, double *urad);
double background_isrf_cmb_Teff(void);
#endif
double slab_averaging_function(double x);
double blackbody_lum_frac(double E_lower, double E_upper, double T_eff);
double stellar_lum_in_band(int i, double E_lower, double E_upper, struct particle_data *pp, struct gas_cell_data *cell);
double rt_irband_egydensity_in_band(int i, double E_lower, double E_upper, struct gas_cell_data *cell);

#ifdef RT_CHEM_PHOTOION
/* rt_return_photon_number_density is now cell[i].rt_photon_number_density(k) */
double rt_photoion_chem_return_temperature(int i, double internal_energy);
void rt_update_chemistry(void);
void rt_get_sigma(void);
void rt_write_chemistry_stats(void);
#endif

#endif
double rt_kappa_adaptive_IR_band(int i, double T_dust, double Trad, int do_emission_absorption_scattering_opacity, int dust_or_gas_opacity_only_flag, struct particle_data *pp, struct gas_cell_data *cell);


void find_block(char *label,FILE *fd);


#if defined(COMPUTE_TIDAL_TENSOR_IN_GRAVTREE)
#ifdef PMGRID
void pmtidaltensor_periodic_diff(void);
int pmtidaltensor_nonperiodic_diff(int grnr);
void pmtidaltensor_periodic_fourier(int component);
int pmtidaltensor_nonperiodic_fourier(int component, int grnr);
#endif
#endif

int ags_gravity_kernel_shared_BITFLAG(short int particle_type_primary);
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
void ags_setup_smoothinglengths(void);
void ags_density(void);
int ags_density_isactive(int i);
double ags_return_maxsoft(int i);
double ags_return_minsoft(int i);
void AGSForce_calc(void);
#endif

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
void advect_mesh_point(int i, double dt);
double calculate_face_area_for_cartesian_mesh(const Vec3<double>& dp, double rinv, double l_side, Vec3<double>& Face_Area_Vec);
#endif

#if (SINGLE_STAR_TIMESTEPPING > 0)
void subtract_companion_gravity(int i);
#endif


void hydro_gradient_calc(void);
void construct_gradient(Vec3<MyDouble>& grad, int i);
#ifdef MHD_MODIFIED_GRADIENT
void mg_gradient_correction_calc(void);
#endif
inline void construct_gradient(double *grad, int i) { Vec3<MyDouble> v{grad[0],grad[1],grad[2]}; construct_gradient(v,i); grad[0]=v[0]; grad[1]=v[1]; grad[2]=v[2]; }
void local_slopelimiter(double *grad, double valmax, double valmin, double alim, double h, double shoot_tol, int pos_preserve, double d_max, double val_cen);
inline void local_slopelimiter(Vec3<double>& grad, double valmax, double valmin, double alim, double h, double shoot_tol, int pos_preserve, double d_max, double val_cen) { local_slopelimiter(grad.data, valmax, valmin, alim, h, shoot_tol, pos_preserve, d_max, val_cen); }

#ifdef TURB_DIFF_DYNAMIC
void dynamic_diff_vel_calc(void);
void dynamic_diff_calc(void);
#endif

#ifdef PARTICLE_EXCISION
void apply_excision();
#endif

#ifdef DM_SIDM
/* prob_of_interaction / g_geo / calculate_interact_kick moved to
   sidm/sidm_helper_functions.h (KOKKOS_INLINE_FUNCTION). Only the
   host-only tabulation and initialization helpers remain here. */
void init_geofactor_table(void);
double geofactor_integ(double x, void * params);
double geofactor_angle_integ(double u, void * params);
void init_self_interactions();
#ifdef GRAIN_COLLISIONS
/* return_grain_cross_section_per_unit_mass / prob_of_grain_interaction moved
   to solids/grain_helper_functions.h. */
#endif
#endif

#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES >= 4)
void special_rt_feedback_injection(void);
#endif

#ifdef CBE_INTEGRATOR
void do_cbe_initialization(void);
void do_cbe_drift_kick(int i, double dt);
/* do_cbe_flux_computation moved to sidm/cbe_integrator_functions.h
   (KOKKOS_INLINE_FUNCTION). */
void do_postgravity_cbe_calcs(int i);
#endif

#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
double do_cbe_nvt_inversion_for_faces(int i);
#endif

#ifdef DM_FUZZY
void do_dm_fuzzy_initialization(void);
void do_dm_fuzzy_drift_kick(int pindex, double dt_entr, int mode);
void DMGrad_gradient_calc(void);
/* do_dm_fuzzy_flux_computation / _old / dm_fuzzy_reconstruct_and_slopelimit
   moved to sidm/dm_fuzzy_functions.h (KOKKOS_INLINE_FUNCTION). */
#endif


#ifdef SINGLE_STAR_TIMESTEPPING
void kepler_timestep(int i, double dt, Vec3<double>& kick_dv, Vec3<double>& drift_dx, int mode);
void odeint_super_timestep(int i, double dt_super, Vec3<double>& kick_dv, Vec3<double>& drift_dx, int mode);
double gravfac(double r, double mass);
double gravfac2(double r, double mass);
void grav_accel_jerk(double mass, Vec3<double>& dx, Vec3<double>& dv, Vec3<double>& accel, Vec3<double>& jerk);
double eccentric_anomaly(double mean_anomaly, double ecc);
#endif

/* Kokkos/GPU lifecycle and sync functions (defined in cooling/cooling.cc).
   Declared unconditionally — the call sites in main.cc/begrun.cc/run.cc are
   retired (Step 5 C4); these sync functions are always compiled
   offload is active.  Declarations are harmless without definitions. */
void gizmo_kokkos_initialize(int argc, char *argv[]);
void gizmo_kokkos_finalize(void);
void gizmo_gpu_sync_all(void);
