/*! \file blackhole.h
 *  \brief routine declarations for gas accretion onto black holes, and black hole mergers
 */
/*
 * This file is largely written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 *   It was based on a similar file in GADGET3 by Volker Springel (volker.springel@h-its.org),
 *   but the physical modules for black hole accretion and feedback have been
 *   replaced, and the algorithm for their coupling is new to GIZMO.  This file was modified
 *   on 1/9/15 by Paul Torrey (ptorrey@mit.edu) for clarity by parsing the existing code into
 *   smaller files and routines.  Some communication and black hole structures were modified
 *   to reduce memory usage. Cleanup, de-bugging, and consolidation of routines by Xiangcheng Ma
 *   (xchma@caltech.edu) followed on 05/15/15; re-integrated by PFH.
 */

#ifndef gizmo_blackhole_h
#define gizmo_blackhole_h


#ifndef BH_CSND_FRAC_BH_MERGE
/* Relative velocity fraction (in units of soundspeed) for merging blackholes, default=1.0 */
#define BH_CSND_FRAC_BH_MERGE 1.0
#endif



#if defined(BLACK_HOLES)
/* blackhole_utils.c */
void blackhole_start(void);
void blackhole_end(void);

/* blackholes.c */
void blackhole_properties_loop(void);
void blackhole_final_operations(void);

/* blackhole_environment.c */
void blackhole_environment_loop(void);
int blackhole_environment_evaluate(int target, int mode, int *nexport, int *nSend_local);
#ifdef BH_GRAVACCRETION
void blackhole_environment_second_loop(void);
int blackhole_environment_second_evaluate(int target, int mode, int *nexport, int *nSend_local);
#endif

/* blackhole_swallow_and_kick.c */
void blackhole_swallow_and_kick_loop(void);
int blackhole_swallow_and_kick_evaluate(int target, int mode, int *nexport, int *nSend_local);

/* blackhole_feed.c */
void blackhole_feed_loop(void);
int blackhole_feed_evaluate(int target, int mode, int *nexport, int *nSend_local);


void out2particle_blackhole(struct blackhole_temp_particle_data *out, int target, int mode);

//void check_for_bh_merger(int j, MyIDType id);
double bh_eddington_mdot(double bh_mass);
double bh_lum_bol(double mdot, double mass, long id);
int bh_check_boundedness(int j, double vrel, double vesc, double dr_code, double sink_radius);
double bh_vesc(int j, double mass, double r_code, double bh_softening);
void normalize_temp_info_struct(int i);
void set_blackhole_mdot(int i, int n, double dt);
void set_blackhole_new_mass(int i, int n, double dt);
#if defined(BH_DRAG) || defined(BH_DYNFRICTION)
void set_blackhole_drag(int i, int n, double dt);
#endif
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)
void set_blackhole_long_range_rp(int i, int n);
#endif

#endif


#endif
