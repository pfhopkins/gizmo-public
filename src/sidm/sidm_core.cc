#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "../declarations/gizmo_quadrature.h"

#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"

#define GSLWORKSIZE 100000

/*! \file sidm_core.cc
 *  \brief Host-only startup initialization for the DM_SIDM module:
 *         per-particle
 *         scatter-counter zeroing. Called once each from core/begrun.cc
 *         and core/init.cc.
 *
 *         The per-pair scatter physics (prob_of_interaction, isotropic
 *         kick) lives in sidm/sidm_helper_functions.h as
 *         KOKKOS_INLINE_FUNCTION and is invoked from the AgsForceSpec
 *         pair kernel via sidm/sidm_core_flux_functions.h. Nothing in
 *         the per-step hot path is host-only any more.
 *
 *         Originally written by Miguel Rocha, rocham@uci.edu (Oct 2010);
 *         updated 2014 and re-written by PFH March 2018; GPU migration
 *         of the hot path completed with the AgsForceSpec port.
 */

#ifdef DM_SIDM

/*! This routine initializes the table that will be used to get the geometrical factor
 *  as a function of the two particle separations. It populates a table with the results of the numerical integration */



/*! This function simply initializes some variables to prevent memory errors */
void init_self_interactions() {int i; for(i = 0; i < NumPart; i++) {P[i].dtime_sidm = 0; P[i].NInteractions = 0;}}

#endif
