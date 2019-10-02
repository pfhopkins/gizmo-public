#include "../../allvars.h"
/*
* This file was originally part of the GADGET3 code developed by Volker Springel.
* It has been updated significantly by PFH for basic compatibility with GIZMO,
* as well as code cleanups, and accommodating new GIZMO functionality for various
* other operations. See notes in subfind.c and GIZMO User Guide for details.
*/

#ifdef SUBFIND

#include "subfind.h"
#include "../fof.h"

int Ncollective;
int MaxNsubgroups;
int Nsubgroups;
int TotNsubgroups;
struct subgroup_properties *SubGroup;
struct nearest_r2_data *R2Loc;
struct nearest_ngb_data *NgbLoc;
struct r2data *R2list;
//struct nearest_ngb_data *NgbLoc;
//struct nearest_r2_data *R2Loc;
double *Dist2list;

#endif
