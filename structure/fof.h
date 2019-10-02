#ifndef FOF_H
#define FOF_H
#include "../allvars.h"

/* variables for fof.c code */
/*
* This file was originally part of the GADGET3 code developed by Volker Springel.
* It has been updated significantly by PFH for basic compatibility with GIZMO,
* as well as code cleanups, and accommodating new GIZMO functionality for various
* other operations.
*/

void fof_exchange_group_data(void);
int fof_compare_FOF_PList_MinID(const void *a, const void *b);
int fof_compare_FOF_GList_MinID(const void *a, const void *b);
int fof_compare_FOF_GList_MinIDTask(const void *a, const void *b);
int fof_compare_FOF_GList_MinIDTask_MinID(const void *a, const void *b);
int fof_compare_FOF_GList_LocCountTaskDiffMinID(const void *a, const void *b);
int fof_compare_FOF_GList_ExtCountMinID(const void *a, const void *b);
int fof_compare_Group_GrNr(const void *a, const void *b);
int fof_compare_Group_MinIDTask(const void *a, const void *b);
int fof_compare_Group_MinID(const void *a, const void *b);
int fof_compare_ID_list_GrNrID(const void *a, const void *b);
int fof_compare_Group_MinIDTask_MinID(const void *a, const void *b);
int fof_compare_Group_Len(const void *a, const void *b);
void fof_finish_group_properties(void);
int compare_group_mass_ID(const void *a, const void *b);
void fof_assign_HostHaloMass(void);
extern int Ngroups, TotNgroups;
extern long long TotNids;




typedef struct
{
  MyIDType MinID;
  MyIDType MinIDTask;
  int LocCount;
  int ExtCount;
#ifdef FOF_DENSITY_SPLIT_TYPES
  int LocDMCount;
  int ExtDMCount;
#endif
  int GrNr;
} fof_group_list; 

typedef struct
{
  int Len;
  unsigned int Offset;
  MyIDType MinID;
  MyIDType MinIDTask;
  int GrNr;
  int LenType[6];
  MyOutputFloat MassType[6];
  MyOutputFloat Mass;
  MyOutputFloat CM[3];
  MyOutputFloat Vel[3];
  MyDouble FirstPos[3];
#ifdef GALSF
  double Sfr;
#endif
#ifdef BLACK_HOLES
  MyOutputFloat BH_Mass;
  MyOutputFloat BH_Mdot;
#ifdef BH_SEED_FROM_FOF
  MyOutputFloat MinPot;
  int index_maxdens, task_maxdens;
#endif
#endif
    
#ifdef SUBFIND
  int Nsubs;
  int FirstSub;
  MyDouble Pos[3];
  int ContaminationLen;
  MyOutputFloat ContaminationMass;
  struct Subfind_DensityOtherPropsEval_data_out SubHaloProps_vsDelta[SUBFIND_ADDIO_NUMOVERDEN];
#endif

} group_properties;
extern group_properties *Group;

typedef struct
{
  MyIDType ID;
  unsigned int GrNr;
} fof_id_list;

#endif // FOF_H
