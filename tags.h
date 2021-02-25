/* defines tag names for labeling MPI messages to prevent collisions */

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * in part by Phil Hopkins (phopkins@caltech.edu) for GIZMO (new variables,
 * and different naming conventions for some old variables)
 */


#define TAG_N             10      /*!< Various tags used for labeling MPI messages */
#define TAG_HEADER        11
#define TAG_PDATA         12
#define TAG_SPHDATA       13
#define TAG_KEY           14
#define TAG_DMOM          15
#define TAG_NODELEN       16
#define TAG_HMAX          17
#define TAG_GRAV_A        18
#define TAG_GRAV_B        19
#define TAG_DIRECT_A      20
#define TAG_DIRECT_B      21
#define TAG_HYDRO_A       22 
#define TAG_HYDRO_B       23
#define TAG_NFORTHISTASK  24
#define TAG_PERIODIC_A    25
#define TAG_PERIODIC_B    26
#define TAG_PERIODIC_C    27
#define TAG_PERIODIC_D    28
#define TAG_NONPERIOD_A   29 
#define TAG_NONPERIOD_B   30
#define TAG_NONPERIOD_C   31
#define TAG_NONPERIOD_D   32
#define TAG_POTENTIAL_A   33
#define TAG_POTENTIAL_B   34
#define TAG_DENS_A        35
#define TAG_DENS_B        36
#define TAG_LOCALN        37

#define TAG_BH_A          38
#define TAG_BH_B          39
#define TAG_BH_C          40
#define TAG_BH_D          41
#define TAG_BH_E          42
#define TAG_BH_F          43
#define TAG_BH_G          44
#define TAG_BH_H          45
#define TAG_BH_I          46

#define TAG_FOF_A         47
#define TAG_FOF_B         48
#define TAG_FOF_C         49
#define TAG_FOF_D         50
#define TAG_FOF_E         51
#define TAG_FOF_F         52
#define TAG_FOF_G         53
#define TAG_FOF_H         54
#define TAG_FOF_I         55
#define TAG_FOF_J         56

#define TAG_SWAP          57
#define TAG_PM_FOLD       58

#define TAG_AGS_DENS_A    59
#define TAG_AGS_DENS_B    60

#define TAG_GRDENS_A      61
#define TAG_GRDENS_B      62

#define TAG_DMDENS_A      63
#define TAG_DMDENS_B      64

#define TAG_PSRT_A        65
#define TAG_PDATA_SPH     66
#define TAG_KEY_SPH       67

#define TAG_RT_A          72
#define TAG_RT_B          73
#define TAG_RT_C          74
#define TAG_RT_D          75

#define TAG_FBLOOP_1A     80
#define TAG_FBLOOP_1B     81
#define TAG_FBLOOP_2A     82
#define TAG_FBLOOP_2B     83
#define TAG_FBLOOP_3A     84
#define TAG_FBLOOP_3B     85
#define TAG_FBLOOP_4A     86
#define TAG_FBLOOP_4B     87
#define TAG_FBLOOP_5A     88
#define TAG_FBLOOP_5B     89

#define TAG_GRADLOOP_A    90
#define TAG_GRADLOOP_B    91
#define TAG_GRADLOOP_C    92

#define TAG_PDATA_STARS   93
#define TAG_KEY_STARS     94
#define TAG_PDATA_BHS     95
#define TAG_KEY_BHS       96
#define TAG_BHDATA        97

#define TAG_DYNSMAGLOOP_A 98
#define TAG_DYNSMAGLOOP_B 99
#define TAG_DYNSMAGLOOP_C 100

#ifdef CHIMES
#define TAG_CHIMESDATA 101
#define TAG_ABUNDATA 102
#endif

#define TAG_MPI_GENERIC_COM_BUFFER_A 103
#define TAG_MPI_GENERIC_COM_BUFFER_B 104

