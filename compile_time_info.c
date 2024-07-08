#include <stdio.h>
void output_compile_time_options(void)
{
printf(
"        BOX_SPATIAL_DIMENSION=1\n"
"        BOX_PERIODIC\n"
"        HYDRO_MESHLESS_FINITE_MASS\n"
"        HYDRO_PRESSURE_SPH\n"
"        EOS_GAMMA=(5.0/3.0)\n"
"        SELFGRAVITY_OFF\n"
"        OUTPUT_IN_DOUBLEPRECISION\n"
"        INPUT_IN_DOUBLEPRECISION\n"
"\n");
}
