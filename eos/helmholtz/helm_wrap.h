/* wrappers for the Timmes EOS Fortran code. 
    written by David Radice
    WARNING: because of the memory layout used by the original Fortran code 
    which was not thread safe, the Helmholtz wrapper might perform poorly or 
    throw an exception if a large number of OpenMP threads is used. 
    Eventually, someone will have to make a thread safe version of the 
    Helmholtz EOS if you want to run on things like the new KNL on Stampede2
 */

#ifndef EOS_HELM_WRAP_H
#define EOS_HELM_WRAP_H

void helm_read_table_c(
        char const * tab_file_name);

void helm_range_rho_ye_c(
        double * rho_ye_min,
        double * rho_ye_max);
void helm_range_temp_c(
        double * temp_min,
        double * temp_max);

void helm_range_eps_c(
        int const * rank,
        double const * rho,
        double const * abar,
        double const * ye,
        double * eps_min,
        double * eps_max,
        int * eosfail);

void helm_eos_t_c(
        int const * rank,
        double const * rho,
        double const * temp,
        double const * abar,
        double const * ye,
        double * press,
        double * eps,
        double * entropy,
        double * csound,
        double * dedt,
        int * eosfail);

void helm_eos_e_c(
        int const * rank,
        double const * rho,
        double const * eps,
        double const * abar,
        double const * ye,
        double * temp,
        double * press,
        double * entropy,
        double * csound,
        double * dedt,
        int * eosfail);

#endif
