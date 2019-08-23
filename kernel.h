/* this file contains the definitions for the kernel functions */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * substantially by Phil Hopkins (phopkins@caltech.edu) for GIZMO (added the 
 * calls to this for other uses of the kernel in-code in areas like gravity;
 * also added new kernels here)
 */

#if !defined(KERNEL_FUNCTION)
#define KERNEL_FUNCTION 3  // default to cubic spline kernel
#endif


#if (KERNEL_FUNCTION == 1) // linear ramp (not a good kernel for many numerical reasons; here for testing purposes)
#define KERNEL_CORE_SIZE (2.0/3.0)
#if (NUMDIMS==1)
#define  KERNEL_NORM  (1.0)        	         /*!< For 1D-normalized kernel */
#define  KERNEL_NMIN  2 // minimum number of neighbors for this kernel and dimension
#define  KERNEL_NMAX  8 // maximum number of neighbors for this kernel and dimension
#elif (NUMDIMS==2)
#define  KERNEL_NORM  (3.0/M_PI)	         /*!< For 2D-normalized kernel */
#define  KERNEL_NMIN  12
#define  KERNEL_NMAX  24
#else
#define  KERNEL_NORM  (3.0/M_PI)                  /*!< For 3D-normalized kernel */
#define  KERNEL_NMIN  30
#define  KERNEL_NMAX  64
#endif
#endif



#if (KERNEL_FUNCTION == 2) // quadratic 'peak' kernel
#define KERNEL_CORE_SIZE (1.0/2.0)
#if (NUMDIMS==1)
#define  KERNEL_NORM  (3.0/2.0)        	         /*!< For 1D-normalized kernel */
#define  KERNEL_NMIN  2 // minimum number of neighbors for this kernel and dimension
#define  KERNEL_NMAX  8 // maximum number of neighbors for this kernel and dimension
#elif (NUMDIMS==2)
#define  KERNEL_NORM  (6.0/M_PI)	         /*!< For 2D-normalized kernel */
#define  KERNEL_NMIN  12
#define  KERNEL_NMAX  24
#else
#define  KERNEL_NORM  (15.0/(2.0*M_PI))                  /*!< For 3D-normalized kernel */
#define  KERNEL_NMIN  30
#define  KERNEL_NMAX  64
#endif
#endif


#if (KERNEL_FUNCTION == 3) // Morris 1996 cubic spline kernel
#define KERNEL_CORE_SIZE (1.0/2.0)
#if (NUMDIMS==1)
#define  KERNEL_NORM  (4.0/3.0)        	         /*!< For 1D-normalized kernel */
#define  KERNEL_NMIN  2 // minimum number of neighbors for this kernel and dimension
#define  KERNEL_NMAX  8 // maximum number of neighbors for this kernel and dimension
#elif (NUMDIMS==2)
#define  KERNEL_NORM  (40.0/(7.0*M_PI))	         /*!< For 2D-normalized kernel */
#define  KERNEL_NMIN  12
#define  KERNEL_NMAX  24
#else
#define  KERNEL_NORM  (8.0/M_PI)                  /*!< For 3D-normalized kernel */
#define  KERNEL_NMIN  30
#define  KERNEL_NMAX  64
#endif
#endif


#if (KERNEL_FUNCTION == 4) // Morris 1996 quartic spline kernel
#define KERNEL_CORE_SIZE (2.0/5.0)
#if (NUMDIMS==1)
#define  KERNEL_NORM  (243.0/40.0)       	    /*!< For 1D-normalized kernel */
#define  KERNEL_NMIN  3
#define  KERNEL_NMAX  8
#elif (NUMDIMS==2)
#define  KERNEL_NORM  (15309.0/(478.0*M_PI))	    /*!< For 2D-normalized kernel */
#define  KERNEL_NMIN  16
#define  KERNEL_NMAX  32
#else
#define  KERNEL_NORM  (2187.0/(40.0*M_PI))	    /*!< For 3D-normalized kernel */
#define  KERNEL_NMIN  40
#define  KERNEL_NMAX  125
#endif
#endif


#if (KERNEL_FUNCTION == 5) // Morris 1996 quintic spline kernel
#define KERNEL_CORE_SIZE (1.0/3.0)
#if (NUMDIMS==1)
#define  KERNEL_NORM  (3125.0/768.0)        	    /*!< For 1D-normalized kernel */
#define  KERNEL_NMIN  4
#define  KERNEL_NMAX  12
#elif (NUMDIMS==2)
#define  KERNEL_NORM  (46875.0/(2398.0*M_PI))	    /*!< For 2D-normalized kernel */
#define  KERNEL_NMIN  20
#define  KERNEL_NMAX  50
#else
#define  KERNEL_NORM  (15625.0/(512.0*M_PI))	    /*!< For 3D-normalized kernel */
#define  KERNEL_NMIN  60
#define  KERNEL_NMAX  180
#endif
#endif


#if (KERNEL_FUNCTION == 6) // Wendland C2 kernel
#define KERNEL_CORE_SIZE (2.0/5.0)
#if (NUMDIMS==1)
#define  KERNEL_NORM  (5.0/4.0)       	    /*!< For 1D-normalized kernel */
#define  KERNEL_NMIN  3
#define  KERNEL_NMAX  8
#elif (NUMDIMS==2)
#define  KERNEL_NORM  (7.0/M_PI)	    /*!< For 2D-normalized kernel */
#define  KERNEL_NMIN  16
#define  KERNEL_NMAX  32
#else
#define  KERNEL_NORM  (21.0/(2.0*M_PI))	    /*!< For 3D-normalized kernel */
#define  KERNEL_NMIN  40
#define  KERNEL_NMAX  125
#endif
#endif


#if (KERNEL_FUNCTION == 7) // Wendland C4 kernel
#define KERNEL_CORE_SIZE (1.0/3.0)
#if (NUMDIMS==1)
#define  KERNEL_NORM  (3.0/2.0)        	    /*!< For 1D-normalized kernel */
#define  KERNEL_NMIN  4
#define  KERNEL_NMAX  12
#elif (NUMDIMS==2)
#define  KERNEL_NORM  (9.0/M_PI)	    /*!< For 2D-normalized kernel */
#define  KERNEL_NMIN  20
#define  KERNEL_NMAX  50
#else
#define  KERNEL_NORM  (495.0/(32.0*M_PI))	    /*!< For 3D-normalized kernel */
#define  KERNEL_NMIN  60
#define  KERNEL_NMAX  180
#endif
#endif


#if (KERNEL_FUNCTION == 8) // quadratic 2-step kernel ('magic kernel' for image filtering is this in 2d with KERNEL_U0=1/3)
#define KERNEL_CORE_SIZE (1.0/2.0)
#if !defined(KERNEL_U0)
#define KERNEL_U0        (0.2) // set to anything between 0 (most 'peaked') to 1 ('flattest'): values ~0.55-0.59 minimize second derivative
#endif
#if (NUMDIMS==1)
#define  KERNEL_NORM  (3.0/(2.0 + KERNEL_U0)) /*!< For 1D-normalized kernel */
#define  KERNEL_NMIN  2 // minimum number of neighbors for this kernel and dimension
#define  KERNEL_NMAX  8 // maximum number of neighbors for this kernel and dimension
#elif (NUMDIMS==2)
#define  KERNEL_NORM  (6.0/(M_PI*(1 + KERNEL_U0 + KERNEL_U0*KERNEL_U0))) /*!< For 2D-normalized kernel */
#define  KERNEL_NMIN  12
#define  KERNEL_NMAX  24
#else
#define  KERNEL_NORM  (15.0/(2.0*M_PI*(1+KERNEL_U0)*(1+KERNEL_U0*KERNEL_U0))) /*!< For 3D-normalized kernel */
#define  KERNEL_NMIN  30
#define  KERNEL_NMAX  64
#endif
#endif





static inline void kernel_hinv(double h, double *hinv, double *hinv3, double *hinv4)
{
  *hinv = 1.0 / h;

#if (NUMDIMS==1)
    *hinv3 = *hinv;
#elif (NUMDIMS==2)
    *hinv3 = *hinv * *hinv;
#else // THREEDIMS
    *hinv3 = *hinv * *hinv * *hinv;
#endif

  *hinv4 = *hinv3 * *hinv;

  return;
} 

/* Attention: Here we assume that kernel is only called 
   with range 0..1 for u as done in hydra or density !! 
   Call with mode 0 to calculate dwk and wk
   Call with mode -1 to calculate only wk
   Call with mode +1 to calculate only dwk */

static inline void kernel_main(double u, double hinv3, double hinv4, 
        double *wk, double *dwk, int mode)
{
#if (KERNEL_FUNCTION == 1) /* linear ramp */
    if(mode >= 0)
        *dwk = -1;
    if(mode <= 0)
        *wk = 1-u;
#endif

#if (KERNEL_FUNCTION == 2) /* quadratic */
    double t1 = 1-u;
    if(mode >= 0)
        *dwk = -2*t1;
    if(mode <= 0)
        *wk = t1*t1;
#endif
    

#if (KERNEL_FUNCTION == 3) /* cubic spline */
  if(u < 0.5)
    {
      if(mode >= 0) 
          *dwk = u * (18.0 * u - 12.0);
      if(mode <= 0) 
          *wk = (1.0 + 6.0 * (u - 1.0) * u * u);
    }
  else
    {
      double t1 = (1.0 - u);
      double t2 = t1 * t1;
      if(mode >= 0) 
          *dwk = -6.0 * t2;
      if(mode <= 0) 
          *wk = 2.0 * t2 * t1;
    }
#endif /* cubic spline */

#if (KERNEL_FUNCTION == 4) /* quartic spline */
  double t1 = (1.0 - u);
  double t2 = t1 * t1;
  double t4 = t2 * t2;

  if(mode >= 0) 
      *dwk = -5.0 * t4;
  if(mode <= 0) 
      *wk = t4 * t1;

  if (u < 2.0/3.0)
    {
      t1 = (2.0/3.0 - u);
      t2 = t1 * t1;
      t4 = t2 * t2;
      if(mode >= 0) 
          *dwk += 30.0 * t4;
      if(mode <= 0) 
          *wk -= 6.0 * t4 * t1;
    }
  if (u < 1.0/3.0)
    {
      t1 = (1.0/3.0 - u);
      t2 = t1 * t1;
      t4 = t2 * t2;
      if(mode >= 0) 
          *dwk -= 75.0 * t4;
      if(mode <= 0) 
          *wk += 15.0 * t4 * t1;
    }
#endif /* quartic spline */

#if (KERNEL_FUNCTION == 5) /* quintic spline */
    double t1 = (1.0 - u);
    double t2 = t1 * t1;
    
    if(mode >= 0)
        *dwk = -4.0 * t2 * t1;
    if(mode <= 0)
        *wk = t2 * t2;
    
    if (u < 0.6)
    {
        t1 = (0.6 - u);
        t2 = t1 * t1;
        if(mode >= 0)
            *dwk += 20.0 * t2 * t1;
        if(mode <= 0)
            *wk -= 5.0 * t2 * t2;
    }
    if (u < 0.2)
    {
        t1 = (0.2 - u);
        t2 = t1 * t1;
        if(mode >= 0)
            *dwk -= 40.0 * t2 * t1;
        if(mode <= 0)
            *wk += 10.0 * t2 * t2;
    }
#endif /* quintic spline */
    

#if (KERNEL_FUNCTION == 6) /* Wendland C2 */
    double t1 = (1 - u);
    double t3 = t1*t1*t1;
#if (NUMDIMS == 1)
    if(mode >= 0)
        *dwk = -12.0 * u * t1*t1;
    if(mode <= 0)
        *wk = t3 * (1.0 + 3.0*u);
#else
    if(mode >= 0)
        *dwk = -20.0 * u * t3;
    if(mode <= 0)
        *wk = t3 * t1 * (1.0 + 4.0*u);
#endif
#endif

    
#if (KERNEL_FUNCTION == 7) /* Wendland C4 */
    double t1 = (1 - u);
    double t5 = t1*t1; t5 *= t5*t1;
#if (NUMDIMS == 1)
    if(mode >= 0)
        *dwk = -14.0 * (t5/t1) * u * (1.0 + 4.0*u);
    if(mode <= 0)
        *wk = t5 * (1.0 + 5.0*u + 8.0*u*u);
#else
    if(mode >= 0)
        *dwk = -(56.0/3.0) * t5 * u * (1.0 + 5.0*u);
    if(mode <= 0)
        *wk = t5 * t1 * (1.0 + 6.0*u + (35.0/3.0)*u*u);
#endif
#endif

    
#if (KERNEL_FUNCTION == 8) /* quadratic '2-part' kernel */
    if(u < KERNEL_U0)
    {
        if(mode >= 0)
            *dwk = -2*u/KERNEL_U0;
        if(mode <= 0)
            *wk = 1-u*u/KERNEL_U0;
    } else {
        if(mode >= 0)
            *dwk = -2*(1-u)/(1-KERNEL_U0);
        if(mode <= 0)
            *wk = (1-u)*(1-u)/(1-KERNEL_U0);
    }
#endif
    
    
  if(mode >= 0) 
      *dwk *= KERNEL_NORM * hinv4;
  if(mode <= 0) 
      *wk *= KERNEL_NORM * hinv3;

  return;
}


/* this defines the kernel for the short-range gravitational softening, 
 which does not have to correspond to that of the gas (although for the gas
 itself, it should). dwk is the force kernel, wk is the potential kernel
  Call with mode 0 to calculate only dphi_dh (for zeta correction)
  Call with mode -1 to calculate only phi
  Call with mode +1 to calculate only (1/u) * dphi_du */


static inline double kernel_gravity(double u, double hinv, double hinv3, int mode)
{
    /* here everything is newtonian, add this as a check just in case */
    if(u >= 1)
    {
        if(mode ==  0) return 0;
        if(mode ==  1) return hinv3/(u*u*u);
        if(mode == -1) return -hinv/u;
    }
    double wk;


#if (KERNEL_FUNCTION == 1) /* linear ramp */
    if(mode == 1)
    {
        wk = 4 - 3*u;
        return wk * hinv3;
    }
    else if(mode == -1)
    {
        wk = -2 + u*u*(2-u);
        return wk * hinv;
    }
    else
    {
        wk = 2*(1-u)*(1-u)*(1+2*u);
        return wk * hinv * hinv;
    }
#endif

    
#if (KERNEL_FUNCTION == 2) /* quadratic */
    if(mode == 1)
    {
        wk = 10. + 3.*u * (-5. + 2.*u);
        return wk * hinv3;
    }
    else if(mode == -1)
    {
        double uu = u*u;
        wk = -2.5 + 5.*uu*(1.-u) + 1.5*uu*uu;
        return wk * hinv;
    }
    else
    {
        double uu = 1 - u;
        wk = 2.5 * uu*uu*uu * (1. + 3.*u);
        return wk * hinv * hinv;
    }
#endif


#if (KERNEL_FUNCTION == 3) /* cubic spline */
    if(mode == 1)
    {
        if(u < 0.5)
            wk = (10.666666666667 + u * u * (32.0 * u - 38.4));
        else
            wk = (21.333333333333 - 48.0 * u +
                  38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
        return wk * hinv3;
    }
    else if(mode == -1)
    {
        if(u < 0.5)
            wk = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
        else
            wk = -3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
                                                      u * (-16.0 + u * (9.6 - 2.133333333333 * u)));
        return wk * hinv;
    }
    else
    {
        if(u < 0.5)
            wk = 2.8 + 16.0 * u * u * (-1.0 + 3.0 * u * u * (1.0 - 0.8 * u));
        else
            wk = 3.2 + 32.0 * u * u * (-1.0 + u * (2.0 - 1.5 * u + 0.4 * u * u));
        return wk * hinv * hinv;
    }
#endif /* cubic spline */

    
#if (KERNEL_FUNCTION == 4) /* quartic spline */
    double u2=u*u;
    if(mode == 1)
    {
        if(u < 0.2)
            wk = 125.*(161. - 630.*u2 + 1125.*u2*u2) / 1344.;
        else
            if(u < 0.6)
                wk = (1. - 625.*u2*u*(-154. + 5.*u*(-21. + 2.*u*(126. + 25.*u*(-7. + 3.*u))))) / (6720.*u2*u);
            else
                wk = (-437. + 3125.*u2*u*(35. + u*(-105. + u*(126. + 5.*u*(-14. + 3.*u))))) / (2688.*u2*u);
        return wk * hinv3;
    }
    else if(mode == -1)
    {
        if(u < 0.2)
            wk = (-8393. + 125.*u2 * (161. - 315.*u2 + 375.*u2*u2)) / 2688.;
        else
            if(u < 0.6)
                wk = -(1. + 5.*u*(4193. + 125.*u2*(-77. + 5.*u*(-7. + u*(63. + 5.*u*(-14. + 5.*u)))))) / (6720.*u);
            else
                wk = (874. + 3125.*u*(-7. + u2*(35. + u*(-70. + u*(63. + u*(-28. + 5.*u)))))) / (5376.*u);
        return wk * hinv;
    }
    else
    {
        if(u < 0.2)
            wk = (1199. - 375.*u2*(23. - 75.*u2 + 125.*u2*u2)) / 384.;
        else
            if(u < 0.6)
                wk = (599. + 125.*u2*(-33. + 5.*u*(-4. + 5.*u*(9. + u*(-12. + 5.*u))))) / 192.;
            else
            {
                double um=1.0-u; u2=um*um;
                wk = (3125./768.) * u2*u2*um * (1. + 5.*u);
            }
        return wk * hinv * hinv;
    }
#endif /* quartic spline */
    
    
    
#if (KERNEL_FUNCTION == 5) /* quintic spline */
    double u2=u*u;
    if(mode == 1)
    {
        if(u < 1./3.)
            wk = -(9./280.)*(-616.+27.*u2*(112.+45.*u2*(-8.+7.*u)));
        else
            if(u < 2./3.)
                wk = (5.+27.*u2*u*(952.+9.*u*(350.+3.*u*(-784.+5.*u*(280.+9.*u*(-24.+7.*u)))))) / (1680.*u2*u);
            else
                wk = -(169.+729.*u2*u*(-56.+u*(210.+u*(-336.+u*(280.+3.*u*(-40.+7.*u)))))) / (560.*u2*u);
        return wk * hinv3;
    }
    else if(mode == -1)
    {
        if(u < 1./3.)
            wk = (-956.-9.*u2*(-308.+27.*u2*(28.+15.*u2*(-4.+3.*u)))) / 280.;
        else
            if(u < 2./3.)
                wk = (-5.+3.*u*(-1892.+9.*u2*(476.+3.*u*(350.+9.*u*(-196.+5.*u*(56.+9.*u*(-4.+u))))))) / (1680.*u);
            else
                wk = (169.-729*u*(4.+(-2.+u)*u2*(14.+u*(-28.+u*(28.+u*(-14.+3.*u)))))) / (560.*u);
        return wk * hinv;
    }
    else
    {
        if(u < 1./3.)
            wk = (239.+27.*u2*(-77.+45.*u2*(7.+3.*u2*(-7.+6.*u)))) / 70.;
        else
            if(u < 2./3.)
                wk = (473.-27.*u2*(119.+5.*u*(70.+9*u*(-49.+3*u*(28.+3.*u*(-7.+2.*u)))))) / 140.;
            else
            {
                double um=1.0-u; um*=um*um; um*=um;
                wk = (729./140.) * um * (1.+6.*u);
            }
        return wk * hinv * hinv;
    }
#endif /* quintic spline */

    
#if (KERNEL_FUNCTION == 6) /* Wendland C2 */
    if(mode == 1)
    {
        wk = 14. + u*u * (-84. + u*(140. + 3.*u*(-30. + 7.*u)));
        return wk * hinv3;
    }
    else if(mode == -1)
    {
        double uu = u*u;
        wk = -3. + uu * (7. + uu * (-21. + u * (28. + 3. * (-5. + u) * u)));
        return wk * hinv;
    }
    else
    {
        double t1 = 1 - u;
        double t2 = t1*t1; t2 *= t2; t2 *= t1;
        wk = 3. * t2 * (1. + u * (5. + 8. * u));
        return wk * hinv * hinv;
    }
#endif

    
#if (KERNEL_FUNCTION == 7) /* Wendland C4 */
    if(mode == 1)
    {
        double uu = u*u;
        wk = 0.125 * (165 + uu*(-924 + uu*(4950 + u*(-9240 + u*(7700 - 3168*u + 525*uu)))));
        return wk * hinv3;
    }
    else if(mode == -1)
    {
        double uu = u*u;
        wk = 0.0625 * (-55 + uu*(165 + uu*(-462 + uu*(1650 + u*(-2640 + u*(1925 + u*(-704 + 105*u)))))));
        return wk * hinv;
    }
    else
    {
        double t1 = 1 - u;
        double t2 = t1*t1*t1; t2 *= t2; t2 *= t1;
        wk = (55./16.) * t2 * (1 + 3*u) * (1 + u*(4 + 7*u));
        return wk * hinv * hinv;
    }
#endif
    
    
    return 0;
}


