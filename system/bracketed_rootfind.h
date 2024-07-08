/*
Code block containing a Brent 1973 rootfind routine, guaranteed to find the root
provided upper and lower bounds for a continuous function

Required initializations outside this block:
    ROOTFIND_FUNCTION(x) - #define'd macro that evaluates a function of a single
variable whose root we wish to find
    ROOTFIND_X_a, ROOTFIND_X_b - Arguments to ROOTFIND_FUNCTION that are known
to bracket the root.
    ROOTFUNC_a, ROOTFUNC_b - values of the function evaluated at the bracket
points ROOTFIND_REL_X_tol - Tolerance for desired relative error in the root;
stop iterating when this is achieved ROOTFIND_ABS_X_tol - Tolerance for desired
absolute error in the root; stop iterating when this is achieved
*/

// if doing nested rootfinds, need to def a different inner function because
// ROOTFIND_FUNCTION has not yet been undef'd
#ifdef ROOTFIND_FUNCTION_INNER
#define ROOTFUNC ROOTFIND_FUNCTION_INNER
#else
#define ROOTFUNC ROOTFIND_FUNCTION
#endif

if (ROOTFUNC_a * ROOTFUNC_b > 0) {
    PRINT_WARNING("ERROR: Bounds supplied to bracketed_roofind.h block do not bracket "
                  "the root. x_a=%g x_b=%g f_a=%g f_b=%g Expanding region...",
                  ROOTFIND_X_a, ROOTFIND_X_b, ROOTFUNC_a, ROOTFUNC_b);
    double bracket_fac = 1.1;
    int bracket_iter = 0;
    do {
        double tmp = ROOTFIND_X_a; // let a be the lower value
        ROOTFIND_X_a = DMIN(ROOTFIND_X_a, ROOTFIND_X_b) / bracket_fac;
        ROOTFIND_X_b = DMAX(tmp, ROOTFIND_X_b) * bracket_fac;
        ROOTFUNC_a = ROOTFUNC(ROOTFIND_X_a);
        ROOTFUNC_b = ROOTFUNC(ROOTFIND_X_b);
        bracket_iter++;
    } while (ROOTFUNC_a * ROOTFUNC_b > 0 && bracket_iter < MAXITER);
    if ((bracket_iter == MAXITER) || isnan(ROOTFUNC_a) || isnan(ROOTFUNC_b)) {
        PRINT_WARNING("ERROR: Could not bracket root. x_a=%g x_b=%g f_a=%g f_b=%g\n", ROOTFIND_X_a, ROOTFIND_X_b, ROOTFUNC_a, ROOTFUNC_b);
        endrun(234528);
    }
}

if (fabs(ROOTFUNC_a) < fabs(ROOTFUNC_b)) { // in our convention 'a' represents
                                           // the bracket with the larger
                                           // residual
    double tmp = ROOTFUNC_a;
    ROOTFUNC_a = ROOTFUNC_b;
    ROOTFUNC_b = tmp;
    tmp = ROOTFIND_X_a;
    ROOTFIND_X_a = ROOTFIND_X_b;
    ROOTFIND_X_b = tmp;
}

double ROOTFIND_X_c = ROOTFIND_X_a, ROOTFUNC_c = ROOTFUNC_a;
int USED_BISECTION = 1, DO_BISECTION = 0, ROOTFIND_ITER = 0;
double ROOTFIND_X_c_old = ROOTFIND_X_c, ROOTFIND_X_new, ROOTFUNC_new = ROOTFUNC_c, ROOTFIND_X_error = 1e100, DELTA_TOL = 0.;

/* now we do a Brent 1973 method root-find */
do {
    ROOTFIND_X_new = 0;
    if ((ROOTFUNC_a != ROOTFUNC_c) && (ROOTFUNC_b != ROOTFUNC_c)) { // inverse quadratic interpolation
        ROOTFIND_X_new += ROOTFIND_X_a * ROOTFUNC_c * ROOTFUNC_b / (ROOTFUNC_a - ROOTFUNC_b) / (ROOTFUNC_a - ROOTFUNC_c);
        ROOTFIND_X_new += ROOTFIND_X_b * ROOTFUNC_c * ROOTFUNC_a / (ROOTFUNC_b - ROOTFUNC_a) / (ROOTFUNC_b - ROOTFUNC_c);
        ROOTFIND_X_new += ROOTFIND_X_c * ROOTFUNC_a * ROOTFUNC_b / (ROOTFUNC_c - ROOTFUNC_a) / (ROOTFUNC_c - ROOTFUNC_b);
    } else { // secant method
        ROOTFIND_X_new = (ROOTFIND_X_a * ROOTFUNC_b - ROOTFIND_X_b * ROOTFUNC_a) / (ROOTFUNC_b - ROOTFUNC_a);
    }
    DELTA_TOL = DMAX(fabs(ROOTFIND_ABS_X_tol), ROOTFIND_REL_X_tol * fabs(ROOTFIND_X_new)); // new absolute tolerance, incorporating relative tolerance

    DO_BISECTION = 0;
    double ROOTFIND_X_midpoint_a = 0.25 * (3 * ROOTFIND_X_a + ROOTFIND_X_b);
    if ((ROOTFIND_X_new < DMIN(ROOTFIND_X_midpoint_a, ROOTFIND_X_b)) || (ROOTFIND_X_new > DMAX(ROOTFIND_X_midpoint_a, ROOTFIND_X_b))) {
        DO_BISECTION = 1;
    } else { // accept the interpolation and bug out if it looks converged
        if (fabs(ROOTFIND_X_new - ROOTFIND_X_b) < DELTA_TOL) {
            break;
        }
    }
    if (USED_BISECTION) {
        if (fabs(ROOTFIND_X_new - ROOTFIND_X_b) >= 0.5 * fabs(ROOTFIND_X_c - ROOTFIND_X_b)) {
            DO_BISECTION = 1;
        }
        if (ROOTFIND_X_b != ROOTFIND_X_c) {
            if (fabs(ROOTFIND_X_b - ROOTFIND_X_c) < DELTA_TOL) {
                DO_BISECTION = 1;
            }
        }
    } else {
        if (fabs(ROOTFIND_X_new - ROOTFIND_X_b) >= 0.5 * fabs(ROOTFIND_X_c_old - ROOTFIND_X_c)) {
            DO_BISECTION = 1;
        }
        if (ROOTFIND_X_c_old != ROOTFIND_X_c) {
            if (fabs(ROOTFIND_X_c_old - ROOTFIND_X_c) < DELTA_TOL) {
                DO_BISECTION = 1;
            }
        }
    }
    if (DO_BISECTION) {
        // bisection in log space can help convergence for typical use cases; do
        // this if possible
        if ((ROOTFIND_X_b > 0) && (ROOTFIND_X_a > 0)) {
            ROOTFIND_X_new = sqrt(ROOTFIND_X_b * ROOTFIND_X_a);
        } else {
            ROOTFIND_X_new = 0.5 * (ROOTFIND_X_b + ROOTFIND_X_a);
        }
        USED_BISECTION = 1;
    } // bisection
    else {
        USED_BISECTION = 0;
    }
    ROOTFUNC_new = ROOTFUNC(ROOTFIND_X_new);
    if (ROOTFUNC_new == 0) {
        break;
    }

    ROOTFIND_X_c_old = ROOTFIND_X_c;
    ROOTFIND_X_c = ROOTFIND_X_b;
    ROOTFUNC_c = ROOTFUNC_b;
    if (ROOTFUNC_a * ROOTFUNC_new < 0) {
        ROOTFIND_X_b = ROOTFIND_X_new;
        ROOTFUNC_b = ROOTFUNC_new;
    } else {
        ROOTFIND_X_a = ROOTFIND_X_new;
        ROOTFUNC_a = ROOTFUNC_new;
    }

    if (fabs(ROOTFUNC_a) < fabs(ROOTFUNC_b)) {
        double tmp = ROOTFUNC_a;
        ROOTFUNC_a = ROOTFUNC_b;
        ROOTFUNC_b = tmp;
        tmp = ROOTFIND_X_a;
        ROOTFIND_X_a = ROOTFIND_X_b;
        ROOTFIND_X_b = tmp;
    }
    ROOTFIND_X_error = fabs(ROOTFIND_X_b - ROOTFIND_X_a);
    ROOTFIND_ITER++;
    if (ROOTFIND_ITER > MAXITER) {
        break;
    }
} while (ROOTFIND_X_error > DELTA_TOL);

#undef ROOTFUNC

#ifdef ROOTFIND_FUNCTION_INNER
#undef ROOTFIND_FUNCTION_INNER
#else
#undef ROOTFIND_FUNCTION
#endif
