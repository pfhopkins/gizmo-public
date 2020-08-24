#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

#if (SINGLE_STAR_TIMESTEPPING > 0)
// wraps around angle to the interval [0, 2pi)
double wrap_angle(double angle){
    if (angle > 2*M_PI)	return fmod(angle, 2*M_PI);
    else if (angle < 0) return 2*M_PI + fmod(angle, 2*M_PI);
    else return angle;
}

// Solve Kepler's equation to convert mean anomaly into eccentric anomaly
double eccentric_anomaly(double mean_anomaly, double ecc){
    double x0 = mean_anomaly;
    double err = 1e100;
    int iterations = 0;
    double twopi = 2*M_PI;
    while(fabs(err/twopi) > 1e-14 && iterations < 20){ // do Newton iterations
        err = (x0 - ecc*sin(x0) - mean_anomaly)/(1 - ecc*cos(x0)); 
        x0 -= err;
        x0 = wrap_angle(x0);
        iterations += 1;
    }
    return x0;
}


/* 
Advances the binary by timestep dt 
mode 0 - Just fill out the particle's kick and drift for the timestep, without doing the update
mode 1 - Actually update the binary separation and relative velocity. This should be done on the full-step drift.
*/

void kepler_timestep(int i, double dt, double kick_dv[3], double drift_dx[3], int mode){
    double h[3]; //Specific angular momentum vector
    double dr = sqrt(P[i].comp_dx[0]*P[i].comp_dx[0] + P[i].comp_dx[1]*P[i].comp_dx[1] + P[i].comp_dx[2]*P[i].comp_dx[2]); 
    double dv = sqrt(P[i].comp_dv[0]*P[i].comp_dv[0] + P[i].comp_dv[1]*P[i].comp_dv[1] + P[i].comp_dv[2]*P[i].comp_dv[2]);

    double dx_normalized[3] = {P[i].comp_dx[0]/dr, P[i].comp_dx[1]/dr, P[i].comp_dx[2]/dr};
    double dx_new[3], dv_new[3];
    double n_x[3]; // normalized Laplace-Runge-Lenz vector, just to get the unit vector along the major axis of the binary
    double n_y[3]; // normalized unit vector along the minor axis of the binary
    double norm, true_anomaly, mean_anomaly, ecc_anomaly, cos_true_anomaly,sin_true_anomaly;
    double x = 0, y =0, vx =0, vy = 0; // Coordinates in the frame aligned with the binary
    int k,l,m;
    double Mtot = P[i].Mass + P[i].comp_Mass;

    double specific_energy = .5*dv*dv - All.G * Mtot / dr;
    double semimajor_axis = -All.G * Mtot / (2*specific_energy);

    for(k=0; k<3; k++){ // dx cross dv to get specific angular momentum vector
        l = (k+1)%3; m = (k+2)%3;
        h[k] = P[i].comp_dx[l]*P[i].comp_dv[m] - P[i].comp_dx[m]*P[i].comp_dv[l];
    }
    
    double h2 = h[0]*h[0] + h[1]*h[1] + h[2]*h[2];
    double ecc = sqrt(1 + 2 * specific_energy * h2 / (All.G*All.G*Mtot*Mtot)); 
    for(k=0; k<3; k++){ // Get the LRL vector dv x h - GM dx/r
        l = (k+1)%3; m = (k+2)%3;
        n_x[k] = P[i].comp_dv[l]*h[m] - P[i].comp_dv[m]*h[l] - All.G * Mtot * dx_normalized[k]; // Worry about cancellation error for low eccentricity?
    }

    norm = sqrt(n_x[0]*n_x[0] + n_x[1]*n_x[1] + n_x[2]*n_x[2]);    
    for(k=0; k<3; k++) n_x[k] /= norm; // direction should be so that x points from periapsis to apoapsis

    for(k=0; k<3; k++, l=(k+1)%3, m=(k+2)%3){ // cross product of n_x with angular momentum to get a vector along the minor axis
        l = (k+1)%3; m = (k+2)%3;
        n_y[k] = n_x[l] * h[m] - n_x[m] * h[l];
        n_y[k] /= sqrt(h2);
    }

    // Transform to coordinates in the plane of the ellipse
    for(k=0; k<3; k++){
        x += P[i].comp_dx[k]*n_x[k];
        y += P[i].comp_dx[k]*n_y[k];
    }
    //printf("Kepler transform stuff x %g y %g nx %g %g %g ny %g %g %g h %g %g %g\n", x,y,n_x[0],n_x[1],n_x[2],n_y[0],n_y[1],n_y[2],h[0],h[1],h[2]);
    
    true_anomaly = wrap_angle(atan2(y,x));
    ecc_anomaly = wrap_angle(atan2(sqrt(1 - ecc*ecc) * sin(true_anomaly), ecc + cos(true_anomaly)));
    mean_anomaly = wrap_angle(ecc_anomaly - ecc * sin(ecc_anomaly));
    //printf("Kepler x %g y %g dr orig %g dv orig %g ecc_anomaly %g mean_anomaly %g true anomaly %g change in mean anomaly %g ID %d \n", x, y, dr, dv, ecc_anomaly, mean_anomaly, true_anomaly, (dt/P[i].min_bh_t_orbital * 2 * M_PI),P[i].ID);
    //Changes mean anomaly as time passes
    mean_anomaly -= dt/P[i].min_bh_t_orbital * 2 * M_PI;
    mean_anomaly = wrap_angle(mean_anomaly);
    //Get eccentric anomaly for new position
    ecc_anomaly = eccentric_anomaly(mean_anomaly, ecc);
    //Get sine and cosine of new true anomaly (we don't actually need the new value)
    sin_true_anomaly = sqrt(1.0-ecc*ecc)*sin(ecc_anomaly)/( 1 - ecc*cos(ecc_anomaly) );
    cos_true_anomaly = ( cos(ecc_anomaly) - ecc )/( 1 - ecc*cos(ecc_anomaly) );
    dr = semimajor_axis * (1-ecc*ecc)/(1+ecc*cos_true_anomaly);
    x = dr * cos_true_anomaly;
    y = dr * sin_true_anomaly;
   
    
    dv = sqrt(All.G * Mtot * (2/dr - 1/semimajor_axis)); // We conserve energy exactly

    double v_phi = -sqrt(h2) / dr; // conserving angular momentum
    double v_r = sqrt(DMAX(0, dv*dv - v_phi*v_phi));
    if(ecc_anomaly < M_PI) v_r = -v_r; // if radius is decreasing, make sure v_r is negative
    
    //relative velocities in the frame aligned with the ellipse:
    vx = v_phi * (-y/dr) + v_r * x/dr;
    vy = v_phi * (x/dr) + v_r * y/dr;

    // transform back to global coordinates
    double two_body_factor=-P[i].comp_Mass/Mtot;
    //printf("Kepler comp_dx %g %g %g  comp_dv %g %g %g ID %d \n", P[i].comp_dx[0],P[i].comp_dx[1],P[i].comp_dx[2],P[i].comp_dv[0],P[i].comp_dv[1], P[i].comp_dv[2], P[i].ID);
    for(k=0; k<3; k++){
        dx_new[k] = x * n_x[k] + y * n_y[k];
        dv_new[k] = vx * n_x[k] + vy * n_y[k];
        drift_dx[k] = (dx_new[k] - P[i].comp_dx[k]) * two_body_factor;
        kick_dv[k] = (dv_new[k] - P[i].comp_dv[k]) * two_body_factor;
        if(mode==1){ // if we want to do the actual self-consistent binary update
            P[i].comp_dx[k] = dx_new[k];
            P[i].comp_dv[k] = dv_new[k];
        }
    }
}

// Quantity needed for gravitational acceleration and jerk; mass / r^3 in Newtonian gravity
double gravfac(double r, double mass){
    if(r < All.ForceSoftening[5]) {
	double u = r / All.ForceSoftening[5];
	double h_inv = 1/All.ForceSoftening[5];
	return mass * kernel_gravity(u, h_inv, h_inv*h_inv*h_inv, 1);
    } else return mass / (r*r*r);
}

// quantity needed for the jerk, 3* mass/r^5 in Newtonian gravity
double gravfac2(double r, double mass)
{
    double hinv = 1. / All.ForceSoftening[5];
    return mass * kernel_gravity(r*hinv, hinv, hinv*hinv*hinv, 2);
}

// Computes the gravitational acceleration of a body at separation dx from a mass, accounting for softening
void grav_accel(double mass, double dx[3], double accel[3]){
    double fac, r = 0;
    int k;
    for(k=0; k<3; k++) r += dx[k]*dx[k];
    r = sqrt(r);
    fac = gravfac(r, mass); // mass / r^3 for Newtonian gravity
    for(k=0; k<3; k++) accel[k] = -dx[k] * fac;
}

// Computes the gravitational acceleration and time derivative of acceleration (the jerk) of a body at separation dx and relative velocity dv from a mass, accounting for softening
void grav_accel_jerk(double mass, double dx[3], double dv[3], double accel[3], double jerk[3]){
    double fac, fac2, r = 0, dv_dot_dx = 0;
    int k;
    for(k=0; k<3; k++) {
	r += dx[k]*dx[k];
	dv_dot_dx += dv[k] * dx[k];
    }
    
    r = sqrt(r);
    
    fac = gravfac(r, mass); // mass / r^3 for Newtonian gravity
    fac2 = gravfac2(r, mass);
    for(k=0; k<3; k++) {
	accel[k] = All.G * (-dx[k] * fac);
	jerk[k] = All.G * (-dv[k] * fac + dv_dot_dx * fac2 * dx[k]);
    }
}

// Perform a 4th order Hermite timestep for the softened Kepler problem, evolving the orbital separation dx and relative velocity dv
void hermite_step(double mass, double dx[3], double dv[3], double dt){
    double old_accel[3], old_jerk[3], old_dx[3], old_dv[3], accel[3], jerk[3];
    double dt2 = dt*dt, dt3 = dt2 * dt;
    int k;
    grav_accel_jerk(mass, dx, dv, old_accel, old_jerk);

    // Predictor step
    for(k=0; k<3; k++){
	old_dx[k] = dx[k];
	old_dv[k] = dv[k];
	dx[k] += dv[k] * dt + 0.5*old_accel[k]*dt2 + old_jerk[k]*dt3/6;
	dv[k] += old_accel[k] * dt + 0.5*old_jerk[k]*dt2;
    }

    grav_accel_jerk(mass, dx, dv, accel, jerk);

    for(k=0; k<3; k++){
	dv[k] = old_dv[k] + 0.5*(accel[k] + old_accel[k]) * dt + (old_jerk[k] - jerk[k]) * dt2/12;
	dx[k] = old_dx[k] + 0.5*(dv[k] + old_dv[k]) * dt + (old_accel[k] - accel[k])*dt2/12;
    }
}

/* 
Advances the binary by timestep dt
mode 0 - Just fill out the particle's kick and drift for the timestep, without doing the update
mode 1 - Actually update the binary separation and relative velocity. This should be done on the full-step drift.
*/
void odeint_super_timestep(int i, double dt_super, double kick_dv[3], double drift_dx[3], int mode)
{
    double t = 0, total_mass = P[i].comp_Mass + P[i].Mass, dt;
    double dx_old[3], dv_old[3], dx[3], dv[3], vSqr, rSqr;
    int k;
    for(k=0; k<3; k++){
	dx_old[k] = P[i].comp_dx[k];  // saving the old separation
	dx[k] = -P[i].comp_dx[k];  // this one gets evolved; note sign change from comp_dx to the effective 1-body problem
	dv_old[k] = P[i].comp_dv[k];  // saving the old separation 
	dv[k] = -P[i].comp_dv[k];  // this one gets evolved; note sign change from comp_dx to the effective 1-body problem
    }

    while(t < dt_super){
	vSqr = rSqr = 0;
	// Determine timestep adaptively; tuned here to give 1% energy error over 10^5 orbits for a 0.9 eccentricty binary
	for(k=0; k<3; k++){
	    vSqr += dv[k]*dv[k];
	    rSqr += dx[k]*dx[k];
	}
    double r_effective = KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER * All.ForceSoftening[5];
    rSqr += r_effective*r_effective;
	dt = DMIN(0.1/(sqrt(vSqr/rSqr) + sqrt((All.G * total_mass)/sqrt(rSqr*rSqr*rSqr))), dt_super-t); // harmonic mean of approach time and freefall time
    // could swap in any integration scheme you want here; default to Hermite
	hermite_step(total_mass, dx, dv, dt);
	t += dt;
    }

    double two_body_factor=-P[i].comp_Mass/total_mass;

    for(k=0; k<3; k++){
        drift_dx[k] = (-dx[k] - P[i].comp_dx[k]) * two_body_factor;
        kick_dv[k] = (-dv[k] - P[i].comp_dv[k]) * two_body_factor;
        if(mode==1){ // if we want to do the actual self-consistent binary update
            P[i].comp_dx[k] = -dx[k];
            P[i].comp_dv[k] = -dv[k];
        }
    }
    
}

#endif
