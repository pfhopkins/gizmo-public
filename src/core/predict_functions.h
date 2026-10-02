/* predict_functions.h — Canonical KOKKOS_INLINE_FUNCTION implementation of
 * evaluate_NH_from_GradRho.  Single source of truth for both CPU and GPU.
 *
 * Proto.h has an inline Vec3<MyFloat> wrapper that forwards to this.
 *
 * Include order: after allvars.h (for All, MyFloat).  This header includes nothing of
 * its own: core/predict.cc re-includes it with non-inline linkage to provide the host
 * symbols, and anything included from here would be re-included that way too. */
#pragma once

#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif

KOKKOS_INLINE_FUNCTION
double evaluate_NH_from_GradRho(MyFloat gradrho[3], double rkern, double rho, double numngb_ndim, double include_h, int target, struct particle_data *pp)
{
    double gradrho_mag=0;
    if(rho>0)
    {
#ifdef RT_USE_TREECOL_FOR_NH
        gradrho_mag = include_h * rho * rkern / numngb_ndim; if(target>=0) {gradrho_mag += pp[target].SigmaEff;}
#else
        gradrho_mag = sqrt(gradrho[0]*gradrho[0]+gradrho[1]*gradrho[1]+gradrho[2]*gradrho[2]);
        if(gradrho_mag > 0) {gradrho_mag = rho*rho/gradrho_mag;} else {gradrho_mag=0;}
        if(include_h > 0) if(numngb_ndim > 0) gradrho_mag += include_h * rho * rkern / numngb_ndim;
#endif
    }
    return gradrho_mag * All.cf_a2inv;
}

/* calculate_face_area_for_cartesian_mesh — migrated from predict.cc to fix
 * #20011-D (host-only fn called from KOKKOS_INLINE_FUNCTION
 * compute_finitevol_faces template under HYDRO_REGULAR_GRID). Body uses only
 * All.cf_atime (mirror-safe), std::max, fabs — all device-callable. */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
KOKKOS_INLINE_FUNCTION
double calculate_face_area_for_cartesian_mesh(const Vec3<double>& dp, double rinv, double l_side, Vec3<double>& Face_Area_Vec)
{
    Face_Area_Vec = {}; double Face_Area_Norm;
#if (NUMDIMS==1)
    Face_Area_Norm = 1; Face_Area_Vec[0] = Face_Area_Norm * dp[0]/fabs(dp[0]);
#elif (NUMDIMS==2)
    if(fabs(dp[0]) > fabs(dp[1])) {Face_Area_Vec[0] = Face_Area_Norm = DMAX(0.,l_side-fabs(dp[1])) * dp[0]/fabs(dp[0]) * All.cf_atime;} else {Face_Area_Vec[1] = Face_Area_Norm = DMAX(0.,l_side-fabs(dp[0])) * dp[1]/fabs(dp[1]) * All.cf_atime;}
#else
    Vec3<double> dp_abs = {fabs(dp[0]), fabs(dp[1]), fabs(dp[2])};
    int kdir;
    if((dp_abs[0]>=dp_abs[1])&&(dp_abs[0]>=dp_abs[2])) {kdir=0;} else if ((dp_abs[1]>=dp_abs[0])&&(dp_abs[1]>=dp_abs[2])) {kdir=1;} else {kdir=2;}
    Face_Area_Norm=1; for(int k=0;k<3;k++) {if(k!=kdir) {Face_Area_Norm *= DMAX(0.,l_side-dp_abs[k]) * All.cf_atime*All.cf_atime;}}
    Face_Area_Vec[kdir] = Face_Area_Norm * dp[kdir]/fabs(dp[kdir]);
#endif
    return fabs(Face_Area_Norm);
}
#endif

/* Get_Particle_Expected_Area — migrated from predict.cc to fix #20011-D
 * (host-only fn called from KOKKOS_INLINE_FUNCTION compute_finitevol_faces
 * under SLOPE_LIMITER_TOLERANCE==0.
 * Pure-compute function of `h`, dimension-dependent. */
KOKKOS_INLINE_FUNCTION
double Get_Particle_Expected_Area(double h)
{
#if (NUMDIMS == 1)
    return 2;
#endif
#if (NUMDIMS == 2)
    return 2 * M_PI * h;
#endif
#if (NUMDIMS == 3)
    return 4 * M_PI * h * h;
#endif
}


#ifdef DIVBCLEANING_DEDNER
KOKKOS_INLINE_FUNCTION
double Get_Gas_PhiField_P(int i_particle_id, struct particle_data *pp, struct gas_cell_data *cell)
{
    //return cell[i_particle_id].PhiPred * cell[i_particle_id].Density / pp[i_particle_id].Mass; // volumetric phy-flux (requires extra term compared to mass-based flux)
    return cell[i_particle_id].PhiPred / pp[i_particle_id].Mass; // mass-based phi-flux
}
#endif /* DIVBCLEANING_DEDNER */


#ifdef DIVBCLEANING_DEDNER
KOKKOS_INLINE_FUNCTION
double Get_Gas_PhiField_DampingTimeInv_P(int i_particle_id, struct particle_data *pp, struct gas_cell_data *cell)
{
    /* this timescale should always be returned as a -physical- time */
#ifdef HYDRO_SPH
    /* PFH: add simple damping (-phi/tau) term */
    double damping_tinv = 0.5 * All.DivBcleanParabolicSigma * (cell[i_particle_id].MaxSignalVel / (All.cf_atime*pp[i_particle_id].Get_Particle_Size()));
#else
    double damping_tinv;
#ifdef SELFGRAVITY_OFF
    damping_tinv = All.DivBcleanParabolicSigma * All.FastestWaveSpeed / (All.cf_atime*pp[i_particle_id].Get_Particle_Size()); // fastest wavespeed has units of [vphys]
    //double damping_tinv = All.DivBcleanParabolicSigma * All.FastestWaveDecay * All.cf_a2inv; // no improvement over fastestwavespeed; decay has units [vphys/rphys]
#else
    // only see a small performance drop from fastestwavespeed above to maxsignalvel below, despite the fact that below is purely local (so allows more flexible adapting to high dynamic range)
    damping_tinv = 0.0;

    if(pp[i_particle_id].KernelRadius > 0)
    {
        double h_eff = pp[i_particle_id].Get_Particle_Size();
        double vsig2 = 0.5 * fabs(cell[i_particle_id].MaxSignalVel);
        double phi_B_eff = 0.0;
        if(vsig2 > 0) {phi_B_eff = Get_Gas_PhiField_P(i_particle_id, pp, cell) / (All.cf_atime * vsig2);}
        double vsig1 = 0.0;
        if(cell[i_particle_id].Density > 0)
        {
            vsig1 = sqrt( cell[i_particle_id].effective_soundspeed()*cell[i_particle_id].effective_soundspeed() +
                 (1. / All.cf_atime) *
                 (cell[i_particle_id].Bfield().norm_sq() +
                  phi_B_eff*phi_B_eff) / cell[i_particle_id].Density );
        }
        vsig1 = DMAX(vsig1, vsig2);
        vsig2 = 0.0;
        vsig2 = cell[i_particle_id].Gradients.Velocity.frobenius_norm();
        vsig2 = 3.0 * h_eff * DMAX( vsig2, fabs(pp[i_particle_id].Particle_DivVel)) / All.cf_atime;
        double prefac_fastest = 0.1;
        double prefac_tinv = 0.5;
        double area_0 = 0.1;
#ifdef MHD_CONSTRAINED_GRADIENT
        prefac_fastest = 1.0;
        prefac_tinv = 2.0;
        area_0 = 0.05;
        vsig2 *= 5.0;
        if(cell[i_particle_id].FlagForConstrainedGradients <= 0) prefac_tinv *= 30;
#endif
        prefac_tinv *= sqrt(1. + cell[i_particle_id].ConditionNumber/100.);
        double area = fabs(cell[i_particle_id].Face_Area[0]) + fabs(cell[i_particle_id].Face_Area[1]) + fabs(cell[i_particle_id].Face_Area[2]);
        area /= Get_Particle_Expected_Area(pp[i_particle_id].KernelRadius);
        prefac_tinv *= (1. + area/area_0)*(1. + area/area_0);

        double vsig_max = DMAX( DMAX(vsig1,vsig2) , prefac_fastest * All.FastestWaveSpeed );
        damping_tinv = prefac_tinv * All.DivBcleanParabolicSigma * (vsig_max / (All.cf_atime * h_eff));
    }
#endif
#endif
    return damping_tinv;
}
#endif /* DIVBCLEANING_DEDNER */
/* How many times box_wrap_position_to_primary_image folded the x axis, in the order it folded it: first up
 * from below zero, then down from at or above the box.  Kept as two counts, not one net count, because a
 * position just below zero can fold up onto the box length itself and then straight back down, and the
 * shearing offsets a caller applies per fold do not cancel exactly in floating point. */
struct box_wrap_x_folds {int up_from_below; int down_from_above;};

/* The primary-box image of a position, folded exactly as the box wrapping folds a particle: axis by axis,
 * x first, each axis by whole box lengths into [0, box).  In a shearing box each x fold also moves the
 * shearing coordinate by the boundary's position offset (BOX_SHEARING > 1) before that axis is itself
 * folded, and changes the velocity across that face by the boundary's velocity offset, which this does not
 * apply: it returns the x folds, and the caller applies the offset once per fold, the up folds first, to
 * whatever velocities it carries.  Outside a periodic build the position is returned unchanged with no
 * folds.  The position must be finite. */
KOKKOS_INLINE_FUNCTION
struct box_wrap_x_folds box_wrap_position_to_primary_image(Vec3<MyDouble> &pos)
{
    struct box_wrap_x_folds folds = {0, 0};
#ifdef BOX_PERIODIC
    const double boxsize[3] = {boxSize_X, boxSize_Y, boxSize_Z};
    for(int j = 0; j < 3; j++)
    {
        while(pos[j] < 0)
        {
            pos[j] += boxsize[j];
            if(j == 0)
            {
                folds.up_from_below++;
#if defined(BOX_SHEARING) && (BOX_SHEARING > 1)
                pos[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Pos_Offset;
#endif
            }
        }
        while(pos[j] >= boxsize[j])
        {
            pos[j] -= boxsize[j];
            if(j == 0)
            {
                folds.down_from_above++;
#if defined(BOX_SHEARING) && (BOX_SHEARING > 1)
                pos[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Pos_Offset;
#endif
            }
        }
    }
#endif
    return folds;
}
