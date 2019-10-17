/* here is where we call the core of the SIDM calculation for DM particle-particle interactions */
#ifdef DM_SIDM
{
    /* check if target+neighbor are an SIDM candidate, and against self-interaction */
    if( ((1 << local.Type) & (DM_SIDM)) && ((1 << P[j].Type) & (DM_SIDM)) && (local.ID != P[j].ID) && (local.dt_step <= P[j].dt_step))
    {
        if((local.dt_step==P[j].dt_step) && (local.ID > P[j].ID)) continue; // ensures interaction will only be calculated once for each pair //
        double h_si = 0.5*(kernel.h_i + kernel.h_j), m_si = 0.5*(local.Mass + P[j].Mass);
#ifdef GRAIN_COLLISIONS
        double prob = prob_of_grain_interaction(local.Grain_CrossSection_PerUnitMass , local.Mass, kernel.r, h_si, kernel.dv, local.dt_step, j);
#else
        double prob = prob_of_interaction(m_si, kernel.r, h_si, kernel.dv, local.dt_step);
#endif
        if(prob > 0.2) {out.dt_step_sidm = DMIN(out.dt_step_sidm , local.dt_step*(0.2/prob));} // timestep condition not being met as desired, warn code to lower timestep next turn //
        if (gsl_rng_uniform(random_generator) < prob)
        {
#ifdef WAKEUP
            if(!(TimeBinActive[P[j].TimeBin])) {if(WAKEUP*local.dt_step < P[j].dt_step) {PPPZ[j].wakeup=1; NeedToWakeupParticles_local = 1;}}
#endif
            double kick[3]; calculate_interact_kick(kernel.dv, kick, m_si);
            int k; for(k=0;k<3;k++) {P[j].Vel[k] += (local.Mass/m_si)*kick[k]; out.sidm_kick[k] -= (P[j].Mass/m_si)*kick[k];}
            out.si_count++; P[j].NInteractions++;
        }
    } // if((1 << ptype) & (DM_SIDM))
}
#endif
