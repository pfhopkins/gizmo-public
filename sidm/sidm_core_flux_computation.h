/* here is where we call the core of the SIDM calculation for DM particle-particle interactions */
#ifdef DM_SIDM
{
    /* check if target+neighbor are an SIDM candidate, and against self-interaction */
    if( ((1 << local.Type) & (DM_SIDM)) && ((1 << P[j].Type) & (DM_SIDM)) && (local.ID != P[j].ID) && (local.dt_step <= P[j].dt_step))
    {
        if((local.dt_step==P[j].dt_step) && (local.ID > P[j].ID)) continue; // ensures interaction will only be calculated once for each pair //
        double h_si = 0.5*(kernel.h_i + kernel.h_j), m_si = 0.5*(local.Mass + P[j].Mass);
        double prob = prob_of_interaction(m_si, kernel.r, h_si, kernel.dv, local.dt_step);
        if(prob > 0.2)
        {
            integertime dt_step_guess = DMIN(local.dt_step, out.dt_step_sidm);
            while(prob_of_interaction(m_si,kernel.r,h_si,kernel.dv,dt_step_guess) > 0.2) {dt_step_guess/=2;}
            out.dt_step_sidm = dt_step_guess;
        }
        if (gsl_rng_uniform(random_generator) < prob)
        {
            double kick[3]; calculate_interact_kick(kernel.dv, kick, m_si);
            int k; for(k=0;k<3;k++) {P[j].Vel[k] += (local.Mass/m_si)*kick[k]; out.sidm_kick[k] -= (P[j].Mass/m_si)*kick[k];}
            out.si_count++; P[j].NInteractions++;
        }
    } // if((1 << ptype) & (DM_SIDM))
}
#endif
