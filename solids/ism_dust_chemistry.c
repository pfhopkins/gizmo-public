#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/* This module collects the live ism dust chemistry modules developed by Caleb Choban in Choban et al., 2022.
    Written by C. Choban, reorganized and collected by PFH.
 */

#if defined(GALSF_ISMDUSTCHEM_MODEL)


#define GALSF_ISMDUSTCHEM_VAR_IRON_INCL_FRAC 0.7 /* assumed fraction of iron dust mass locked as inclusions in silicates, this scales with the total fraction of silicate formed vs maximum amount of possible silicate dust */

/* routine to give yields for dust for different types of SNe (Ia & II) followed in-code */
void ISMDustChem_get_SNe_dust_yields(double *yields, int i, double t_gyr, int SNeIaFlag, double Msne)
{
    double dust_yields[NUM_ISMDUSTCHEM_ELEMENTS]={0}, sources_yields[NUM_ISMDUSTCHEM_SOURCES]={0}, species_yields[NUM_ISMDUSTCHEM_SPECIES]={0}; double SNeIa_age = 0.03753; int k,source_key=1;
    if(t_gyr < SNeIa_age) {source_key=2;} // 1=1a, 2=II
    for(k=0;k<NUM_ISMDUSTCHEM_ELEMENTS+NUM_ISMDUSTCHEM_SOURCES+NUM_ISMDUSTCHEM_SPECIES;k++) {yields[k+NUM_METAL_SPECIES]=0;} // initialize yields to null
    if(GALSF_ISMDUSTCHEM_MODEL & 1) {
        double C_condens_eff = 0.5, other_condens_eff = 0.8;
        dust_yields[2] = C_condens_eff * yields[2];         // C
        dust_yields[6] = other_condens_eff * yields[6];     // Mg
        dust_yields[7] = other_condens_eff * yields[7];     // Si
        dust_yields[10] = other_condens_eff * yields[10];   // Fe
        dust_yields[4] = 16 * (dust_yields[6]/All.ISMDustChem_AtomicMassTable[6] + dust_yields[7]/All.ISMDustChem_AtomicMassTable[7] + dust_yields[10]/All.ISMDustChem_AtomicMassTable[10]); // O
        if(dust_yields[4]>yields[4]) {dust_yields[4]=yields[4];} // Just in case there's not enough O
        for(k=2;k<NUM_ISMDUSTCHEM_ELEMENTS;k++)  dust_yields[0] += dust_yields[k]; // Fraction of yields that is dust
        for(k=0;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {yields[k+NUM_METAL_SPECIES]=dust_yields[k];}
        yields[NUM_METAL_SPECIES+NUM_ISMDUSTCHEM_ELEMENTS+source_key] = dust_yields[0]; // total yield goes to the source term of this type
        return; // all done, if only using this model
    } // below follows species model, will be default if above not set
    
    double SNeII_sil_cond = 0.00035, SNeII_C_cond = 0.15, SNeII_SiC_cond = 0.0003, SNeII_Fe_cond = 0.001, SNeI_Fe_cond = 0.005, sil_elem_abund[GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES];
    int missing_element = 0, key_elem = 0;
    // For each dust species find the key element and condense a fraction of that element into dust
    if(t_gyr < SNeIa_age)
    {
        /******** SILICATE ********/
        // first check that there are non-zero amounts of all elements required to make dust species
        for(k=0;k<GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES;k++) {if(yields[All.ISMDustChem_SilicateMetallicityFieldIndexTable[k]] <= 0) {missing_element = 1;}}
        if(!missing_element)
        {
            // Find the key element for silicate dust
            for (k=0;k<GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES;k++) {sil_elem_abund[k] = yields[All.ISMDustChem_SilicateMetallicityFieldIndexTable[k]] / (All.ISMDustChem_AtomicMassTable[All.ISMDustChem_SilicateMetallicityFieldIndexTable[k]] * All.ISMDustChem_SilicateNumberOfAtomsTable[k]);}
            for (k=1;k<GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES;k++) {if (sil_elem_abund[key_elem]>sil_elem_abund[k]) key_elem = k;}
            species_yields[0] = SNeII_sil_cond * yields[All.ISMDustChem_SilicateMetallicityFieldIndexTable[key_elem]] * All.ISMDustChem_EffectiveSilicateDustAtomicWeight / (All.ISMDustChem_SilicateNumberOfAtomsTable[key_elem] * All.ISMDustChem_AtomicMassTable[All.ISMDustChem_SilicateMetallicityFieldIndexTable[key_elem]]);
            for (k=0;k<GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES;k++) {dust_yields[All.ISMDustChem_SilicateMetallicityFieldIndexTable[k]] += species_yields[0] * All.ISMDustChem_SilicateNumberOfAtomsTable[k] * All.ISMDustChem_AtomicMassTable[All.ISMDustChem_SilicateMetallicityFieldIndexTable[k]] / All.ISMDustChem_EffectiveSilicateDustAtomicWeight;}
        }
        /******** CARBONACEOUS ********/
        species_yields[1] = SNeII_C_cond * yields[2];
        dust_yields[2] += species_yields[1];
        /******** SILICONE CARBIDE ********/
        if (yields[2]>0 && yields[7]>0)
        {
            if (yields[7]/All.ISMDustChem_AtomicMassTable[7] < yields[2]/All.ISMDustChem_AtomicMassTable[2]) key_elem = 7;
            else key_elem = 2;
            species_yields[2] = SNeII_SiC_cond * yields[key_elem] * ((All.ISMDustChem_AtomicMassTable[2] + All.ISMDustChem_AtomicMassTable[7]) / All.ISMDustChem_AtomicMassTable[key_elem]);
            dust_yields[2] += species_yields[2] * All.ISMDustChem_AtomicMassTable[2] / (All.ISMDustChem_AtomicMassTable[2] + All.ISMDustChem_AtomicMassTable[7]);
            dust_yields[7] += species_yields[2] * All.ISMDustChem_AtomicMassTable[7] / (All.ISMDustChem_AtomicMassTable[2] + All.ISMDustChem_AtomicMassTable[7]);
        }
        /******** METALLIC IRON ********/
        species_yields[3] = SNeII_Fe_cond * yields[10];
        dust_yields[10] += species_yields[3];
    }
    else
    {
        // Only a little bit of metallic iron dust from SNIa
        species_yields[3] = SNeI_Fe_cond * yields[10];
        dust_yields[10] = species_yields[3];
    }
    for(k=2;k<NUM_ISMDUSTCHEM_ELEMENTS;k++)  dust_yields[0] += dust_yields[k]; // Fraction of yields that is dust
    for(k=0;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {yields[k+NUM_METAL_SPECIES]=dust_yields[k];}
    yields[NUM_METAL_SPECIES+NUM_ISMDUSTCHEM_ELEMENTS+source_key] = dust_yields[0]; // total yield goes to the source term of this type
    for(k=0;k<NUM_ISMDUSTCHEM_SPECIES;k++) {yields[k+NUM_METAL_SPECIES+NUM_ISMDUSTCHEM_ELEMENTS+NUM_ISMDUSTCHEM_SOURCES]=species_yields[k];}
}


/* routine to give the dust yields for AGB winds (currently no dust yield assumed for stars younger than AGB age from continuous mass-loss, i.e. O/B winds) */
void ISMDustChem_get_wind_dust_yields(double *yields, int i)
{
    double dust_yields[NUM_ISMDUSTCHEM_ELEMENTS]={0}, sources_yields[NUM_ISMDUSTCHEM_SOURCES]={0}, species_yields[NUM_ISMDUSTCHEM_SPECIES]={0}; int k,source_key=3;
    for(k=0;k<NUM_ISMDUSTCHEM_ELEMENTS+NUM_ISMDUSTCHEM_SOURCES+NUM_ISMDUSTCHEM_SPECIES;k++) {yields[k+NUM_METAL_SPECIES]=0;} // initialize yields to null
    double transition_age = 0.03753, star_age = evaluate_stellar_age_Gyr(i); // Assume AGB dust production stars at SNe II to SNe Ia transition. This limits AGB stars with mass < ~8 solar masses
    if(star_age <= transition_age) {return;} // no yield here if too young, otherwise continue
    if(GALSF_ISMDUSTCHEM_MODEL & 1) {
        double condens_eff = 0.8;
        if((yields[2]/All.ISMDustChem_AtomicMassTable[2])/(yields[4]/All.ISMDustChem_AtomicMassTable[4]) > 1.0) // AGB stars with abundace ratio C/O > 1 only produce carbonacous dust
        {
            dust_yields[2] = yields[2] - 0.75*yields[4]; dust_yields[0] = dust_yields[2]; // C
        } else { // AGB stars with abundance C/O < 1 produce general silicate dust
            dust_yields[6] = condens_eff * yields[6]; // Mg
            dust_yields[7] = condens_eff * yields[7]; // Si
            dust_yields[10] = condens_eff * yields[10]; // Fe
            dust_yields[4] = 16 * (dust_yields[6]/All.ISMDustChem_AtomicMassTable[6] + dust_yields[7]/All.ISMDustChem_AtomicMassTable[7] + dust_yields[10]/All.ISMDustChem_AtomicMassTable[10]); // O
            // Check to make sure we dont produce too much O dust given the leftover O dust not in CO
            if (dust_yields[4] > yields[4]-(4./3.*yields[2])) {dust_yields[4] = yields[4]-(4./3.*yields[2]);}
            for(k=2;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {dust_yields[0]+=dust_yields[k];}
        }
        for(k=0;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {yields[k+NUM_METAL_SPECIES]=dust_yields[k];}
        yields[NUM_METAL_SPECIES+NUM_ISMDUSTCHEM_ELEMENTS+source_key] = dust_yields[0]; // total yield goes to the source term of this type
        return; // end routine
    } // below follows species model, and will be default if above not set
    double dt,Z,elem_yield,wind_rate;
    dt=GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i)*UNIT_TIME_IN_GYR;
    Z = Z_for_stellar_evol(i);
    
    // Take difference in cumulative dust production between start and end time to get estimate of instantaneous dust injection rate (M_solar/Gyr)
    for (k=0;k<NUM_ISMDUSTCHEM_SPECIES;k++) {
        species_yields[k] = (cumulative_AGB_dust_returns(k,(star_age+dt)*1E3,Z)-cumulative_AGB_dust_returns(k,star_age*1E3,Z))/dt;
        species_yields[k] = DMAX(0.,species_yields[k]); // Deal with unphysical values which can result near boundaries in fit
    }
    // All done if no dust is produced
    if (species_yields[0]+species_yields[1]+species_yields[2]+species_yields[3]>0.)
    {
        // Now convert from instantaneous dust injection rates to dust yields using instantaneous wind rate
        wind_rate=0.41987*pow(star_age,-1.1)/(12.9-log(star_age));
        if(star_age < 0.033) {wind_rate *= 0.01 + calculate_relative_light_to_mass_ratio_from_imf(star_age,i,1);} // late-time independent of massive stars
        wind_rate *= All.StellarMassLoss_Rate_Renormalization;
        
        for (k=0;k<NUM_ISMDUSTCHEM_SPECIES;k++) {
            species_yields[k] = species_yields[k]/wind_rate;
            species_yields[k] = DMAX(0.,species_yields[k]);
        }
        
        // Now check to make sure there are enough metals in the wind to produce the dust since the metal and dust yields are calculated separately
        // If not renorm dust species which are made up of the element in question
        // Check C
        elem_yield = species_yields[1] + species_yields[2] * All.ISMDustChem_AtomicMassTable[2] / (All.ISMDustChem_AtomicMassTable[2] + All.ISMDustChem_AtomicMassTable[7]);
        if (elem_yield > yields[2])
        {
            species_yields[1] *= yields[2]/elem_yield;
            species_yields[2] *= yields[2]/elem_yield;
        }
        // Check O
        elem_yield = species_yields[0] * All.ISMDustChem_AtomicMassTable[4] * All.ISMDustChem_SilicateNumberOfAtomsTable[0] / All.ISMDustChem_EffectiveSilicateDustAtomicWeight;
        if (elem_yield > yields[4])
        {
            species_yields[0] *= yields[4]/elem_yield;
        }
        // Check Mg
        elem_yield = species_yields[0] * All.ISMDustChem_AtomicMassTable[6] * All.ISMDustChem_SilicateNumberOfAtomsTable[1] / All.ISMDustChem_EffectiveSilicateDustAtomicWeight;
        if (elem_yield > yields[6])
        {
            species_yields[0] *= yields[6]/elem_yield;
        }
        // Check Si
        elem_yield = species_yields[0] * All.ISMDustChem_AtomicMassTable[7] * All.ISMDustChem_SilicateNumberOfAtomsTable[2] / All.ISMDustChem_EffectiveSilicateDustAtomicWeight + species_yields[2] * All.ISMDustChem_AtomicMassTable[7] / (All.ISMDustChem_AtomicMassTable[2] + All.ISMDustChem_AtomicMassTable[7]);
        if (elem_yield > yields[7])
        {
            species_yields[0] *= yields[7]/elem_yield;
            species_yields[2] *= yields[7]/elem_yield;
        }
        // Check Fe
        if(GALSF_ISMDUSTCHEM_MODEL & 4) {
            if (species_yields[3] > yields[10]) {species_yields[3] = yields[10];} // Fe is only in free-flying iron, assume no iron inclusions in stellar dust
        } else { // Fe is in free-flying iron and silicates
            elem_yield = species_yields[0] * All.ISMDustChem_AtomicMassTable[10] * All.ISMDustChem_SilicateNumberOfAtomsTable[3] / All.ISMDustChem_EffectiveSilicateDustAtomicWeight + species_yields[3];
            if(elem_yield > yields[10]) {species_yields[0] *= yields[10]/elem_yield; species_yields[3] *= yields[10]/elem_yield;}
        }
        // Now convert from dust species to dust elemental mass
        // silicates
        for (k=0;k<GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES;k++)
        {
            dust_yields[All.ISMDustChem_SilicateMetallicityFieldIndexTable[k]] += species_yields[0] * All.ISMDustChem_SilicateNumberOfAtomsTable[k] * All.ISMDustChem_AtomicMassTable[All.ISMDustChem_SilicateMetallicityFieldIndexTable[k]] / All.ISMDustChem_EffectiveSilicateDustAtomicWeight;
        }
        // carbonaceous
        dust_yields[2] += species_yields[1];
        // SiC
        dust_yields[2] += species_yields[2] * All.ISMDustChem_AtomicMassTable[2] / (All.ISMDustChem_AtomicMassTable[2] + All.ISMDustChem_AtomicMassTable[7]);
        dust_yields[7] += species_yields[2] * All.ISMDustChem_AtomicMassTable[7] / (All.ISMDustChem_AtomicMassTable[2] + All.ISMDustChem_AtomicMassTable[7]);
        // iron
        dust_yields[10] += species_yields[3];
        
        for (k=2;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {dust_yields[0] += dust_yields[k];}
    }
    for(k=0;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {yields[k+NUM_METAL_SPECIES]=dust_yields[k];}
    yields[NUM_METAL_SPECIES+NUM_ISMDUSTCHEM_ELEMENTS+source_key] = dust_yields[0]; // total yield goes to the source term of this type
    for(k=0;k<NUM_ISMDUSTCHEM_SPECIES;k++) {yields[k+NUM_METAL_SPECIES+NUM_ISMDUSTCHEM_ELEMENTS+NUM_ISMDUSTCHEM_SOURCES]=species_yields[k];}
}



/* Simple fit to cumulative AGB dust production for a Kroupa IMF stellar population only with specific metallicities and stellar ages (assuming stars become AGBs at the ends of the main sequence lifetime) derived from AGB data table in Zhukovska+(08) */
double specific_Z_AGB_dust(int dust_type, double star_age, int z_bound)
{
    /* dust_type: 0 = silicate, 1 = carbon, 2 = silicon carbide, 3 = metallic iron */
    /* z_bound: 0 = 2*Z_solar, 1 = Z_solar, 2 = 0.4*Z_solar, 3 = 0.2*Z_solar, 4 = 0.05*Z_solar*/
    double cum_return;
    double logt = log10(star_age);
    switch(dust_type) {
        case 0:
            if (z_bound == 0) {
                if (star_age < 284) {cum_return = 1.77E-4*logt - 2.87E-4;}
                else if (star_age < 1244) {cum_return = 3.03E-5*(logt - 2.45) + 1.47E-4;}
                else {cum_return = 5.93E-4*(logt - 3.09) + 1.67E-4;}
            }
            else if (z_bound == 1) {
                if (star_age < 295) {cum_return = 4.18E-5*logt - 6.78E-5;}
                else if (star_age < 1808) {cum_return = 1.72E-6*(logt - 2.47) + 3.55E-5;}
                else {cum_return = 1.05E-4*(logt - 3.26) + 3.68E-5;}
            }
            else if (z_bound == 2) {
                if (star_age < 286) {cum_return = 5.00E-6*logt - 8.01E-6;}
                else if (star_age < 3948) {cum_return = 1.85E-9*(logt - 2.46) + 4.25E-6;}
                else {cum_return = 3.35E-7*(logt - 3.60) + 4.26E-6;}
            }
            else if (z_bound == 3) {
                if (star_age < 269) {cum_return = 1.04E-8*logt - 1.64E-08;}
                else if (star_age < 1560) {cum_return = 2.69E-10*(logt - 2.43) + 8.82E-9;}
                else {cum_return = 3.05E-19*(logt - 3.19) + 9.03E-9;}
            }
            else if (z_bound == 4) {
                if (star_age < 147) {cum_return = 5.80E-11*logt - 9.10E-11;}
                else if (star_age < 252) {cum_return = 4.10E-11*(logt - 2.17) + 3.47E-11;}
                else {cum_return = 6.05E-14*(logt - 2.40) + 4.43E-11;}
            }
            break;
            
        case 1:
            if (z_bound == 0) {
                if (star_age < 262) {cum_return = 8.10E-7*logt - 1.54E-6;}
                else if (star_age < 840) {cum_return = 4.00E-4*(logt - 2.42) + 4.23E-7;}
                else {cum_return = 7.37E-6*(logt - 2.92) + 2.02E-4;}
            }
            else if (z_bound == 1) {
                if (star_age < 305) {cum_return = 3.66E-6*logt - 6.75E-6;}
                else if (star_age < 1250) {cum_return = 3.71E-4*(logt - 2.48) +  2.34E-6;}
                else {cum_return = 8.67E-6*(logt - 3.10) + 2.29E-4;}
            }
            else if (z_bound == 2) {
                if (star_age < 367) {cum_return = 1.27E-5*logt - 2.34E-5;}
                else if (star_age < 2329) {cum_return = 4.24E-4*(logt - 2.56) + 9.27E-6;}
                else {cum_return = 7.55E-5*(logt - 3.37) + 3.49E-4;}
            }
            else if (z_bound == 3) {
                if (star_age < 344) {cum_return = 1.40E-5*logt - 2.53E-5;}
                else if (star_age < 3105) {cum_return = 4.44E-4*(logt - 2.54) + 1.03e-5;}
                else {cum_return = 9.59e-5*(logt - 3.49) + 4.35E-4;}
            }
            else if (z_bound == 4) {
                if (star_age < 280) {cum_return = 8.48E-6*logt - 1.43E-5;}
                else if (star_age < 4504) {cum_return = 3.10E-4*(logt - 2.45) + 6.47E-6;}
                else {cum_return = 5.73E-5*(logt - 3.65) + 3.80E-4;}
            }
            break;
            
        case 2:
            if (z_bound == 0) {
                if (star_age < 272) {cum_return = 3.55E-7*logt - 6.64E-7;}
                else if (star_age < 890) {cum_return = 1.03E-4*(logt - 2.44) + 2.00E-7;}
                else {cum_return = 8.64E-7*(logt - 2.95) + 5.31E-5;}
            }
            else if (z_bound == 1) {
                if (star_age < 272) {cum_return = 3.73E-7*logt - 6.33E-7;}
                else if (star_age < 1544) {cum_return = 2.93E-5*(logt - 2.43) +  2.75E-7;}
                else {cum_return = 5.82E-7*(logt - 3.19) + 2.24E-5;}
            }
            else if (z_bound == 2) {
                if (star_age < 235) {cum_return = 1.06E-7*logt - 1.74E-7;}
                else if (star_age < 4812) {cum_return = 1.04E-6*(logt - 2.37) + 7.84E-8;}
                else {cum_return = 2.47E-7*(logt - 3.68) + 1.44E-6;}
            }
            else if (z_bound == 3) {
                if (star_age < 202) {cum_return = 6.38E-9*logt - 1.02E-8;}
                else if (star_age < 3394) {cum_return = 2.48E-8*(logt - 2.31) + 4.52E-9;}
                else {cum_return = 9.51E-9*(logt - 3.53) + 3.48E-8;}
            }
            else if (z_bound == 4) {
                if (star_age < 245) {cum_return = 2.97E-11*logt - 4.81E-11;}
                else if (star_age < 2392) {cum_return = 1.11E-10*(logt - 2.39) + 2.28E-11;}
                else {cum_return = 3.73E-13*(logt - 3.38) + 1.33E-10;}
            }
            break;
            
        case 3:
            if (z_bound == 0) {
                if (star_age < 525) {cum_return = 5.98E-6*logt - 9.97E-6;}
                else if (star_age < 1108) {cum_return = 1.02E-4*(logt - 2.72) + 6.30E-6;}
                else {cum_return = 5.98E-5*(logt - 3.04) + 3.94E-5;}
            }
            else if (z_bound == 1) {
                if (star_age < 339) {cum_return = 2.16E-6*logt - 3.60E-6;}
                else if (star_age < 1074) {cum_return = 4.09E-7*(logt - 2.53) + 1.87E-6;}
                else {cum_return = 1.50E-5*(logt - 3.03) + 2.08E-6;}
            }
            else if (z_bound == 2) {
                if (star_age < 307) {cum_return = 9.86E-7*logt - 1.62E-6;}
                else if (star_age < 4120) {cum_return = 1.13e-9*(logt - 2.49) + 8.34E-7;}
                else {cum_return = 2.23E-7*(logt - 3.61) + 8.36E-07;}
            }
            else if (z_bound == 3) {
                if (star_age < 253) {cum_return = 5.53E-9*logt - 8.71E-9;}
                else if (star_age < 297) {cum_return = 2.50E-9*(logt - 2.40) +  4.57E-9;}
                else {cum_return = 9.50E-12 *(logt - 2.47) + 4.75E-9;}
            }
            else if (z_bound == 4) {
                if (star_age < 128) {cum_return = 3.18E-11*logt - 5.00E-11;}
                else if (star_age < 269) {cum_return = 2.69E-11*(logt - 2.11) + 1.70E-11;}
                else {cum_return = 2.14E-14*(logt - 2.43) + 2.57E-11;}
            }
            break;
            
        default:
            cum_return = 0.;
            break;
    }
    if (cum_return < 0.) {cum_return = 0.;} // catch case were some fits are negative near the beginning of the AGB start time
    return cum_return;
}


/* Simple fit to cumulative AGB dust mass returns for a stellar population*/
double cumulative_AGB_dust_returns(int dust_type, double star_age, double z)
{
    /* dust_type: 0 = silicate, 1 = carbon, 2 = silicon carbide, 3 = metallic iron */
    double cumulative_mass;
    if (z <= 0.05)     {cumulative_mass = specific_Z_AGB_dust(dust_type,star_age,4);}
    else if (z <= 0.2) {cumulative_mass = specific_Z_AGB_dust(dust_type,star_age,4) + (z-0.05)/(0.2-0.05)*(specific_Z_AGB_dust(dust_type,star_age,3) - specific_Z_AGB_dust(dust_type,star_age,4));}
    else if (z <= 0.4) {cumulative_mass = specific_Z_AGB_dust(dust_type,star_age,3) + (z-0.2)/(0.4-0.2)*(specific_Z_AGB_dust(dust_type,star_age,2) - specific_Z_AGB_dust(dust_type,star_age,3));}
    else if (z <= 1.)  {cumulative_mass = specific_Z_AGB_dust(dust_type,star_age,2) + (z-0.4)/(1.-0.4)*(specific_Z_AGB_dust(dust_type,star_age,1) - specific_Z_AGB_dust(dust_type,star_age,2));}
    else if (z <= 2.)  {cumulative_mass = specific_Z_AGB_dust(dust_type,star_age,1) + (z-1.)/(2.-1.)*(specific_Z_AGB_dust(dust_type,star_age,0) - specific_Z_AGB_dust(dust_type,star_age,1));}
    else               {cumulative_mass = specific_Z_AGB_dust(dust_type,star_age,0);}
    return cumulative_mass;
}


/* initialize values of tables and variables for startup of runs */
void Initialize_ISMDustChem_Variables(int i)
{
    int j;
    /* atomic mass for each element in metallicity field, and some other variables. these always need to be initialized */
    All.ISMDustChem_AtomicMassTable[0] = 1.01;    // H
    All.ISMDustChem_AtomicMassTable[1] = 4.0;     // He
    All.ISMDustChem_AtomicMassTable[2] = 12.01;   // C
    All.ISMDustChem_AtomicMassTable[3] = 14;      // N
    All.ISMDustChem_AtomicMassTable[4] = 15.99;   // O
    All.ISMDustChem_AtomicMassTable[5] = 20.2;    // Ne
    All.ISMDustChem_AtomicMassTable[6] = 24.305;  // Mg
    All.ISMDustChem_AtomicMassTable[7] = 28.086;  // Si
    All.ISMDustChem_AtomicMassTable[8] = 32.065;  // S
    All.ISMDustChem_AtomicMassTable[9] = 40.078;  // Ca
    All.ISMDustChem_AtomicMassTable[10] = 55.845; // Fe
    All.ISMDustChem_SNeSputteringShutOffTime = 0.3E-3; // Destruction of dust due to SNe thermal sputtering ends around 0.3 Myr after SNe (from idealized SNe in Hu+2019)
    // Fiducial olivine-pyroxene silicate dust composition with olivine fraction = 0.63 and Mg frac = 0.65. If using iron nanoparticles assume iron is always present for silicate structure in the form of iron inclusions. index in metallicity field for elements which make up silicate dust (O,Mg,Si)
    All.ISMDustChem_SilicateMetallicityFieldIndexTable[0] = 4;
    All.ISMDustChem_SilicateMetallicityFieldIndexTable[1] = 6;
    All.ISMDustChem_SilicateMetallicityFieldIndexTable[2] = 7;
    // number of O, Mg, and Si in one formula unit of silicate dust
    All.ISMDustChem_SilicateNumberOfAtomsTable[0] = 3.63;
    All.ISMDustChem_SilicateNumberOfAtomsTable[1] = 1.06;
    All.ISMDustChem_SilicateNumberOfAtomsTable[2] = 1.;
    if(!(GALSF_ISMDUSTCHEM_MODEL & 4)) {All.ISMDustChem_SilicateMetallicityFieldIndexTable[3] = 10; All.ISMDustChem_SilicateNumberOfAtomsTable[3] = 0.571;} // add Fe as well if not accounting for iron inclusions
    All.ISMDustChem_EffectiveSilicateDustAtomicWeight = 0.; for(j=0;j<GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES;j++) {All.ISMDustChem_EffectiveSilicateDustAtomicWeight += All.ISMDustChem_SilicateNumberOfAtomsTable[j] * All.ISMDustChem_AtomicMassTable[All.ISMDustChem_SilicateMetallicityFieldIndexTable[j]];}

    /* only initialize these on a new run */
    if(RestartFlag == 0) {
        SphP[i].ISMDustChem_DelayTimeSNeSputtering = SphP[i].ISMDustChem_C_in_CO = SphP[i].ISMDustChem_MassFractionInDenseMolecular = 0.;
        if(All.Initial_ISMDustChem_Depletion > 0)
        {
            for(j=0;j<NUM_ISMDUSTCHEM_ELEMENTS;j++) {SphP[i].ISMDustChem_Dust_Metal[j] = 0.;}
            if(GALSF_ISMDUSTCHEM_MODEL & 1) {
                SphP[i].ISMDustChem_Dust_Metal[4] = All.Initial_ISMDustChem_Depletion*P[i].Metallicity[4]; // Silicate dust O
                SphP[i].ISMDustChem_Dust_Metal[6] = All.Initial_ISMDustChem_Depletion*P[i].Metallicity[6]; // Silicate dust Mg
                SphP[i].ISMDustChem_Dust_Metal[7] = All.Initial_ISMDustChem_Depletion*P[i].Metallicity[7]; // Silicate dust Si
                SphP[i].ISMDustChem_Dust_Metal[10] = All.Initial_ISMDustChem_Depletion*P[i].Metallicity[10]; // Silicate dust Fe
                SphP[i].ISMDustChem_Dust_Metal[2] = DMIN(P[i].Metallicity[2],SphP[i].ISMDustChem_Dust_Metal[4]+SphP[i].ISMDustChem_Dust_Metal[6]+SphP[i].ISMDustChem_Dust_Metal[7]+SphP[i].ISMDustChem_Dust_Metal[10]/All.Initial_ISMDustChem_SiliconToCarbonRatio); // Carbonaceous dust
            }
            if(GALSF_ISMDUSTCHEM_MODEL & 2) {
                // Silicate dust
                double sil_mass_frac=0.; SphP[i].ISMDustChem_Dust_Metal[7] = All.Initial_ISMDustChem_Depletion*P[i].Metallicity[7]; // Set Si depletion
                sil_mass_frac+=SphP[i].ISMDustChem_Dust_Metal[7];
                for(j=0;j<GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES;j++) // Set element depletions for all other elements in silicates given initial Si depletion
                {
                    if(j != 2)
                    {
                        SphP[i].ISMDustChem_Dust_Metal[All.ISMDustChem_SilicateMetallicityFieldIndexTable[j]] += SphP[i].ISMDustChem_Dust_Metal[7] / (All.ISMDustChem_SilicateNumberOfAtomsTable[2] * All.ISMDustChem_AtomicMassTable[7]) * (All.ISMDustChem_SilicateNumberOfAtomsTable[j] * All.ISMDustChem_AtomicMassTable[All.ISMDustChem_SilicateMetallicityFieldIndexTable[j]]);
                        sil_mass_frac += SphP[i].ISMDustChem_Dust_Metal[All.ISMDustChem_SilicateMetallicityFieldIndexTable[j]];
                    }
                }
                SphP[i].ISMDustChem_Dust_Species[0] = sil_mass_frac;
                // Carbonaceous dust
                SphP[i].ISMDustChem_Dust_Metal[2] = DMIN(P[i].Metallicity[2],sil_mass_frac/All.Initial_ISMDustChem_SiliconToCarbonRatio); SphP[i].ISMDustChem_Dust_Species[1] = SphP[i].ISMDustChem_Dust_Metal[2];
                if(GALSF_ISMDUSTCHEM_MODEL & 4) { // Metallic Iron Nanoparticles
                    SphP[i].ISMDustChem_Dust_Metal[10] = All.Initial_ISMDustChem_Depletion*P[i].Metallicity[10];
                    SphP[i].ISMDustChem_Dust_Species[3] = (1.-GALSF_ISMDUSTCHEM_VAR_IRON_INCL_FRAC)*SphP[i].ISMDustChem_Dust_Metal[10];
                    SphP[i].ISMDustChem_Dust_Species[5] = GALSF_ISMDUSTCHEM_VAR_IRON_INCL_FRAC*SphP[i].ISMDustChem_Dust_Metal[10];
                }
            }
            for (j=1;j<NUM_ISMDUSTCHEM_ELEMENTS;j++) {SphP[i].ISMDustChem_Dust_Metal[0] += SphP[i].ISMDustChem_Dust_Metal[j];}
            for (j=0;j<NUM_ISMDUSTCHEM_SOURCES;j++) {SphP[i].ISMDustChem_Dust_Source[j] = 0.;}
            SphP[i].ISMDustChem_Dust_Source[2] = SphP[i].ISMDustChem_Dust_Metal[0];  // Assume initial dust population is from SNe II
        }
        else
        {
            for (j=0;j<NUM_ISMDUSTCHEM_ELEMENTS;j++) {SphP[i].ISMDustChem_Dust_Metal[j] = 0.;}
            for (j=0;j<NUM_ISMDUSTCHEM_SOURCES;j++) {SphP[i].ISMDustChem_Dust_Source[j] = 0.;}
            for (j=0;j<NUM_ISMDUSTCHEM_SPECIES;j++) {SphP[i].ISMDustChem_Dust_Species[j] = 0.;}
        }
    }
}


/* simple indexing routine to return the value we need when looping over yields and the like */
double return_ismdustchem_species_of_interest_for_diffusion_and_yields(int i, int k)
{
    k -= NUM_METAL_SPECIES;
    if(k<NUM_ISMDUSTCHEM_ELEMENTS) {return SphP[i].ISMDustChem_Dust_Metal[k];}
    k -= NUM_ISMDUSTCHEM_ELEMENTS;
    if(k<NUM_ISMDUSTCHEM_SOURCES) {return SphP[i].ISMDustChem_Dust_Source[k];}
    k -= NUM_ISMDUSTCHEM_SOURCES;
    if(k<NUM_ISMDUSTCHEM_SPECIES) {return SphP[i].ISMDustChem_Dust_Species[k];}
    return 0;
}


/* Approximate dust cooling via electron-dust collisions for MRN sized dust in plasmas from Dwek(1987)+Dewk&Werner(1981). Should surpass metal-line cooling for >10^6 K, but this will also overpredict dust cooling for <10^7 K since cooling is dominated by small grains which should be destroyed via sputtering */
double Lambda_Dust_HighTemperature_Gas_ISM(int target, double T, double n_elec)
{
    if(target<0 || T<1.e5) {return 0;}
    if(SphP[target].ISMDustChem_Dust_Metal[0] <= 0) {return 0;}
    double rho_c=3., a3=2.21e-18, DG, Havg; // rho_c=gm cm^-3 grain solid density (intermediate between silicate and carbonaceous), a3=cm^3 average grain volume for MRN grain size distribution
    if (T>=7.17E7) {Havg=1.43E-11;} else if (T>=2.4E7) {Havg=-2.07E-12+1.23E-16*pow(T,0.745)+2.10E-17*pow(T,0.75)-1.07E-17*pow(T,0.88);}
    else if (T>=4.55E6) {Havg=-2.07E-12+1.70E-17*pow(T,0.745)+3.96E-17*pow(T,0.75)-5.44E-23*pow(T,1.5);}
    else if (T>=1.52E6) {Havg=-1.06E-16*pow(T,0.745)+1.86E-17*pow(T,0.75)+1.56E-17*pow(T,0.88)-5.44E-23*pow(T,1.5);} else {Havg=3.76E-22*pow(T,1.5);}
    double coolrate = (3./(4.*M_PI)*(PROTONMASS_CGS/HYDROGEN_MASSFRAC))*(SphP[target].ISMDustChem_Dust_Metal[0]/rho_c*n_elec*(Havg/a3));
    if(!isfinite(coolrate)) {coolrate=0;}
    return coolrate;
}



/* return the mass fraction we will assume of dust destroyed in surrounding gas due to SNe shock, taken from McKee 1989 and Cioffi 1988 */
double ISMDustChem_Return_Mass_Fraction_Where_Dust_Destroyed(double rho_cell_in_code_units, double Esne51_into_cell, double mass_preshock_in_code_units)
{
    double dest_eff=0.4, vs7=1., local_n0=rho_cell_in_code_units*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS; // dust destruction efficiency, minimum gas shock velocity in 10^7 cm/s which destroys dust, and number density around SNe
    double mass_shocked = dest_eff * 2460 * Esne51_into_cell / (pow(local_n0, 0.1) * pow(vs7, 9./7.) * UNIT_MASS_IN_SOLAR); // mass shocked to 100 km/s which destroys dust. use the weights to distribute shocked mass across the neighboring gas particles
    return DMIN(1., mass_shocked / mass_preshock_in_code_units); // mass fraction destroyed
}


/* Subroutine to update the dust abundances after enrichment in the mechanical feedback subroutine */
void update_ISMDustChem_after_mechanical_injection(int j, double massfrac_destroyed, double m0, double mf, double *Z_injected)
{
    // If SNe events happened need to first destroy the appropriate amount of dust if there is any dust
    int k;
    if((massfrac_destroyed > 0) && (SphP[j].ISMDustChem_Dust_Metal[0] > 0))
    {
        SphP[j].ISMDustChem_DelayTimeSNeSputtering = All.ISMDustChem_SNeSputteringShutOffTime; // update thermal sputtering delay time due to SNe
        if (massfrac_destroyed >= 1.) // destroy all dust
        {
            for(k=0;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {SphP[j].ISMDustChem_Dust_Metal[k]=0.;}
            for(k=0;k<NUM_ISMDUSTCHEM_SOURCES;k++)  {SphP[j].ISMDustChem_Dust_Source[k]=0.;}
            for(k=0;k<NUM_ISMDUSTCHEM_SPECIES;k++)  {SphP[j].ISMDustChem_Dust_Species[k]=0.;}
        }
        else
        {
            double protected_frac = 0.; // Fraction of dust protected from destruction (only iron inclusions are currently considered)
            if(GALSF_ISMDUSTCHEM_MODEL & 4) { // Take out the iron inclusions protected in silicate dust and then add it back in later
                protected_frac = SphP[j].ISMDustChem_Dust_Species[NUM_ISMDUSTCHEM_SPECIES-1]/SphP[j].ISMDustChem_Dust_Metal[0];
                SphP[j].ISMDustChem_Dust_Metal[10] -= SphP[j].ISMDustChem_Dust_Species[NUM_ISMDUSTCHEM_SPECIES-1]; // Assume all dust species are destroyed evenly but leave out iron inclusions
            }
            for(k=0;k<NUM_ISMDUSTCHEM_SPECIES-1;k++) {SphP[j].ISMDustChem_Dust_Species[k] *= 1.-massfrac_destroyed;} // Assume all dust species are destroyed evenly
            // Assume all dust sources are destroyed evenly and take into account protected dust
            for(k=0;k<NUM_ISMDUSTCHEM_SOURCES;k++) {SphP[j].ISMDustChem_Dust_Source[k] *= (1.-(1.-protected_frac)*massfrac_destroyed);}
            SphP[j].ISMDustChem_Dust_Metal[0] = 0.0;
            for(k=1;k<NUM_ISMDUSTCHEM_ELEMENTS;k++)
            {
                SphP[j].ISMDustChem_Dust_Metal[k] *= 1.-massfrac_destroyed;
                SphP[j].ISMDustChem_Dust_Metal[0] += SphP[j].ISMDustChem_Dust_Metal[k];
            }
            if(GALSF_ISMDUSTCHEM_MODEL & 4) { // Add the protected iron dust back in
                SphP[j].ISMDustChem_Dust_Metal[10] += SphP[j].ISMDustChem_Dust_Species[NUM_ISMDUSTCHEM_SPECIES-1];
                SphP[j].ISMDustChem_Dust_Metal[0] += SphP[j].ISMDustChem_Dust_Species[NUM_ISMDUSTCHEM_SPECIES-1];
                // Update amount of free-flying iron and iron inclusions since some of the inclusions are released from silicate. Assume this leads to constant fraction of iron inclusions that scales with amount of silicate dust
                int key_elem = 0; double key_mass, key_num_atoms, frac_of_max_sil, incl_frac, sil_elem_abunds[GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES];
                for(k=0;k<GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES;k++)
                {
                    int index = All.ISMDustChem_SilicateMetallicityFieldIndexTable[k];
                    sil_elem_abunds[k] = P[j].Metallicity[index] / (All.ISMDustChem_AtomicMassTable[index] * PROTONMASS_CGS);
                    if (sil_elem_abunds[key_elem] / All.ISMDustChem_SilicateNumberOfAtomsTable[key_elem] > sil_elem_abunds[k] / All.ISMDustChem_SilicateNumberOfAtomsTable[k]) key_elem = k;
                }
                key_mass = All.ISMDustChem_AtomicMassTable[All.ISMDustChem_SilicateMetallicityFieldIndexTable[key_elem]];
                key_num_atoms = All.ISMDustChem_SilicateNumberOfAtomsTable[key_elem];
                key_elem = All.ISMDustChem_SilicateMetallicityFieldIndexTable[key_elem];
                frac_of_max_sil = SphP[j].ISMDustChem_Dust_Species[0] / (P[j].Metallicity[key_elem] * All.ISMDustChem_EffectiveSilicateDustAtomicWeight/(key_num_atoms * key_mass));
                incl_frac = DMAX(DMIN(GALSF_ISMDUSTCHEM_VAR_IRON_INCL_FRAC*frac_of_max_sil,GALSF_ISMDUSTCHEM_VAR_IRON_INCL_FRAC),0.);
                SphP[j].ISMDustChem_Dust_Species[3] = (1.-incl_frac) * SphP[j].ISMDustChem_Dust_Metal[10];
                SphP[j].ISMDustChem_Dust_Species[NUM_ISMDUSTCHEM_SPECIES-1] = incl_frac * SphP[j].ISMDustChem_Dust_Metal[10];
            }
        }
    }
    for(k=0;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {SphP[j].ISMDustChem_Dust_Metal[k]   = (m0/mf)*SphP[j].ISMDustChem_Dust_Metal[k]   + (1./mf)*DMAX(0.,Z_injected[k+NUM_METAL_SPECIES]);}
    for(k=0;k<NUM_ISMDUSTCHEM_SOURCES;k++)  {SphP[j].ISMDustChem_Dust_Source[k]  = (m0/mf)*SphP[j].ISMDustChem_Dust_Source[k]  + (1./mf)*DMAX(0.,Z_injected[k+NUM_METAL_SPECIES+NUM_ISMDUSTCHEM_ELEMENTS]);}
    for(k=0;k<NUM_ISMDUSTCHEM_SPECIES;k++)  {SphP[j].ISMDustChem_Dust_Species[k] = (m0/mf)*SphP[j].ISMDustChem_Dust_Species[k] + (1./mf)*DMAX(0.,Z_injected[k+NUM_METAL_SPECIES+NUM_ISMDUSTCHEM_ELEMENTS+NUM_ISMDUSTCHEM_SOURCES]);}
}



/* subroutine to update dust masses from growth via gas-dust accretion and destruction via thermal sputtering */
void update_dust_acc_and_sput(int i, double dtime_gyr)
{
    int k; double ne=1, nh0=0, nHe0, nHepp, nhp, nHeII, temp, mu_meanwt=1, rho=SphP[i].Density*All.cf_a3inv, u0=SphP[i].InternalEnergyPred;
    temp = ThermalProperties(u0, rho, i, &mu_meanwt, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp);
    rho*=UNIT_DENSITY_IN_CGS;
    
    /* Need to calculate H2 fraction to determine whether gas is either in the CNM/diffuse MC or dense MC phase since
     * gas-dust accretion has Coloumb enhancing in CNM/diffuse MC due to dust grain charge and ionized metal species.
     * Also used to determine when C is locked up in CO which need to be taken into account since CO rapidly forms
     * once the gas is sufficently molecular.
     */
    double fH2=0., new_ISMDustChem_MassFractionInDenseMolecular=0.; // mass fraction of gas that is H2 and gas in dense MC phase
    double NH2 = 1.5E21; // cm^-2 Column density of H2 needed to be in dense MC phase
    double l_depth, x_dens; // depth into cloud to reach NH2 and radial fraction of cloud in dense MC phase
    double surface_density = evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1,i) * UNIT_SURFDEN_IN_CGS; // converts to cgs
    // shielding length giving effective radius of gas particle
    double l_shield = surface_density / rho;
    fH2 = Get_Gas_Molecular_Mass_Fraction(i, temp, nh0, ne, 0.);
    if (fH2 > 0)
    {
        double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS_CGS;
        l_depth = 2.*NH2 / (fH2*nHcgs);
        x_dens = (l_shield-l_depth)/l_shield;
        x_dens = DMIN(1,DMAX(0,x_dens));
        new_ISMDustChem_MassFractionInDenseMolecular = pow(x_dens,3);
        new_ISMDustChem_MassFractionInDenseMolecular = DMIN(new_ISMDustChem_MassFractionInDenseMolecular,fH2); // Maximum dense molecular fraction set by total molecular fraction
    }
    // Only need to update CO if there is C present
    if (P[i].Metallicity[2]>0)
    {
        // If dense MC has shrunk, reduce the C in CO by the fraction it has shrunk
        if (new_ISMDustChem_MassFractionInDenseMolecular < SphP[i].ISMDustChem_MassFractionInDenseMolecular) {SphP[i].ISMDustChem_C_in_CO *= new_ISMDustChem_MassFractionInDenseMolecular/SphP[i].ISMDustChem_MassFractionInDenseMolecular;}
        // If dense MC has grown, increase the C in CO by the newly add volume of remaining gas-phase C if any is left
        else
        {
            if (P[i].Metallicity[2]-SphP[i].ISMDustChem_Dust_Metal[2]-SphP[i].ISMDustChem_C_in_CO > 0.)
            {
                SphP[i].ISMDustChem_C_in_CO += (new_ISMDustChem_MassFractionInDenseMolecular-SphP[i].ISMDustChem_MassFractionInDenseMolecular) * ((P[i].Metallicity[2]-SphP[i].ISMDustChem_Dust_Metal[2])-SphP[i].ISMDustChem_C_in_CO) / (1.-SphP[i].ISMDustChem_MassFractionInDenseMolecular);
            }
        }
    }
    else
    {
        SphP[i].ISMDustChem_C_in_CO = 0.;
    }
    SphP[i].ISMDustChem_MassFractionInDenseMolecular = new_ISMDustChem_MassFractionInDenseMolecular;
    
    if (SphP[i].ISMDustChem_Dust_Metal[0] <= 0) {return;} // No dust so nothing more to do
    
    double dF; // change in fraction of element condensed into dust
    double growth_timescale, sputter_timescale, t_ref, T_ref, avg_grain_radius;
    double dust_yields[NUM_ISMDUSTCHEM_ELEMENTS] = {0.0};
    int source = 0;
    
    // now accrete and sputter dust //
#if (GALSF_ISMDUSTCHEM_MODEL & 1)
    SphP[i].ISMDustChem_Dust_Metal[0] = 0.; // First renorm dust due to building numerical error that can arise from stellar feedback. This may no longer be necessary.
    for (k=2;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {SphP[i].ISMDustChem_Dust_Metal[0] += SphP[i].ISMDustChem_Dust_Metal[k];}
    double total = SphP[i].ISMDustChem_Dust_Source[0]+SphP[i].ISMDustChem_Dust_Source[1]+SphP[i].ISMDustChem_Dust_Source[2]+SphP[i].ISMDustChem_Dust_Source[3];
    for (k=0;k<NUM_ISMDUSTCHEM_SOURCES;k++) SphP[i].ISMDustChem_Dust_Source[k] = DMAX(0,SphP[i].ISMDustChem_Dust_Metal[0]/total*SphP[i].ISMDustChem_Dust_Source[k]);
    
    
    /* First the newly created dust mass due to accretion is added with the creation source updated, then the dust destroyed
     due to thermal sputtering is subtracted from the overall dust mass with the creation source being unchanged since we
     assume the dust is destroyed evenly regardless of creation source */
    double rho_ref = PROTONMASS_CGS; // 1 H atom cm^-3
    T_ref = 20.; avg_grain_radius = 0.032; /* um */ t_ref = 0.2; /* Gyr */
    growth_timescale = t_ref * (rho_ref / rho) * pow((T_ref / temp), .5);
    // Calculate the fraction of mass of a certain element to be added to dust due to accretion
    for (k=2;k<NUM_ISMDUSTCHEM_ELEMENTS;k++)
    {
        double in_mol_frac; // fraction of element in molecular form and unable to accrete onto dust (CO is the only molecule considered)
        if (k==2) {in_mol_frac = SphP[i].ISMDustChem_C_in_CO;}
        else if (k==4) {in_mol_frac = SphP[i].ISMDustChem_C_in_CO * All.ISMDustChem_AtomicMassTable[4] / All.ISMDustChem_AtomicMassTable[2];}
        else {in_mol_frac = 0.;}
        // If no dust, metals, or all metals in dust then no accretion
        if (P[i].Metallicity[k] == 0. || SphP[i].ISMDustChem_Dust_Metal[k] == 0. || (P[i].Metallicity[k] - SphP[i].ISMDustChem_Dust_Metal[k]) <= 0) {dF = 0.;}
        else
        {
            dF = dtime_gyr * (1. - SphP[i].ISMDustChem_Dust_Metal[k] / (P[i].Metallicity[k] - in_mol_frac)) * (SphP[i].ISMDustChem_Dust_Metal[k] / growth_timescale);
            // Check in case we use up the rest of the remaining metal in the gas phase and deal with unphysical values
            dF = DMIN(P[i].Metallicity[k] - SphP[i].ISMDustChem_Dust_Metal[k] - in_mol_frac,DMAX(0.,dF));
            dust_yields[k] = dF;
            dust_yields[0] += dust_yields[k];
        }
    }
    // Update dust yields and creation source
    if (dust_yields[0] != 0.)
    {
        SphP[i].ISMDustChem_Dust_Source[source] += dust_yields[0];
        for (k=0;k< NUM_ISMDUSTCHEM_ELEMENTS;k++) {SphP[i].ISMDustChem_Dust_Metal[k] += dust_yields[k];}
    }
    
    // Check if sputtering is delayed due to recent SNe
    if(SphP[i].ISMDustChem_DelayTimeSNeSputtering > 0) {SphP[i].ISMDustChem_DelayTimeSNeSputtering = DMAX(0,SphP[i].ISMDustChem_DelayTimeSNeSputtering-dtime_gyr);} // count off clock since last SNe
    else // Now determine amount of dust destroyed by thermal sputtering
    {
        T_ref = 2E6; avg_grain_radius = 0.032; /* um */ t_ref = 0.17; /* Gyr */
        sputter_timescale = t_ref * (avg_grain_radius / 0.1) / (rho*1E27) * (pow((T_ref/ temp), 2.5) + 1.);
        for (k=0;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {dust_yields[k] = 0.;}
        // Calculate the fraction of mass of a certain element to be destroyed due to thermal sputtering
        for (k=2;k<NUM_ISMDUSTCHEM_ELEMENTS;k++)
        {
            if (SphP[i].ISMDustChem_Dust_Metal[k] <= 0.) {dF = 0.;}
            else {dF = - dtime_gyr * (SphP[i].ISMDustChem_Dust_Metal[k] / (sputter_timescale / 3.));}
            // can't destroy more dust then there is available and deal with unphysical values
            dF = DMAX(-SphP[i].ISMDustChem_Dust_Metal[k],DMIN(0,dF));
            dust_yields[k] = dF;
            dust_yields[0] += dF;
        }
        
        // Update dust yields and sources
        if (dust_yields[0] != 0.)
        {
            // Assume all dust sources are destroyed evenly
            for(k=0;k<NUM_ISMDUSTCHEM_SOURCES;k++) {SphP[i].ISMDustChem_Dust_Source[k] *= (1.+dust_yields[0]/SphP[i].ISMDustChem_Dust_Metal[0]);}
            for (k=0;k<NUM_ISMDUSTCHEM_ELEMENTS;k++)
            {
                SphP[i].ISMDustChem_Dust_Metal[k] += dust_yields[k];
            }
            // Deal with rounding error causing total dust to not equal zero
            int no_dust = 1;
            for (k=2; k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {if (SphP[i].ISMDustChem_Dust_Metal[k] > 0.) {no_dust = 0; break;}}
            if (no_dust)
            {
                SphP[i].ISMDustChem_Dust_Metal[0] = 0.;
                // if all dust is destroyed need to zero creation sources
                for(k=0;k<NUM_ISMDUSTCHEM_SOURCES;k++) {SphP[i].ISMDustChem_Dust_Source[k] = 0.;}
            }
        }
    }
#endif // model == 1, elemental model
#if (GALSF_ISMDUSTCHEM_MODEL & 2)
    // First renorm dust due to building numerical error that can arise from stellar feedback. This may no longer be necessary.
    for (k=0;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {SphP[i].ISMDustChem_Dust_Metal[k]=0.;}
    // silicate
    for (k=0;k<GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES;k++)
    {
        SphP[i].ISMDustChem_Dust_Metal[All.ISMDustChem_SilicateMetallicityFieldIndexTable[k]] += SphP[i].ISMDustChem_Dust_Species[0] * All.ISMDustChem_SilicateNumberOfAtomsTable[k] * All.ISMDustChem_AtomicMassTable[All.ISMDustChem_SilicateMetallicityFieldIndexTable[k]] / All.ISMDustChem_EffectiveSilicateDustAtomicWeight;
    }
    // carbonaceous
    SphP[i].ISMDustChem_Dust_Metal[2] += SphP[i].ISMDustChem_Dust_Species[1];
    // SiC
    SphP[i].ISMDustChem_Dust_Metal[2] += SphP[i].ISMDustChem_Dust_Species[2] * All.ISMDustChem_AtomicMassTable[2] / (All.ISMDustChem_AtomicMassTable[2] + All.ISMDustChem_AtomicMassTable[7]);
    SphP[i].ISMDustChem_Dust_Metal[7] += SphP[i].ISMDustChem_Dust_Species[2] * All.ISMDustChem_AtomicMassTable[7] / (All.ISMDustChem_AtomicMassTable[2] + All.ISMDustChem_AtomicMassTable[7]);
    if(GALSF_ISMDUSTCHEM_MODEL & 4) {SphP[i].ISMDustChem_Dust_Metal[10] += SphP[i].ISMDustChem_Dust_Species[3]+SphP[i].ISMDustChem_Dust_Species[5];} else {SphP[i].ISMDustChem_Dust_Metal[10] += SphP[i].ISMDustChem_Dust_Species[3];} // metallic iron
    if(GALSF_ISMDUSTCHEM_MODEL & 8) {SphP[i].ISMDustChem_Dust_Metal[4] += SphP[i].ISMDustChem_Dust_Species[4];} // oxygen reservoir
    for (k=2;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {SphP[i].ISMDustChem_Dust_Metal[0] += SphP[i].ISMDustChem_Dust_Metal[k];}
    double total = SphP[i].ISMDustChem_Dust_Source[0]+SphP[i].ISMDustChem_Dust_Source[1]+SphP[i].ISMDustChem_Dust_Source[2]+SphP[i].ISMDustChem_Dust_Source[3];
    if (total<=0.) {for (k=0;k<NUM_ISMDUSTCHEM_SOURCES;k++) SphP[i].ISMDustChem_Dust_Source[k] = 0.;}
    else {for (k=0;k<NUM_ISMDUSTCHEM_SOURCES;k++) SphP[i].ISMDustChem_Dust_Source[k] = DMAX(0,SphP[i].ISMDustChem_Dust_Metal[0]/total*SphP[i].ISMDustChem_Dust_Source[k]);}
    
    /* Restrict accretion only to MC environments by assuming sticking efficiency of 1 for
     * T <= 300K and 0 otherwise. Check the three main dust species that can form through
     * accretion in the ISM silicates, carbon, and metallic iron.
     */
    double bulk_dens; //  mass density of the dust condensed phase (g cm^-3)
    double dust_formula_mass; // atomic weight of one formula unit of dust species
    double max_num_dens; // max number density of key element assuming all of element is in gas phase
    double species_yields[NUM_ISMDUSTCHEM_SPECIES] = {0.0};
    int key_elem;  // index of least abundant element needed to create the dust species under consideration
    double key_mass, key_num_atoms; // atomic mass of key element number of atoms in dust species
    double num_dens[NUM_ISMDUSTCHEM_ELEMENTS]; // number density of all elements in gas only
    int missing_element; // check if any elements are missing from the gas phase
    // reference accretion timescales for ionized (with Coulomb enhancment) and neutral (no enhancement) gas-phase metals
    double t_ref_CNM, t_ref_MC;
    // Assuming sticking efficiency of zero for T > 300 K
    if (temp <= 300)
    {
        // Determine number density of each element in the gas phase, use this to determine the key element for each dust species
        num_dens[0] = rho * (1. - P[i].Metallicity[0] - P[i].Metallicity[1]) / (All.ISMDustChem_AtomicMassTable[0] * PROTONMASS_CGS);
        for (k=1;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) num_dens[k] = rho * (P[i].Metallicity[k] - SphP[i].ISMDustChem_Dust_Metal[k])/ (All.ISMDustChem_AtomicMassTable[k] * PROTONMASS_CGS);
        
        /******** SILICATE ********/
        t_ref_CNM = 0.252E-3;   // Gyr
        t_ref_MC = 1.38E-3;     // Gyr
        // Calculate effective timescale assuming produced dust is redistributed throughout the gas
        t_ref = (t_ref_CNM * t_ref_MC) / (SphP[i].ISMDustChem_MassFractionInDenseMolecular * t_ref_CNM + (1.-SphP[i].ISMDustChem_MassFractionInDenseMolecular) * t_ref_MC);
        // check that all the elements required to make silicates are available
        missing_element = 0;
        for(k=0;k<GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES;k++) {if(num_dens[All.ISMDustChem_SilicateMetallicityFieldIndexTable[k]] <= 0.) {missing_element = 1;}}
        if(!missing_element)
        {
            bulk_dens = 3.13; dust_formula_mass = All.ISMDustChem_EffectiveSilicateDustAtomicWeight; key_elem = 0;
            // determine key element which sets the accretion rate
            for (k=1;k<GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES;k++)
            {
                if (num_dens[All.ISMDustChem_SilicateMetallicityFieldIndexTable[key_elem]]/All.ISMDustChem_SilicateNumberOfAtomsTable[key_elem] > num_dens[All.ISMDustChem_SilicateMetallicityFieldIndexTable[k]]/All.ISMDustChem_SilicateNumberOfAtomsTable[k]) {key_elem = k;}
            }
            key_mass = All.ISMDustChem_AtomicMassTable[All.ISMDustChem_SilicateMetallicityFieldIndexTable[key_elem]];
            key_num_atoms = All.ISMDustChem_SilicateNumberOfAtomsTable[key_elem];
            key_elem = All.ISMDustChem_SilicateMetallicityFieldIndexTable[key_elem];
            max_num_dens = rho * P[i].Metallicity[key_elem] / (key_mass * PROTONMASS_CGS);
            
            growth_timescale = t_ref * (key_num_atoms * sqrt(key_mass) / dust_formula_mass) * bulk_dens / max_num_dens / sqrt(temp);
            // change in dust condensation for key element
            dF = dtime_gyr * (1. - SphP[i].ISMDustChem_Dust_Metal[key_elem] / P[i].Metallicity[key_elem]) * SphP[i].ISMDustChem_Dust_Metal[key_elem] / growth_timescale;
            // Check in case we use up the rest of the remaining metal in the gas phase and deal with unphysical values
            dF = DMAX(0.,DMIN(P[i].Metallicity[key_elem]-SphP[i].ISMDustChem_Dust_Metal[key_elem],dF));
            
            // change in dust condensation for all elements in silicates
            for (k=0;k<GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES;k++)
            {
                dust_yields[All.ISMDustChem_SilicateMetallicityFieldIndexTable[k]] += dF * ((All.ISMDustChem_SilicateNumberOfAtomsTable[k] * All.ISMDustChem_AtomicMassTable[All.ISMDustChem_SilicateMetallicityFieldIndexTable[k]])/ (key_num_atoms * key_mass));
            }
            species_yields[0] = dF * (dust_formula_mass / (key_num_atoms * key_mass));
        }
        
        /******** CARBONACEOUS ********/
        // Since the transition between C+ -> C -> CO is quick, assume C+ -> CO so carbon dust only grows in CNM environments
        // Also need to take into account C in CO reduces the maximum amount of carbon dust which can be formed
        if (SphP[i].ISMDustChem_MassFractionInDenseMolecular < 1.)
        {
            t_ref_CNM = 1.54E-3; // Gyr
            t_ref = t_ref_CNM / (1.-SphP[i].ISMDustChem_MassFractionInDenseMolecular);
            key_elem = 2; key_mass = All.ISMDustChem_AtomicMassTable[key_elem]; key_num_atoms = 1.; bulk_dens = 2.25; dust_formula_mass = key_mass;
            if (num_dens[key_elem] > 0)
            {
                max_num_dens = rho * P[i].Metallicity[key_elem]/ (key_mass * PROTONMASS_CGS);
                growth_timescale = t_ref * sqrt(key_mass) / dust_formula_mass * bulk_dens / max_num_dens / sqrt(temp);
                dF = dtime_gyr * (1. - SphP[i].ISMDustChem_Dust_Metal[key_elem] / (P[i].Metallicity[key_elem] - SphP[i].ISMDustChem_C_in_CO)) * SphP[i].ISMDustChem_Dust_Metal[key_elem] / growth_timescale;
                // Check in case we use up the rest of the remaining metal in the gas phase and deal with unphysical values
                dF = DMAX(0,DMIN(P[i].Metallicity[key_elem]-SphP[i].ISMDustChem_C_in_CO-SphP[i].ISMDustChem_Dust_Metal[key_elem],dF));
                dust_yields[key_elem] += dF;
                species_yields[1] = dF;
            }
        }
        
        /******** METALLIC IRON ********/
        if(GALSF_ISMDUSTCHEM_MODEL & 4) {t_ref_CNM = 1.66E-6; t_ref_MC = 0.139E-3;} else {t_ref_CNM = 0.252E-3; t_ref_MC = 1.38E-3;} // Gyr
        // Calculate effective timescale assuming produced dust is redistributed throughout the gas
        t_ref = (t_ref_CNM * t_ref_MC) / (SphP[i].ISMDustChem_MassFractionInDenseMolecular * t_ref_CNM + (1.-SphP[i].ISMDustChem_MassFractionInDenseMolecular) * t_ref_MC);
        key_elem = 10; key_mass = All.ISMDustChem_AtomicMassTable[key_elem]; key_num_atoms = 1.; bulk_dens = 7.86; dust_formula_mass = All.ISMDustChem_AtomicMassTable[key_elem];
        if (num_dens[key_elem] > 0)
        {
            max_num_dens = rho * P[i].Metallicity[key_elem]/(key_mass * PROTONMASS_CGS);
            growth_timescale = t_ref * sqrt(key_mass) / dust_formula_mass * bulk_dens / max_num_dens / sqrt(temp);
            dF = dtime_gyr * (1. - SphP[i].ISMDustChem_Dust_Metal[key_elem] / P[i].Metallicity[key_elem]) * SphP[i].ISMDustChem_Dust_Species[3] / growth_timescale;
            // Check in case we use up the rest of the remaining metal in the gas phase and deal with unphysical values
            dF = DMAX(0.,DMIN(P[i].Metallicity[key_elem]-SphP[i].ISMDustChem_Dust_Metal[key_elem],dF));
            dust_yields[key_elem] += dF;
            species_yields[3] = dF;
        }
        
        // Add up all the dust elements
        for (k=2;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) dust_yields[0] += dust_yields[k];
        
        // Update dust yields and creation source
        if (dust_yields[0] != 0.)
        {
            // update dust source
            SphP[i].ISMDustChem_Dust_Source[source] += dust_yields[0];
            for (k=0;k<NUM_ISMDUSTCHEM_SPECIES;k++) SphP[i].ISMDustChem_Dust_Species[k] += species_yields[k];
            // update new dust mass
            for (k = 0; k < NUM_ISMDUSTCHEM_ELEMENTS; k++) {SphP[i].ISMDustChem_Dust_Metal[k] += dust_yields[k];}
            if(GALSF_ISMDUSTCHEM_MODEL & 4) { // Update amount of free-flying iron and iron inclusions since some of the free-flying particles become inclusions in silicates. Scales with local amount of silicates
                int key_elem = 0; double key_mass, key_num_atoms, frac_of_max_sil, incl_frac, sil_elem_abunds[GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES];
                for(k=0;k<GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES;k++)
                {
                    int index = All.ISMDustChem_SilicateMetallicityFieldIndexTable[k]; sil_elem_abunds[k] = P[i].Metallicity[index] / (All.ISMDustChem_AtomicMassTable[index] * PROTONMASS_CGS);
                    if (sil_elem_abunds[key_elem] / All.ISMDustChem_SilicateNumberOfAtomsTable[key_elem] > sil_elem_abunds[k] / All.ISMDustChem_SilicateNumberOfAtomsTable[k]) {key_elem = k;}
                }
                key_mass = All.ISMDustChem_AtomicMassTable[All.ISMDustChem_SilicateMetallicityFieldIndexTable[key_elem]];
                key_num_atoms = All.ISMDustChem_SilicateNumberOfAtomsTable[key_elem];
                key_elem = All.ISMDustChem_SilicateMetallicityFieldIndexTable[key_elem];
                frac_of_max_sil = SphP[i].ISMDustChem_Dust_Species[0] / (P[i].Metallicity[key_elem] * All.ISMDustChem_EffectiveSilicateDustAtomicWeight/(key_num_atoms * key_mass));
                incl_frac = DMAX(DMIN(GALSF_ISMDUSTCHEM_VAR_IRON_INCL_FRAC*frac_of_max_sil,GALSF_ISMDUSTCHEM_VAR_IRON_INCL_FRAC),0.);
                SphP[i].ISMDustChem_Dust_Species[3] = (1.-incl_frac) * SphP[i].ISMDustChem_Dust_Metal[10];
                SphP[i].ISMDustChem_Dust_Species[NUM_ISMDUSTCHEM_SPECIES-1] = incl_frac * SphP[i].ISMDustChem_Dust_Metal[10];
            }
        }
    } // if (temp <= 300)
    
    /* Observed O depletions (Jenkins 2009) cannot be explained by silicate dust alone. So
     * throw extra oxygen into a reservoir to better match observations given O depletions vs
     * number density from Whittet (2010). We assume that this reservoir only holds as much O
     * as would be needed to match this trend scaled with the amount of silicate dust vs
     * the maximum allowable amount of silicate dust in the gas. So if the maximum amount
     * of silicate dust has formed than the O depletions should exactly match with
     * observations. This scaling allows for some variability in bursty environments.
     */
    if((GALSF_ISMDUSTCHEM_MODEL & 8) && (temp <= 300)) // We only add to the O reservoir when in an environment where dust can grow
    {
        double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS_CGS;    /* hydrogen number dens in cgs units */
        double D_O = 1. - 0.65441 / pow(nHcgs,0.103725);        /* expected fractional O depletion (upper limit) */
        double max_O_in_sil;                                    /* max O depletion due to silicates */
        double extra_O;                                         /* extra O that needs to be depleted to match observations */
        double frac_of_sil;                                     /* fraction of maximum amount of silicate present in gas */
        double O_in_CO;                                         /* mass fraction of O in CO, sets max for D_O */
        O_in_CO = SphP[i].ISMDustChem_C_in_CO * All.ISMDustChem_AtomicMassTable[4] / All.ISMDustChem_AtomicMassTable[2] / P[i].Metallicity[4];
        D_O = DMAX(0.,DMIN(D_O, 1.-O_in_CO)); // set depletion upper limit to O in CO
        
        // Now determine maximum possible silicate dust based on the least abundant element
        // This roughly scales with the fraction of the key element (usually Si) depleted into dust
        key_elem = 0;
        for(k=0;k<GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES;k++)
        {
            int index = All.ISMDustChem_SilicateMetallicityFieldIndexTable[k];
            num_dens[index] = rho * P[i].Metallicity[index] / (All.ISMDustChem_AtomicMassTable[index] * PROTONMASS_CGS);
            if (num_dens[All.ISMDustChem_SilicateMetallicityFieldIndexTable[key_elem]] / All.ISMDustChem_SilicateNumberOfAtomsTable[key_elem] > num_dens[index] / All.ISMDustChem_SilicateNumberOfAtomsTable[k]) key_elem = k;
        }
        key_mass = All.ISMDustChem_AtomicMassTable[All.ISMDustChem_SilicateMetallicityFieldIndexTable[key_elem]];
        key_num_atoms = All.ISMDustChem_SilicateNumberOfAtomsTable[key_elem];
        key_elem = All.ISMDustChem_SilicateMetallicityFieldIndexTable[key_elem];
        frac_of_sil = SphP[i].ISMDustChem_Dust_Species[0] / (P[i].Metallicity[key_elem] * All.ISMDustChem_EffectiveSilicateDustAtomicWeight/(key_num_atoms * key_mass));
        max_O_in_sil = P[i].Metallicity[key_elem] * ((All.ISMDustChem_SilicateNumberOfAtomsTable[0] * All.ISMDustChem_AtomicMassTable[4])/(key_num_atoms * key_mass));
        extra_O = frac_of_sil * D_O * P[i].Metallicity[4] - max_O_in_sil - SphP[i].ISMDustChem_Dust_Species[4];
        // If needed O depletion can't be attributed to silicate dust and what's already in the oxygen reservoir throw more oxygen into the reservoir
        if (extra_O > 0)
        {
            // Update creation source
            SphP[i].ISMDustChem_Dust_Source[source] += extra_O;
            // Now add the dust
            SphP[i].ISMDustChem_Dust_Metal[0] += extra_O;
            SphP[i].ISMDustChem_Dust_Metal[4] += extra_O;
            SphP[i].ISMDustChem_Dust_Species[4] += extra_O;
        }
    } // if using o reservoir and (temp <= 300)
    
    // Check if sputtering is delayed due to recent SNe
    if(SphP[i].ISMDustChem_DelayTimeSNeSputtering > 0) {SphP[i].ISMDustChem_DelayTimeSNeSputtering = DMAX(0,SphP[i].ISMDustChem_DelayTimeSNeSputtering-dtime_gyr);} // count off clock since last SNe
    else // Now determine amount of dust destroyed by thermal sputtering
    {
        T_ref = 2E6; avg_grain_radius = 0.032; /* um */ t_ref = 0.17; /* Gyr */
        sputter_timescale = t_ref * (avg_grain_radius / 0.1) / (rho*1E27) * (pow((T_ref/ temp), 2.5) + 1.);
        
        for (k=0;k<NUM_ISMDUSTCHEM_SPECIES;k++) species_yields[k] = 0.;
        for (k=0;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) dust_yields[k] = 0.;
        
        /******** SILICATE ********/
        dust_formula_mass = All.ISMDustChem_EffectiveSilicateDustAtomicWeight;
        dF = - dtime_gyr * (SphP[i].ISMDustChem_Dust_Species[0] / (sputter_timescale / 3.));
        dF = DMAX(-SphP[i].ISMDustChem_Dust_Species[0],DMIN(0,dF)); // can't destroy more dust then there is available
        species_yields[0] += dF;
        for (k=0;k<GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES;k++)
        {
            dust_yields[All.ISMDustChem_SilicateMetallicityFieldIndexTable[k]] += dF * (All.ISMDustChem_SilicateNumberOfAtomsTable[k] * All.ISMDustChem_AtomicMassTable[All.ISMDustChem_SilicateMetallicityFieldIndexTable[k]]) / dust_formula_mass;
        }
        /******** CARBONACEOUS ********/
        dF = - dtime_gyr * (SphP[i].ISMDustChem_Dust_Species[1] / (sputter_timescale / 3.));
        dF = DMAX(-SphP[i].ISMDustChem_Dust_Species[1],DMIN(0,dF)); // can't destroy more dust then there is available
        species_yields[1] += dF;
        dust_yields[2] += dF;
        /******** SILICONE CARBIDE ********/
        dust_formula_mass = All.ISMDustChem_AtomicMassTable[2] + All.ISMDustChem_AtomicMassTable[7];
        dF = - dtime_gyr * (SphP[i].ISMDustChem_Dust_Species[2] / (sputter_timescale / 3.));
        dF = DMAX(-SphP[i].ISMDustChem_Dust_Species[2],DMIN(0,dF)); // can't destroy more dust then there is available
        species_yields[2] += dF;
        dust_yields[2] += dF * All.ISMDustChem_AtomicMassTable[2] / dust_formula_mass;
        dust_yields[7] += dF * All.ISMDustChem_AtomicMassTable[7] / dust_formula_mass;
        /******** O RESERVOIR ********/
        if(GALSF_ISMDUSTCHEM_MODEL & 8) {
            dF = - dtime_gyr * (SphP[i].ISMDustChem_Dust_Species[4] / (sputter_timescale / 3.));
            dF = DMAX(-SphP[i].ISMDustChem_Dust_Species[4],DMIN(0,dF)); // can't destroy more dust then there is available
            species_yields[4] += dF;
            dust_yields[4] += dF;
        }
        /******** METALLIC IRON ********/
        if(GALSF_ISMDUSTCHEM_MODEL & 4) {avg_grain_radius = 0.0032; sputter_timescale = t_ref * (avg_grain_radius / 0.1) / (rho*1E27) * (pow((T_ref/ temp), 2.5) + 1.);} /* radius in um */
        dF = - dtime_gyr * (SphP[i].ISMDustChem_Dust_Species[3] / (sputter_timescale / 3.));
        dF = DMAX(-SphP[i].ISMDustChem_Dust_Species[3],DMIN(0,dF)); // can't destroy more dust then there is available
        species_yields[3] += dF;
        dust_yields[10] += dF;
        
        for(k=1;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) dust_yields[0] += dust_yields[k];
        
        // Update dust yields and creation source and deal with rounding errors when all dust is destroyed
        if (dust_yields[0] != 0.)
        {
            // Assume all dust sources are destroyed evenly
            for(k=0;k<NUM_ISMDUSTCHEM_SOURCES;k++) {SphP[i].ISMDustChem_Dust_Source[k] *= (1.+dust_yields[0]/SphP[i].ISMDustChem_Dust_Metal[0]);}
            // Update new dust mass
            for (k=0;k<NUM_ISMDUSTCHEM_SPECIES;k++) SphP[i].ISMDustChem_Dust_Species[k] += species_yields[k];
            // If all dust (silicates, carbonaceous, SiC, and free-flying iron) is destroyed zero everything to avoid rounding error
            int all_dest = 1;
            for (k=0;k<4;k++) {if(SphP[i].ISMDustChem_Dust_Species[k]>0) {all_dest = 0; break;}}
            if (all_dest)
            {
                for(k=0;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) dust_yields[k] = -SphP[i].ISMDustChem_Dust_Metal[k];
                for(k=0;k<NUM_ISMDUSTCHEM_SOURCES;k++) {SphP[i].ISMDustChem_Dust_Source[k] = 0.;}
                for(k=0;k<NUM_ISMDUSTCHEM_SPECIES;k++) {SphP[i].ISMDustChem_Dust_Species[k] = 0.;}
            }
            for (k=0;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {SphP[i].ISMDustChem_Dust_Metal[k] = DMAX(0,SphP[i].ISMDustChem_Dust_Metal[k]+dust_yields[k]);}
            if(GALSF_ISMDUSTCHEM_MODEL & 4) { // Update amount of free-flying iron and iron inclusions since some of the inclusions are released as silicates are sputtered. This scales with local amount of silicates. Note if all free-flying dust is destroyed then we assume all iron inclusions are also destroyed
                int key_elem = 0; double key_mass, key_num_atoms, frac_of_max_sil, incl_frac, sil_elem_abunds[GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES];
                for(k=0;k<GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES;k++)
                {
                    int index = All.ISMDustChem_SilicateMetallicityFieldIndexTable[k]; sil_elem_abunds[k] = P[i].Metallicity[index] / (All.ISMDustChem_AtomicMassTable[index] * PROTONMASS_CGS);
                    if (sil_elem_abunds[key_elem] / All.ISMDustChem_SilicateNumberOfAtomsTable[key_elem] > sil_elem_abunds[k] / All.ISMDustChem_SilicateNumberOfAtomsTable[k]) {key_elem = k;}
                }
                key_mass = All.ISMDustChem_AtomicMassTable[All.ISMDustChem_SilicateMetallicityFieldIndexTable[key_elem]];
                key_num_atoms = All.ISMDustChem_SilicateNumberOfAtomsTable[key_elem];
                key_elem = All.ISMDustChem_SilicateMetallicityFieldIndexTable[key_elem];
                frac_of_max_sil = SphP[i].ISMDustChem_Dust_Species[0] / (P[i].Metallicity[key_elem] * All.ISMDustChem_EffectiveSilicateDustAtomicWeight/(key_num_atoms * key_mass));
                incl_frac = DMAX(DMIN(GALSF_ISMDUSTCHEM_VAR_IRON_INCL_FRAC*frac_of_max_sil,GALSF_ISMDUSTCHEM_VAR_IRON_INCL_FRAC),0.);
                SphP[i].ISMDustChem_Dust_Species[3] = (1.-incl_frac) * SphP[i].ISMDustChem_Dust_Metal[10];
                SphP[i].ISMDustChem_Dust_Species[NUM_ISMDUSTCHEM_SPECIES-1] = incl_frac * SphP[i].ISMDustChem_Dust_Metal[10];
            }
        }
    }
#endif // dust species model options
}





#endif // GALSF_ISMDUSTCHEM_MODEL //
