/*
Special routines for computing thermodynamic properties of the hydrogen
molecule, assuming a standard 3:1 ortho:para mixture and accounting for
rotational, vibration, and translational degrees of freedom.

Equations follow Boley 2007, ApJ, 656, L89
*/

#include "../allvars.h"
#include <math.h>

#ifdef EOS_SUBSTELLAR_ISM

void hydrogen_molecule_zrot_mixture(double temp, double result[3]) {
    /*
    Rotational partition function of hydrogen molecule and derived quantities,
    considering a 3:1 mixture of ortho- and parahydrogen that cannot efficiently
    come into equilibrium.

    Parameters
    ----------
    temp: double
        Temperature in K
    ortho_frac: double
        Fraction of ortho-H2 (default is 3:1 ortho:para mixture)
    result: double[3]
        Stores the partition function value, the average rotational energy per
    molecule, and the heat capacity per molecule at constant volume.
    */

    const double EPSILON = 2.220446049250313e-16;
    const double THETA_ROT = 85.4;  // in K
    const double ortho_frac = 0.75; // 3:1 mixture
    const double para_frac = 1 - ortho_frac;
    const double x = THETA_ROT / temp;
    const double expmx = exp(-x);
    const double expmx4 = pow(expmx, 4);

    double error = 1e100;
    double z[2] = {0}; // index 0 for para, 1 for ortho
    double dz_dtemp[2] = {0};
    double d2z_dtemp2[2] = {0};
    double zterm[2] = {0};

    z[0] = zterm[0] = 1.0;
    z[1] = zterm[1] = 9.0;
    double expterm = expmx4 * expmx * expmx;

    // Summing over rotational levels
    int j = 2;
    double dzterm, d2zterm;
    while (error > EPSILON) {
        int s = j % 2;
        zterm[s] *= (2 * j + 1) * expterm / (2 * j - 3);
        int jjplusone = j * (j + 1);
        if (s == 1) { // ortho
            dzterm = (jjplusone - 2) * x * zterm[1];
            d2zterm = ((jjplusone - 2) * x - 2) * dzterm;
        } else { // para
            dzterm = jjplusone * x * zterm[0];
            d2zterm = (jjplusone * x - 2) * dzterm;
        }
        z[s] += zterm[s];
        dz_dtemp[s] += dzterm;
        d2z_dtemp2[s] += d2zterm;
        double err0 = zterm[0] / z[0];
        double err1 = zterm[1] / z[1];
        if (err1 > err0) {
            error = err1;
        } else {
            error = err0;
        }
        expterm *= expmx4;
        j++;
    }

    result[0] = exp(para_frac * log(z[0]) + ortho_frac * log(z[1]));                                       // partition function
    result[1] = BOLTZMANN_CGS * temp * (para_frac * dz_dtemp[0] / z[0] + ortho_frac * dz_dtemp[1] / z[1]); // mean energy per molecule
    result[2] = BOLTZMANN_CGS * (ortho_frac * (2 * dz_dtemp[1] + d2z_dtemp2[1] - dz_dtemp[1] * dz_dtemp[1] / z[1]) / z[1] +
                                 para_frac * (2 * dz_dtemp[0] + d2z_dtemp2[0] - dz_dtemp[0] * dz_dtemp[0] / z[0]) / z[0]); // heat capacity
}

void hydrogen_molecule_zvib(double temp, double result[3]) {
    /*
    Vibrational partition function of hydrogen molecule and derived quantities.

    Parameters
    ----------
    temp: double
        Temperature in K
    result: double[3]
        Stores the partition function value, the average rotational energy per
    molecule, and the heat capacity per molecule at constant volume.
    */
    const double THETA_VIB = 6140;
    const double x = THETA_VIB / temp;
    result[0] = -1.0 / expm1(-x);
    result[1] = BOLTZMANN_CGS * THETA_VIB / expm1(x);
    result[2] = THETA_VIB * result[0] * result[1] / (temp * temp);
}

void hydrogen_molecule_partitionfunc(double temp, double result[3]) {
    /*
    Thermodynamic quantities derived from the partition function of the
    hydrogen molecule.

    Parameters
    ----------
    temp: double
        Temperature in K
    result: double[3]
        Stores the the average rotational energy per molecule in erg,
        the heat capacity per molecule at constan volume in erg/K,
        and the adiabatic index
    */

    double zrot[4], zvib[4];
    hydrogen_molecule_zrot_mixture(temp, zrot);
    hydrogen_molecule_zvib(temp, zvib);
    double etot = 1.5 * BOLTZMANN_CGS * temp; // translation
    double cv = 1.5 * BOLTZMANN_CGS;
    etot += zrot[1]; // rotation
    cv += zrot[2];
    etot += zvib[1]; // vibration
    cv += zvib[2];
    double gamma = (cv / BOLTZMANN_CGS + 1) / (cv / BOLTZMANN_CGS);
    result[0] = etot;
    result[1] = cv;
    result[2] = gamma;
}

double hydrogen_molecule_energy(double temp) {
    /*
    Average energy of a H2 molecule in thermodynamic equilibrium

    Parameters
    ----------
    temp: double
        Temperature in K

    Returns
    -------
    etot: double
        Average energy of a H2 molecule of temperture T in erg
    */

    if (temp < 12.5) {
        return 1.5 * BOLTZMANN_CGS * temp; // only translation
    } else if (temp > 1e5) {
        return 3.5 * BOLTZMANN_CGS * temp; // all DOF excited
    }

    double zrot[3], zvib[3];
    hydrogen_molecule_zrot_mixture(temp, zrot);
    hydrogen_molecule_zvib(temp, zvib);
    double etot = 1.5 * BOLTZMANN_CGS * temp; // translation
    etot += zrot[1];                          // rotation
    etot += zvib[1];                          // vibration
    return etot;
}

double hydrogen_molecule_gamma(double temp) {
    /*
    First adiabatic index of hydrogen molecule assuming
    a 3:1 ortho:para mixture

    Parameters
    ----------
    temp: double
        Temperature in K

    Returns
    -------
    gamma: double
        Adiabatic index
    */

    if (temp < 12.5) {
        return 5. / 3; // only translation
    } else if (temp > 1e5) {
        return 9. / 7; // all DOF excited
    }

    double zrot[3], zvib[3];
    hydrogen_molecule_zrot_mixture(temp, zrot);
    hydrogen_molecule_zvib(temp, zvib);
    double cv = 1.5;               // translation
    cv += zrot[2] / BOLTZMANN_CGS; // rotation
    cv += zvib[2] / BOLTZMANN_CGS; // vibration
    return (cv + 1) / cv;
}

#endif // EOS_SUBSTELLAR_ISM
