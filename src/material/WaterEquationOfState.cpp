#include "WaterEquationOfState.h"

double WaterEquationOfState::rho(double P, double T) const noexcept{
    // METROLOGIE
    // DOSSIER N° 18 : Calcul de la masse volumique de l'eau
    // Projet DensiCal
    T -= 273.15;
    if (T < -12) return rho(P, -12);
    if (T > 500) return rho(P, 500);
    double a1 = -3.983035;
    double a2 = 301.797;
    double a3 = 522528.9;
    double a4 = 69.34881;
    double a5 = 999.974950;
    double rho = a5 * (1 - (T+a1)*(T+a1)*(T+a2)/(a3*(T+a4)));

    P -= 101325;
    double c1 = 5.074e-10;
    double c2 = -3.26e-12;
    double c3 = 4.16e-15;
    rho *= (1.0 + (c1+c2*T+c3*T*T)*P);
    return rho;
}

double WaterEquationOfState::drho_dP(double P, double T) const noexcept{
    T -= 273.15;
    if (T < -12 || T > 500) return 0;
    double a1 = -3.983035;
    double a2 = 301.797;
    double a3 = 522528.9;
    double a4 = 69.34881;
    double a5 = 999.974950;
    double drho = a5 * (1 - (T+a1)*(T+a1)*(T+a2)/(a3*(T+a4)));

    double c1 = 5.074e-10;
    double c2 = -3.26e-12;
    double c3 = 4.16e-15;
    drho *= (c1+c2*T+c3*T*T);
    return drho;
}

double WaterEquationOfState::drho_dT(double P, double T) const noexcept{
    T -= 273.15;
    if (T < -12 || T > 500) return 0;
    double a1 = -3.983035;
    double a2 = 301.797;
    double a3 = 522528.9;
    double a4 = 69.34881;
    double a5 = 999.974950;
    double Tpart = a5 * (1 - (T+a1)*(T+a1)*(T+a2)/(a3*(T+a4)));
    double dTpart = -a5/a3 * ((T+a1)*(3*T+a1+2*a2)*(T+a4) - (T+a1)*(T+a1)*(T+a2))/((T+a4)*(T+a4));

    P -= 101325;
    double c1 = 5.074e-10;
    double c2 = -3.26e-12;
    double c3 = 4.16e-15;
    double Ppart = (1.0 + (c1+c2*T+c3*T*T)*P);
    double dPpart = (c2+2*c3*T)*P;
    return Tpart*dPpart + Ppart*dTpart;
}