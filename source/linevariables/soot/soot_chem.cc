/** Implementation file for class soot_chem
 *
 *      @file soot_chem.cc
 *      @author Victoria B. Lansinger
 */

#include "soot_chem.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
/** Initialization function
 *
 *      @param      p_odtl  \input pointer to odtline object
 */

void soot_chem::init(odtline *p_odtl) {

    //-------------- define constants and values

    odtl             = p_odtl;                                                  ///< pointer to odtline object

    rho_soot         = odtl->odtp->rho_soot;                                    ///< solid soot density, kg/m3
    Cmin             = odtl->odtp->Cmin;                                        ///< minimum number of carbon atoms in a soot particle
    nsvar            = odtl->odtp->nsvar;                                       ///< number of soot variables (i.e. moments, bins)
    b                = odtl->odtp->b_coag;                                      ///< coagulation constant; 1/sqrt(2) < b < 1
    PSD_method       = odtl->odtp->PSD_method;                                  ///< flag for PSD solver method (MONO, QMOM, MOMIC, etc.)
    nucleation_mech  = odtl->odtp->nucleation_mech;                             ///< flag for nucleation chemistry
    growth_mech      = odtl->odtp->growth_mech;                                 ///< flag for surface growth chemistry
    oxidation_mech   = odtl->odtp->oxidation_mech;                              ///< flag for oxidation chemistry
    coagulation_mech = odtl->odtp->coagulation_mech;                            ///< flag for coagulation rate method
    MWc              = odtl->gas->atomicWeight(odtl->gas->elementIndex("C"));   ///< molar weight carbon/soot = 12.01 kg/kmol

    eps_c            = 9.0;                                                      ///< coagulation constant
    Na               = 6.02214086E26;                                            ///< Avogadro's constant, particles/kmol
    kb               = 1.38064852E-23;                                           ///< Boltzmann constant, J/K
    Rg               = 8314.46;                                                  ///< Universal gas constant, J/kmol*K

    //-------------- populate list of gas species indices

    int isp;

    isp = odtl->gas->speciesIndex("C2H2");                                      // Copy/paste to add species to this list
    isp = (isp >= 0) ? isp : odtl->gas->speciesIndex("c2h2");
    i_c2h2 = isp;

    isp = odtl->gas->speciesIndex("O2");
    isp = (isp >= 0) ? isp : odtl->gas->speciesIndex("o2");
    i_o2 = isp;

    isp = odtl->gas->speciesIndex("H");
    isp = (isp >= 0) ? isp : odtl->gas->speciesIndex("h");
    i_h = isp;

    isp = odtl->gas->speciesIndex("H2");
    isp = (isp >= 0) ? isp : odtl->gas->speciesIndex("h2");
    i_h2 = isp;

    isp = odtl->gas->speciesIndex("OH");
    isp = (isp >= 0) ? isp : odtl->gas->speciesIndex("oh");
    i_oh = isp;

    isp = odtl->gas->speciesIndex("H2O");
    isp = (isp >= 0) ? isp : odtl->gas->speciesIndex("h2o");
    i_h2o = isp;

    isp = odtl->gas->speciesIndex("CO");
    isp = (isp >= 0) ? isp : odtl->gas->speciesIndex("co");
    i_co = isp;

    //-------------- TO DO: test that the species present are sufficient for the desired soot mechanism

}

////////////////////////////////////////////////////////////////////////////////
/*! getNucleationRate function
 *
 *      Calls appropriate function for nucleation chemistry
 *      based on the nucleation_mech flag.
 *
 *      @param ipt      /input grid point to evaluate at
 */

double soot_chem::getNucleationRate(int ipt) {

    if (nucleation_mech == "NONE") {
        return 0;
    }
    else if (nucleation_mech == "LL") {
        return nucleation_LL(odtl->rho->d[ipt], odtl->temp->d[ipt],
                             odtl->ysp[i_c2h2]->d[ipt], odtl->gas->molecularWeight(i_c2h2));
    }
    else {
        cout << endl << "ERROR: Invalid soot nucleation mechanism." << endl;
        exit(0);
    }

    return 0;

}

////////////////////////////////////////////////////////////////////////////////
/*! getGrowthRate function
 *
 *      Calls appropriate function for surface growth chemistry
 *      based on the growth_mech flag.
 *
 *      @param ipt      /input grid point to evaluate at
 */

double soot_chem::getGrowthRate(int ipt) {

    if (growth_mech == "NONE") {
        return 0;
    }
    else if (growth_mech == "LL") {
        return growth_LL(odtl->rho->d[ipt], odtl->temp->d[ipt],
                         odtl->ysp[i_c2h2]->d[ipt], odtl->gas->molecularWeight(i_c2h2),
                         odtl->svar[0]->d[ipt] * odtl->rho->d[ipt], odtl->svar[1]->d[ipt] * odtl->rho->d[ipt]);
    }
    else if (growth_mech == "LIN") {
        return growth_Lindstedt(odtl->rho->d[ipt], odtl->temp->d[ipt],
                                odtl->ysp[i_c2h2]->d[ipt], odtl->gas->molecularWeight(i_c2h2));
    }
    else if (growth_mech == "HACA") {
        return growth_HACA(odtl->rho->d[ipt], odtl->temp->d[ipt], odtl->odtp->pres,
                           odtl->ysp[i_c2h2]->d[ipt], odtl->gas->molecularWeight(i_c2h2), odtl->ysp[i_o2 ]->d[ipt], odtl->gas->molecularWeight(i_o2 ),
                           odtl->ysp[i_h   ]->d[ipt], odtl->gas->molecularWeight(i_h   ), odtl->ysp[i_h2 ]->d[ipt], odtl->gas->molecularWeight(i_h2 ),
                           odtl->ysp[i_oh  ]->d[ipt], odtl->gas->molecularWeight(i_oh  ), odtl->ysp[i_h2o]->d[ipt], odtl->gas->molecularWeight(i_h2o),
                           odtl->svar[0]->d[ipt] * odtl->rho->d[ipt], odtl->svar[1]->d[ipt] * odtl->rho->d[ipt]);
    }
    else {
        cout << endl << "ERROR: Invalid soot growth mechanism." << endl;
        exit(0);
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
/*! getOxidationRate function
 *
 *      Calls appropriate function for oxidation chemistry
 *      based on the oxidation_mech flag.
 *
 *      @param ipt      /input grid point to evaluate at
 */

double soot_chem::getOxidationRate(int ipt) {

    if (oxidation_mech == "NONE") {
        return 0;
    }
    else if (oxidation_mech == "LL") {
        return oxidation_LL(odtl->rho->d[ipt], odtl->temp->d[ipt],
                            odtl->ysp[i_o2]->d[ipt], odtl->gas->molecularWeight(i_o2));
    }
    else if (oxidation_mech == "HACA") {
        return oxidation_HACA(odtl->rho->d[ipt], odtl->temp->d[ipt], odtl->odtp->pres,
                              odtl->ysp[i_c2h2]->d[ipt], odtl->gas->molecularWeight(i_c2h2), odtl->ysp[i_o2 ]->d[ipt], odtl->gas->molecularWeight(i_o2 ),
                              odtl->ysp[i_h   ]->d[ipt], odtl->gas->molecularWeight(i_h   ), odtl->ysp[i_h2 ]->d[ipt], odtl->gas->molecularWeight(i_h2 ),
                              odtl->ysp[i_oh  ]->d[ipt], odtl->gas->molecularWeight(i_oh  ), odtl->ysp[i_h2o]->d[ipt], odtl->gas->molecularWeight(i_h2o),
                              odtl->svar[0]->d[ipt] * odtl->rho->d[ipt], odtl->svar[1]->d[ipt] * odtl->rho->d[ipt]);
    }
    else if (oxidation_mech == "LEE_NEOH") {
        return oxidation_Lee_Neoh(odtl->rho->d[ipt], odtl->temp->d[ipt], odtl->odtp->pres,
                                  odtl->ysp[i_o2]->d[ipt], odtl->gas->molecularWeight(i_o2), odtl->ysp[i_oh]->d[ipt], odtl->gas->molecularWeight(i_oh), odtl->gas->meanMolecularWeight(),
                                  odtl->svar[0]->d[ipt] * odtl->rho->d[ipt], odtl->svar[1]->d[ipt] * odtl->rho->d[ipt]);
    }
    else {
        cout << endl << "ERROR: Invalid soot oxidation mechanism." << endl;
        exit(0);
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
/*! getCoagulationRate function
 *
 *      Calls appropriate function for coagulation chemistry
 *      based on the coagulation_mech flag. Returns the value of
 *      the collision rate function beta in #/m3*s.
 *
 *      @param ipt      /input  grid point to evaluate at
 *      @param p        /input  abscissa index of m1
 *      @param q        /input  abscissa index of m2
 *      @param wts      /input  vector of weights of the soot PSD
 *      @param absc     /input  vector of abscissas of the soot PSD
 */

double soot_chem::getCoagulationRate(int ipt, int p, int q, vector<double> &wts, vector<double> &absc) {

    if (coagulation_mech == "NONE") {
        return 0;
    }
    else if (coagulation_mech == "LL") {
        return coagulation_LL(absc[p], absc[q], odtl->temp->d[ipt]);
    }
    else if (coagulation_mech == "FUCHS") {
        return coagulation_Fuchs(absc[p], absc[q], odtl->temp->d[ipt], odtl->rho->d[ipt], odtl->dvisc->d[ipt], odtl->gas->meanMolecularWeight());
    }
    else {
        cout << endl << "ERROR: Invalid soot coagulation mechanism." << endl;
        exit(0);
    }

    return 0;

}

////////////////////////////////////////////////////////////////////////////////
/*! Nucleation by Leung_Lindstedt (1991)
 *
 *      Rate from Leung & Lindstedt (1991), Comb. & Flame 87:289-305.
 *      Returns chemical nucleation rate in #/m3*s.
 *
 *      @param rho      /input  local gas density (kg/m3)
 *      @param T        /input  local gas temperature (K)
 *      @param y_c2h2   /input  local gas mass fraction of C2H2
 *      @param MW_c2h2  /input  molecular weight of C2H2 (kg/kmol)
 */

double soot_chem::nucleation_LL(double rho, double T, double y_c2h2, double MW_c2h2) {

    double cC2H2 = rho * y_c2h2 / MW_c2h2;              // kmol/m3
    double Rnuc =  0.1E5 * exp(-21100/T) * cC2H2;       // kmol/m^3*s

    return Rnuc * 2 * Na / Cmin;                        // #/m3*s

}

////////////////////////////////////////////////////////////////////////////////
/*! Growth by Lindstedt (1994)
 *
 *      Rate from Bockhorn (1994) pg. 417, "Simplified Soot Nucleation and Surface Growth Steps..."
 *      Returns chemical surface growth rate in kg/m2*s.
 *
 *      @param rho      /input  local gas density (kg/m3)
 *      @param T        /input  local gas temperature (K)
 *      @param y_c2h2   /input  local gas mass fraction of C2H2
 *      @param MW_c2h2  /input  molecular weight of C2H2 (kg/kmol)
 */

double soot_chem::growth_Lindstedt(double rho, double T, double y_c2h2, double MW_c2h2) {

    double cC2H2 = rho * y_c2h2 / MW_c2h2;              // kmol/m3

    return 750 * exp(-12100/T) * cC2H2 * 2*MWc;         // kg/m^2*s

}

////////////////////////////////////////////////////////////////////////////////
/*! Growth by Leung_Lindstedt (1991)
 *
 *      Rate from Leung & Lindstedt (1991), Comb. & Flame 87:289-305.
 *      Returns chemical surface growth rate in kg/m2*s.
 *
 *      @param rho      /input  local gas density (kg/m3)
 *      @param T        /input  local gas temperature (K)
 *      @param y_c2h2   /input  local gas mass fraction of C2H2
 *      @param MW_c2h2  /input  molecular weight of C2H2 (kg/kmol)
 *      @param M0       /input  local soot number density (#/m3)
 *      @param M1       /input  local soot mass density (kg/m3)
 */

double soot_chem::growth_LL(double rho, double T, double y_c2h2, double MW_c2h2, double M0, double M1) {

    double Am2m3 = M_PI * pow(6*M1/M_PI/rho_soot,2.0/3.0) * pow(M0,1.0/3.0);      // m^2_soot/m^3_total
    double cC2H2 = rho * y_c2h2 / MW_c2h2;                                                  // kmol/m3

    if (Am2m3 <= 0.0)
        return 0;

    return 0.6E4 * exp(-12100/T) * cC2H2/sqrt(Am2m3) * 2*MWc;                               // kg/m^2*s

}

////////////////////////////////////////////////////////////////////////////////
/*! Oxidation by Leung_Lindstedt (1991)
 *
 *      Rate from Leung & Lindstedt (1991), Comb. & Flame 87:289-305.
 *      Returns chemical soot oxidation rate in kg/m2*s.
 *
 *      @param rho      /input  local gas density (kg/m3)
 *      @param T        /input  local gas temperature (K)
 *      @param y_o2     /input  local gas mass fraction of O2
 *      @param MW_02    /input  molecular weight of O2 (kg/kmol)
 */

double soot_chem::oxidation_LL(double rho, double T, double y_o2, double MW_o2) {

    double cO2 = rho * y_o2 / MW_o2;                            // kmol/m3

    return 0.1E5 * sqrt(T) * exp(-19680/T) * cO2 * MWc;         // kg/m^2*s

}

////////////////////////////////////////////////////////////////////////////////
/*! Oxidation by Lee et al. + Neoh
 *
 *      Rates from Lee et al. (1962) Comb. & Flame 6:137-145 and Neoh (1981)
 *      "Soot oxidation in flames" in Particulate Carbon Formation During
 *      Combustion book
 *
 *      Returns chemical soot oxidation rate in kg/m2*s.
 *
 *      @param rho      /input  local gas density (kg/m3)
 *      @param T        /input  local gas temperature (K)
 *      @param P        /input  gas pressure (atm)
 *      @param y_o2     /input  local gas mass fraction of O2
 *      @param MW_02    /input  molecular weight of O2 (kg/kmol)
 *      @param y_oh     /input  local gas mass fraction of OH
 *      @param MW_oh    /input  molecular weight of OH (kg/kmol)
 *      @param MW_mean  /input  mean gas molecular weight (kg/kmol)
 *      @param M0       /input  local soot number density (#/m3)
 *      @param M1       /input  local soot mass density (#/m3)
 */

double soot_chem::oxidation_Lee_Neoh(double rho, double T, double P, double y_o2, double MW_o2, double y_oh, double MW_oh,
                                     double MW_mean, double M0, double M1) {

    double pO2 = y_o2 * MW_mean / MW_o2 * P / 101325.0;                                 // partial pressure of O2
    double pOH = y_oh * MW_mean / MW_oh * P / 101325.0;                                 // partial pressure of OH

    double cO2 = rho * y_o2 / MW_o2;                                                    // kmol/m3

    double Am2m3 = M_PI * pow(abs(6*M1/M_PI/rho_soot),2.0/3.0) * pow(abs(M0),1.0/3.0);  // m^2_soot/m^3_total

    return 1.085E4*pO2/sqrt(T)*exp(-1.977824E4/T)/1000.0 + 1290.0*0.13*pOH/sqrt(T);     // kg/m^2*s

}

////////////////////////////////////////////////////////////////////////////////
/*! Coagulation by Leung_Lindstedt
 *
 *      Kfc rate from Leung & Lindstedt (1991), Comb. & Flame 87:289-305
 *      Expression for beta from DOL thesis pg. 60
 *
 *      Notes:
 *      - eps_c is the same as LL Ca.
 *      - Kfc (#/m3*s) = -4*Kf*b * M0**11/6 * M1**1/6, where Kf = eps_c*(6*kb*T/rhos)**1/2 * (3/4/pi/rhos)**1/6
 *      - For monodisperse case: b = 1/sqrt(2) is exact. For all others, 1/sqrt(2) < b < 1.
 *      - This function returns collision rate beta, NOT Kfc
 *      - pow() can't handle zero to a negative power, so we check and return 0 if necessary.
 *
 *      @param m1   \input  first particle size (kg)
 *      @param m2   \input  second particle size (kg)
 *      @param T    \input  local temperature (K)
 */

double soot_chem::coagulation_LL (double m1, double m2, double T) {

    if (m1 == 0.0 || m2 == 0.0) {     // check for zeros so as not to crash the pow function
        return 0;
    }

    double Kf = eps_c * pow(6*kb*T/rho_soot,0.5) * pow(3.0/(4.0*M_PI*rho_soot),1.0/6.0);

    double paren = pow(m1,1.0/6.0) + 2.0*pow(m1,-1.0/6.0)*pow(m2,1.0/3.0) + pow(m1,-1.0/2.0)*pow(m2,2.0/3.0) +
                   pow(m1,2.0/3.0)*pow(m2,-1.0/2.0) + 2.0*pow(m1,1.0/3.0)*pow(m2,-1.0/6.0) + pow(m2,1.0/6.0);


    return Kf * b * paren;

}

////////////////////////////////////////////////////////////////////////////////
/*! Coagulation by Fuchs
 *
 *      Rate comes from Seinfeld and Pandis Atmospheric Chemistry book (2016), pg. 548.
 *      Details and clarification in Fuchs' Mechanics of Aerosols book (1964)
 *
 *      This function returns the collision rate beta in #/m3*s, NOT K values.
 *
 *      @param m1       \input  first particle size (kg)
 *      @param m2       \input  second particle size (kg)
 *      @param T        \input  local temperature (K)
 *      @param rho      \input  local gas density (kg/m3)
 *      @param mu       \input  local gas dynamic viscosity (kg/m*s)
 *      @param MW       \input  average molecular weight of the gas (kg/kmol)
 */

double soot_chem::coagulation_Fuchs(double m1, double m2, double T, double rho, double mu, double MW) {

    double Dp1 = pow(6.0*m1/M_PI/rho_soot, 1.0/3.0);
    double Dp2 = pow(6.0*m2/M_PI/rho_soot, 1.0/3.0);

    double c1 = sqrt(8.0*kb*T/M_PI/m1);
    double c2 = sqrt(8.0*kb*T/M_PI/m2);

    double mfp_g = mu/rho*sqrt(M_PI*MW/2.0/Rg/T);

    double Kn1 = 2.0*mfp_g/Dp1;
    double Kn2 = 2.0*mfp_g/Dp2;

    double Cc1 = 1 + Kn1*(1.257 + 0.4*exp(-1.1/Kn1));
    double Cc2 = 1 + Kn2*(1.257 + 0.4*exp(-1.1/Kn2));

    double D1 = kb*T*Cc1/(3.0*M_PI*mu*Dp1);
    double D2 = kb*T*Cc2/(3.0*M_PI*mu*Dp2);

    double l1 = 8.0*D1/M_PI/c1;
    double l2 = 8.0*D2/M_PI/c2;

    double g1 = sqrt(2.0)/3.0/Dp1/l1*( pow(Dp1+l1,3.0) - pow(Dp1*Dp1 + l1*l1, 3.0/2.0) ) - sqrt(2.0)*Dp1;
    double g2 = sqrt(2.0)/3.0/Dp2/l2*( pow(Dp2+l2,3.0) - pow(Dp2*Dp2 + l2*l2, 3.0/2.0) ) - sqrt(2.0)*Dp2;

    return 2.0*M_PI*(D1+D2)*(Dp1+Dp2) / ((Dp1+Dp2)/(Dp1+Dp2+2.0*sqrt(g1*g1+g2*g2)) + 8.0*(D1+D2)/sqrt(c1*c1+c2*c2)/(Dp1+Dp2));

}

//////////////////////////////////////////////////////////////////////////////////
/*! Growth by HACA
 *
 *      See Appel, Bockhorn, & Frenklach (2000), Comb. & Flame 121:122-136.
 *      For details, see Franklach and Wang (1990), 23rd Symposium, pp. 1559-1566.
 *
 *      Parameters for steric factor alpha updated to those given in Balthasar
 *      and Franklach (2005) Comb. & Flame 140:130-145.
 *
 *      Returns the chemical soot growth rate in kg/m2*s.
 *
 *      @param rho      /input  local gas density (kg/m3)
 *      @param T        /input  local gas temperature (K)
 *      @param P        /input  gas pressure (atm)
 *      @param y_c2h2   /input  local gas mass fraction of C2H2
 *      @param MW_c2h2  /input  molecular weight of C2H2 (kg/kmol)
 *      @param y_o2     /input  local gas mass fraction of O2
 *      @param MW_o2    /input  molecular weight of O2 (kg/kmol)
 *      @param y_h      /input  local gas mass fraction of H
 *      @param MW_h     /input  molecular weight of H (kg/kmol)
 *      @param y_oh     /input  local gas mass fraction of OH
 *      @param MW_oh    /input  molecular weight of OH (kg/kmol)
 *      @param y_h2o    /input  local gas mass fraction of H2O
 *      @param MW_h2o   /input  molecular weight of H2O (kg/kmol)
 *      @param M0       /input  local soot number density (#/m3)
 *      @param M1       /input  local soot mass density (kg/m3)
 */

double soot_chem::growth_HACA(double rho,    double T,    double P,    double y_c2h2, double MW_c2h2, double y_o2,  double MW_o2,
                              double y_h,    double MW_h, double y_h2, double MW_h2,  double y_oh,    double MW_oh, double y_h2o,
                              double MW_h2o, double M0,   double M1) {

    double cC2H2 = rho * y_c2h2 / MW_c2h2;      // kmol/m3
    double cO2   = rho * y_o2   / MW_o2;        // kmol/m3
    double cH    = rho * y_h    / MW_h;         // kmol/m3
    double cH2   = rho * y_h2   / MW_h2;        // kmol/m3
    double cOH   = rho * y_oh   / MW_oh;        // kmol/m3
    double cH2O  = rho * y_h2o  / MW_h2o;       // kmol/m3

    // calculate alpha, other constants
    double RT       = 1.9872036E-3 * T;         // R (=) kcal/mol
    double chi_soot = 2.3E15;                   // (=) sites/cm^2
    double a_param  = 33.167 - 0.0154 * T;      // a parameter for calculating alpha
    double b_param  = -2.5786 + 0.00112 * T;    // b parameter for calculating alpha

    // calculate raw HACA reaction rates
    double fR1 = 4.2E13 * exp(-13.0 / RT) * cH / 1000;
    double rR1 = 3.9E12 * exp(-11.0 / RT) * cH2 / 1000;
    double fR2 = 1.0E10 * pow(T,0.734) * exp(-1.43 / RT) * cOH / 1000;
    double rR2 = 3.68E8 * pow(T,1.139) * exp(-17.1 / RT) * cH2O /1000;
    double fR3 = 2.0E13 * cH / 1000;
    double fR4 = 8.00E7 * pow(T,1.56) * exp(-3.8 / RT) * cC2H2 / 1000;
    double fR5 = 2.2E12 * exp(-7.5 / RT) * cO2 / 1000;
    double fR6 = 1290.0 * 0.13 * P * (cOH/rho*MW_oh) / sqrt(T);                 // gamma = 0.13 from Neoh et al.

    // Steady state calculation of chi for soot radical; see Frenklach 1990 pg. 1561
    double denom = rR1 + rR2 + fR3 + fR4 + fR5;
    double chi_rad = 0.0;
    if(denom != 0.0)
        chi_rad = 2 * chi_soot * (fR1 + fR2 + fR6) / denom;                     // sites/cm^2

    double alpha = 1.0;                                                         // alpha = fraction of available surface sites
    if (M0 != 0.0)
        alpha = tanh(a_param/log10(M1/M0)+b_param);
    if (alpha < 0.0)
        alpha = 1.0;

    double c_soot_H   = alpha * chi_soot * 1E4;                                 // sites/m2-mixture
    double c_soot_rad = alpha * chi_rad  * 1E4;                                 // sites/m2-mixture

    return (fR5*c_soot_rad + fR6*c_soot_H) / Na * 2 * MWc;                      // kg/m2*s
}

//////////////////////////////////////////////////////////////////////////////////
/*! Oxidation by HACA
 *
 *      See Appel, Bockhorn, & Frenklach (2000), Comb. & Flame 121:122-136.
 *      For details, see Franklach and Wang (1990), 23rd Symposium, pp. 1559-1566.
 *
 *      Parameters for steric factor alpha updated to those given in Balthasar
 *      and Franklach (2005) Comb. & Flame 140:130-145.
 *
 *      Returns the chemical soot oxidation rate in kg/m2*s.
 *
 *      @param rho      /input  local gas density (kg/m3)
 *      @param T        /input  local gas temperature (K)
 *      @param P        /input  gas pressure (atm)
 *      @param y_c2h2   /input  local gas mass fraction of C2H2
 *      @param MW_c2h2  /input  molecular weight of C2H2 (kg/kmol)
 *      @param y_o2     /input  local gas mass fraction of O2
 *      @param MW_o2    /input  molecular weight of O2 (kg/kmol)
 *      @param y_h      /input  local gas mass fraction of H
 *      @param MW_h     /input  molecular weight of H (kg/kmol)
 *      @param y_oh     /input  local gas mass fraction of OH
 *      @param MW_oh    /input  molecular weight of OH (kg/kmol)
 *      @param y_h2o    /input  local gas mass fraction of H2O
 *      @param MW_h2o   /input  molecular weight of H2O (kg/kmol)
 *      @param M0       /input  local soot number density (#/m3)
 *      @param M1       /input  local soot mass density (kg/m3)
 */

double soot_chem::oxidation_HACA(double rho,    double T,    double P,    double y_c2h2, double MW_c2h2, double y_o2,  double MW_o2,
                                 double y_h,    double MW_h, double y_h2, double MW_h2,  double y_oh,    double MW_oh, double y_h2o,
                                 double MW_h2o, double M0,   double M1) {

    double cC2H2 = rho * y_c2h2 / MW_c2h2;      // kmol/m3
    double cO2   = rho * y_o2   / MW_o2;        // kmol/m3
    double cH    = rho * y_h    / MW_h;         // kmol/m3
    double cH2   = rho * y_h2   / MW_h2;        // kmol/m3
    double cOH   = rho * y_oh   / MW_oh;        // kmol/m3
    double cH2O  = rho * y_h2o  / MW_h2o;       // kmol/m3

    // calculate alpha, other constants
    double RT       = 1.9872036E-3 * T;         // R (=) kcal/mol
    double chi_soot = 2.3E15;                   // (=) sites/cm^2
    double a_param  = 33.167 - 0.0154 * T;      // a parameter for calculating alpha
    double b_param  = -2.5786 + 0.00112 * T;    // b parameter for calculating alpha

    // calculate raw HACA reaction rates
    double fR1 = 4.2E13 * exp(-13.0 / RT) * cH / 1000;
    double rR1 = 3.9E12 * exp(-11.0 / RT) * cH2 / 1000;
    double fR2 = 1.0E10 * pow(T,0.734) * exp(-1.43 / RT) * cOH / 1000;
    double rR2 = 3.68E8 * pow(T,1.139) * exp(-17.1 / RT) * cH2O /1000;
    double fR3 = 2.0E13 * cH / 1000;
    double fR4 = 8.00E7 * pow(T,1.56) * exp(-3.8 / RT) * cC2H2 / 1000;
    double fR5 = 2.2E12 * exp(-7.5 / RT) * cO2 / 1000;
    double fR6 = 1290.0 * 0.13 * P * (cOH/rho*MW_oh) / sqrt(T);                         // gamma = 0.13 from Neoh et al.

    // Steady state calculation of chi for soot radical; see Frenklach 1990 pg. 1561
    double denom = rR1 + rR2 + fR3 + fR4 + fR5;
    double chi_rad = 0.0;
    if(denom != 0.0)
        chi_rad = 2 * chi_soot * (fR1 + fR2 + fR6) / denom;                             // sites/cm^2

    double alpha = 1.0;                                                                 // alpha = fraction of available surface sites
    if (M0 != 0.0)
        alpha = tanh(a_param/log10(M1/M0)+b_param);
    if (alpha < 0.0)
        alpha = 1.0;

    double c_soot_H   = alpha * chi_soot * 1E4;                                         // sites/m2-mixture
    double c_soot_rad = alpha * chi_rad  * 1E4;                                         // sites/m2-mixture

    double Roxi = -fR1*c_soot_H + rR1*c_soot_rad - fR2*c_soot_H + rR2*c_soot_rad +
                   fR3*c_soot_rad + fR4*c_soot_rad - fR6*c_soot_H;                      // #-available-sites/m2-mix*s
    return Roxi / Na * MWc;                                                             // kg/m2*s

}

////////////////////////////////////////////////////////////////////////////////
/*! factorial function
 *
 *      Calculates the factorial of an integer n.
 *
 *      @param n \input
 */

int soot_chem::factorial(int n) {

    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;

}

////////////////////////////////////////////////////////////////////////////////
/*! binomCoeff function
 *
 *      Calculates binomial coefficient.
 *
 *      (r)      r!
 *      ( ) = ---------
 *      (k)   (r-k)!*k!
 */

double soot_chem::binomCoeff(int r, int k) {

    return 1.0 * factorial(r) / factorial(r-k) / factorial(k);      // multiply by 1.0 to force float division

}
