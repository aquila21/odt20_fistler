/**
 *  @file soot_chem.h
 *    Header file for class soot_chem
 */

#ifndef SOOT_CHEM_H
#define SOOT_CHEM_H

#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////
/** Class implementing soot_chem, no parent object
 *
 *      @author David O. Lignell
 *      @author Victoria B. Lansinger
 */

class soot_chem {

    //////////////////// DATA MEMBERS //////////////////////

    public:

        odtline      *odtl;                ///< pointer to line object

        double        rho_soot;            ///< solid soot density, kg/m3
        int           Cmin;                ///< minimum number of carbon atoms in a soot particle
        double        eps_c;               ///< coagulation constant
        int           nsvar;               ///< number of soot variables, i.e. moments
        double        b;                   ///< Coagulation constant, 1/sqrt(2) < b < 1
        string        PSD_method;          ///< MONO, QMOM, etc.
        string        nucleation_mech;     ///< soot nucleation chemistry flag
        string        growth_mech;         ///< soot growth chemistry flag
        string        oxidation_mech;      ///< soot oxidation chemistry flag
        string        coagulation_mech;    ///< soot coagulation chemistry flag

        double        Na;                  ///< Avogadro's constant
        double        MWc;                 ///< molar weight carbon/soot = 12.01 kg/kmol
        double        kb;                  ///< Boltzmann constant
        double        Rg;                  ///< Universal gas constant

        int           i_c2h2;              // soot gas species indices
        int           i_o2;
        int           i_h;
        int           i_h2;
        int           i_oh;
        int           i_h2o;
        int           i_co;

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        double getNucleationRate  (int ipt);
        double getGrowthRate      (int ipt);
        double getOxidationRate   (int ipt);
        double getCoagulationRate (int ipt, int p, int q, vector<double> &wts, vector<double> &absc);

        int    factorial  (int n);
        double binomCoeff (int r, int k);

    private:

        double nucleation_LL      (double rho, double T, double y_c2h2, double MW_c2h2);

        double growth_Lindstedt   (double rho, double T, double y_c2h2, double MW_c2h2);
        double growth_LL          (double rho, double T, double y_c2h2, double MW_c2h2, double M0, double M1);
        double growth_HACA        (double rho,    double T,    double P,    double y_c2h2, double MW_c2h2, double y_o2,  double MW_o2,
                                   double y_h,    double MW_h, double y_h2, double MW_h2,  double y_oh,    double MW_oh, double y_h2o,
                                   double MW_h2o, double M0,   double M1);

        double oxidation_LL       (double rho, double T, double y_o2, double MW_o2);
        double oxidation_HACA     (double rho,    double T,    double P,    double y_c2h2, double MW_c2h2, double y_o2,  double MW_o2,
                                   double y_h,    double MW_h, double y_h2, double MW_h2,  double y_oh,    double MW_oh, double y_h2o,
                                   double MW_h2o, double M0,   double M1);
        double oxidation_Lee_Neoh (double rho, double T, double P, double y_o2, double MW_o2, double y_oh, double MW_oh,
                                   double MW_mean, double M0, double M1);

        double coagulation_LL     (double m1, double m2, double T);
        double coagulation_Fuchs  (double m1, double m2, double T, double rho, double mu, double MW);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        soot_chem(){};    // constructor
        void init(odtline *p_odtl);
        ~soot_chem(){}
};

#endif
