/**
 * @file lv_soot_LOGN.cc
 * Header file for class lv_soot_LOGN
 * @author Victoria B. Lansinger
 */

#include "lv_soot_LOGN.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! lv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 *
 *    copied from lv_ygas.cc
 */

void lv_soot_LOGN::getRhsSrc(const int ipt) {

    if(!L_transported)
        return;

    rhsSrc.resize(odtl->ngrd, 0.0);                    // resize rhsSrc

    static vector<vector<double> > rrSvar(nsvar);      // temp storage for moment rates [nsvar][ngrd]

    int iS, iE;
    if(ipt==-1) {
        iS = 0;
        iE = (odtl->ngrd)-1;
    }
    else {
        iS = ipt;
        iE = ipt;
    }

    for(int k=0; k<nsvar; k++)                         // for each mom/section,
        rrSvar.at(k).resize(odtl->ngrd);               // resize to match grid size

    if(kMe==0) {                                       // to save cost, compute needed terms for all lv_soot_LOGN objects using this one.

        for(int k=0; k<nsvar; k++)                     // for each mom/section,
            rrSvar.at(k).resize(odtl->ngrd);           // resize to match grid size

        for(int i=iS; i<=iE; i++) {                    // loop over grid points

            // odtl->odtc->setGasStateAtPt(i);         // keep commented for decoupled soot
            odtl->odtc->enforceSootMom();

            //--------- chemical soot rates

            double Jnuc  = chem->getNucleationRate(i);                                              // #/m3*s
            double Koxi  = chem->getOxidationRate(i);                                               // kg/m2*s
            double Kgrw  = chem->getGrowthRate(i);                                                  // kg/m2*s
            double Kfc   = chem->getCoagulationRate(i);                                             // #/m3*s

            //--------- nucleation terms

            double N0 = Jnuc;                                                                       // #/m3*s
            double N1 = Jnuc*chem->Cmin*chem->MWc/chem->Na;                                         // kg/m3*s
            double N2 = Jnuc*pow(chem->Cmin*chem->MWc/chem->Na,2);                                  // kg2/m3*s

            //--------- growth terms

            double G0 = 0.0;                                                                        // zero by definition, #/m3*s
            double G1 = Kgrw * M_PI*pow(6.0/chem->rho_soot/M_PI,2.0/3.0) * Mk(2.0/3.0,i);                   // kg/m3*s
            double G2 = Kgrw * M_PI*pow(6.0/chem->rho_soot/M_PI,2.0/3.0) * Mk(5.0/3.0,i) * 2;               // kg2/m3*s

            //--------- oxidation terms

            double X0 = 0.0;                                                                        // zero by definition, #/m3*s
            double X1 = Koxi * M_PI*pow(6.0/chem->rho_soot/M_PI,2.0/3.0) * Mk(2.0/3.0,i);                   // kg/m3*s
            double X2 = Koxi * M_PI*pow(6.0/chem->rho_soot/M_PI,2.0/3.0) * Mk(5.0/3.0,i) * 2;               // kg2/m3*s

            //--------- coagulation terms

            double b    = odtl->odtp->b_coag;                                                       // this is exact for monodisperse distribution

            double C0    = -Kfc * b * (Mk(0.0,i)*Mk(1.0/6.0,i) + 2.0*Mk(1.0/3.0,i)*Mk(-1.0/6.0,i) +
                                       Mk(2.0/3.0,i)*Mk(-1.0/2.0,i));                                   // #/m3*s
            double C1    = 0.0;                                                                     // zero by definition, kg/m3*s
            double C2    = 2 * Kfc * b * (Mk(1.0,i)*Mk(7.0/6.0,i) + 2*Mk(4.0/3.0,i)*Mk(5.0/6.0,i) +
                                       Mk(5.0/3.0,i)*Mk(1.0/2.0,i));                                   // kg2/m3*s

            //--------- combinine to make source terms

            rrSvar[0][i] = (N0 + G0 - X0 + C0) / odtl->rho->d[i];                                   // (#/m3*s)/rho = #/kg*s
            rrSvar[1][i] = (N1 + G1 - X1 + C1) / odtl->rho->d[i];                                   // (kg-soot/m3*s)/rho = kg-soot/kg*s
            rrSvar[2][i] = (N2 + G2 - X2 + C2) / odtl->rho->d[i];                                   // (kg-soot/m3*s)ˆ2/rho = kg-sootˆ2/kg*s
        }
    }

    for(int i=iS; i<=iE; i++)
        rhsSrc.at(i) = rrSvar.at(kMe).at(i);

    if(odtl->odtp->Lspatial)
        for(int i=iS; i<=iE; i++)
            rhsSrc.at(i) /= odtl->uvel->d.at(i);

}

////////////////////////////////////////////////////////////////////////////////
/*! Mk function
 *    Calculates fractional moments
 *
 *    @param exp     \input     fractional moment to compute, corresponds to exponent
 *    @param wts     \input     weights
 *    @param absc    \input     abscissas
 *    @param Mk    \output fractional moment value
 *
 */

double lv_soot_LOGN::Mk(double k, const int i) {

    double M0 = odtl->svar[0]->d[i] * odtl->rho->d[i];                               // M0 = #/m3
    double M1 = odtl->svar[1]->d[i] * odtl->rho->d[i];                               // M1 = rhoYs = kg/m3
    double M2 = odtl->svar[2]->d[i] * odtl->rho->d[i];                               // M2 = kg2/m3

    double Mk = pow(M0, 1 - 1.5*k + 0.5*pow(k,2)) *
                pow(M1, 2*k - pow(k,2)) *
                pow(M2,0.5*pow(k,2) - 0.5*k);

    return Mk;

}
