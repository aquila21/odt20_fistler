/**
 * @file lv_soot_MONO.cc
 * Header file for class lv_soot_MONO
 * @author Victoria B. Lansinger
 */

#include "lv_soot_MONO.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! lv source term part of the rhs function
 *
 *  @param ipt \input optional point to compute source at.
 *
 *  adapted from lv_ygas.cc
 */

void lv_soot_MONO::getRhsSrc(const int ipt) {

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

    if(kMe==0) {                                       // to save cost, compute needed terms for all lv_soot_MONO objects using this one.

        for(int k=0; k<nsvar; k++)                     // for each mom/section,
            rrSvar.at(k).resize(odtl->ngrd);           // resize to match grid size

        for(int i=iS; i<=iE; i++) {                    // loop over grid points

            // odtl->odtc->setGasStateAtPt(i);         // keep commented for decoupled soot
            odtl->odtc->enforceSootMom();

            double M0    = odtl->svar[0]->d[ipt] * odtl->rho->d[ipt];                               // M0 = #/m3
            double M1    = odtl->svar[1]->d[ipt] * odtl->rho->d[ipt];                               // M1 = rhoYs = kg/m3

            double Am2m3 = M_PI*pow(abs(6*M1/M_PI/chem->rho_soot),2.0/3.0)*pow(abs(M0),1.0/3.0);    // m^2_soot / m^3_total

            //--------- weights and abcissas for coagulation term

            wts.resize(1);                                                                          // simplification of a 2-moment QMOM case
            absc.resize(1);                                                                         // so we can use one consistent coagulation interface

            if (M0 <= 0.0) {
                wts[0] = 0.0;
                absc[0] = 0.0;
            }
            else {
                wts[0] = M0;                                                                            // defined weights and abscissas for the monodisperse case
                absc[0] = M1/M0;
            }

            //--------- chemical soot rates

            double Jnuc  = chem->getNucleationRate(i);                                              // #/m3*s
            double Koxi  = chem->getOxidationRate(i);                                               // kg/m2*s
            double Kgrw  = chem->getGrowthRate(i);                                                  // kg/m2*s

            //--------- nucleation terms

            double N0 = Jnuc;                                                                       // #/m3*s
            double N1 = Jnuc*chem->Cmin*chem->MWc/chem->Na;                                         // kg/m3*s

            //--------- growth terms

            double G0 = 0.0;                                                                        // zero by definition, #/m3*s
            double G1 = Kgrw*Am2m3;                                                                 // kg/m3*s

            //--------- oxidation terms

            double X0 = 0.0;                                                                        // zero by definition, #/m3*s
            double X1 = Koxi*Am2m3;                                                                 // kg/m3*s

            //--------- coagulation terms

            // double C0    = -4.0*Kfc*b*pow(abs(M0),11.0/6.0)*pow(abs(M1),1.0/6.0);                // #/m3*s; expression from DOL thesis pg. 63
            double C0 = -0.5 * wts[0] * wts[0] * chem->getCoagulationRate(i,0,0,wts,absc);          // #/m3*s; simplification of generalized expression
            double C1 = 0.0;                                                                        // zero by definition, kg/m3*s

            //--------- combinine to make source terms

            rrSvar[0][i] = (N0 + G0 - X0 + C0) / odtl->rho->d[i];                                   // (#/m3*s)/rho = #/kg*s
            rrSvar[1][i] = (N1 + G1 - X1 + C1) / odtl->rho->d[i];                                   // (kg-soot/m3*s)/rho = kg-soot/kg*s
        }
    }

    for(int i=iS; i<=iE; i++)
        rhsSrc.at(i) = rrSvar.at(kMe).at(i);

    if(odtl->odtp->Lspatial)
        for(int i=iS; i<=iE; i++)
            rhsSrc.at(i) /= odtl->uvel->d.at(i);

}
