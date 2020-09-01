/**
 * @file lv_soot_QMOM.cc
 * Header file for class lv_soot_QMOM
 * @author Victoria B. Lansinger
 */

#include "lv_soot_QMOM.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

void pdAlg(int nm, int np, vector<double> &mu, vector<double> &wts, vector<double> &absc );
void wheeler(const vector<double> &m, int N, vector<double> &w, vector<double> &x );
void adaptive_wheeler(const vector<double> &m, int N, const vector<double> &rmin, const double &eabs, int &Nout, vector<double> &w, vector<double> &x );
void adaptiveWheelerAlgorithm(const std::vector<double>& moments, std::vector<double>& w, std::vector<double>& x,const double& rMin, const double& eAbs);

////////////////////////////////////////////////////////////////////////////////
/*! getRhsSrc function
 *
 *      lv source term part of the rhs function. Uses QMOM. Modeled off of
 *      lv_ygas class.
 *
 *      @param ipt \input optional point to compute source at.
 *
 */

void lv_soot_QMOM::getRhsSrc(const int ipt) {

    if(!L_transported)
        return;

    rhsSrc.resize(odtl->ngrd, 0.0);                                             // resize rhsSrc

    static vector<vector<double> > rrSvar(nsvar);                               // temp storage for moment rates [nsvar][ngrd]

    int iS, iE;
    if(ipt==-1) {
        iS = 0;
        iE = odtl->ngrd-1;
    }
    else {
        iS = ipt;
        iE = ipt;
    }

    for(int k=0; k<nsvar; k++)                                                  // for each mom/section,
        rrSvar.at(k).resize(odtl->ngrd);                                        // resize to match grid size

    if(kMe==0) {                                                                // to save cost, compute needed terms for all lv_soot_QMOM objects using this one.

        for(int k=0; k<nsvar; k++)                                              // for each mom/section,
            rrSvar.at(k).resize(odtl->ngrd);                                    // resize to match grid size

        for(int i=iS; i<=iE; i++) {                                             // loop over grid points

            // odtl->odtc->setGasStateAtPt(i);                                  // set gas state; uncomment to couple gas/soot
            odtl->odtc->enforceSootMom();

            // get moments from odtline
            vector <double> M(nsvar,0.0);
            for (int k=0; k<nsvar; k++) {
                M[k] = odtl->svar[k]->d[i] * odtl->rho->d[i];
            }

            // initialize weights and abscissas
            wts.resize(nsvar/2);
            absc.resize(nsvar/2);

            for (int k=0; k<nsvar/2; k++) {
                wts[k] = 0.0;
                absc[k] = 0.0;
            }

            // get weights and abscissas
            getWtsAbs(M, wts, absc, i);                                         // PD and wheeler algorithms called in here

            // get surface area coefficient
            double Acoef = M_PI*pow(abs(6.0/M_PI/chem->rho_soot),2.0/3.0);      // Acoef = kmol^2/3 / kg^2/3

            // get chemical soot rates
            double Jnuc = chem->getNucleationRate(i);                           // #/m3*s
            double Kgrw = chem->getGrowthRate(i);                               // kg/m2*s
            double Koxi = chem->getOxidationRate(i);                            // kg/m2*s

            // calculate QMOM source terms; see DOL thesis pg. 67

            // nucleation terms
            vector<double> Mnuc(nsvar,0.0);

            for (int k=0; k<nsvar; k++) {
                Mnuc[k] = pow(chem->Cmin*chem->MWc/chem->Na,k) * Jnuc;          // Nr = m_min^r * Jnuc
            }

            // growth terms
            vector<double> Mgrw(nsvar,0.0);

            for (int k=0; k<nsvar; k++) {
                Mgrw[k]    = k * Acoef * Mk(k-1.0/3.0) * Kgrw;                 // Gk = r * M_(k-1/3)*Ap*ks
            }

            Mgrw[0] = 0.0;                                                     // by definition

            // oxidation terms
            vector<double> Moxi(nsvar,0.0);

            for (int k=0; k<nsvar; k++) {
                Moxi[k]    = k * Acoef * Mk(k-1.0/3.0) * Koxi;           // Xk = r * M_(k-1/3)*Ap*ks
            }

            Moxi[0] = 0.0;                                                      // by definition

            // coagulation terms
            vector<double> Mcoa(nsvar,0.0);

            for (int r=0; r<nsvar; r++) {                                       // See DOL thesis pg. 60, equation 2.96
                for (int p=0; p<nsvar/2; p++) {
                    for (int q=0; q<nsvar/2; q++) {
                        for (int k=0; k<r+1; k++) {
                            Mcoa[r] += wts[p] * wts[q] * chem->getCoagulationRate(i,p,q,wts,absc) *
                                      (-pow(absc[p],r) + 0.5*chem->binomCoeff(r,k) * pow(absc[p],k) * pow(absc[q],r-k));
                        }
                    }
                }
            }

            Mcoa[1] = 0.0;                                                      // by definition

            // combinine to make source terms
            for (int k=0; k<nsvar; k++) {
                rrSvar[k][i] = (Mnuc[k] + Mgrw[k] - Moxi[k] + Mcoa[k]) / odtl->rho->d[i];
            }

        }
    }

     for(int i=iS; i<=iE; i++)
         rhsSrc.at(i) = rrSvar.at(kMe).at(i);                                   // assign to odtl's rhsSrc variable

     if(odtl->odtp->Lspatial)
         for(int i=iS; i<=iE; i++)
             rhsSrc.at(i) /= odtl->uvel->d.at(i);

}

////////////////////////////////////////////////////////////////////////////////
/*! Mk function
 *
 *      Calculates fractional moments from weights and abscissas.
 *
 *      @param exp  \input  fractional moment to compute, corresponds to exponent
 *      @param wts  \input  weights
 *      @param absc \input  abscissas
 *      @param Mk   \output fractional moment value
 *
 */

double lv_soot_QMOM::Mk(double exp) {

    double Mk = 0;

    for(int k=0; k<nsvar/2; k++) {
        if (wts[k] == 0 || absc[k] == 0)
            return 0;
        else
            Mk += wts[k] * pow(absc[k],exp);
    }

    return Mk;

}

////////////////////////////////////////////////////////////////////////////////
/*! getWtsAbs function
 *
 *      Calculates weights and abscissas from moments using PD algorithm or
 *      wheeler algorithm.
 *
 *      @param M        \input  vector of moments
 *      @param wts      \input  weights
 *      @param absc     \input  abscissas
 *      @param ipt      \input  grid point to evaluate at
 *
 *      Notes:
 *      - Use wheeler over pdalg whenever possible.
 *      - wts and abs DO NOT change size; if we downselect to a smaller number
 *      of moments, the extra values are set at and stay zero
 *      - using w_temp and a_temp means we don't have to resize wts and absc,
 *      which is more convenient when wts and absc are used to reconstitute
 *      moment source terms.
 */

void lv_soot_QMOM::getWtsAbs(vector<double> M, vector<double> &wts, vector<double> &absc, const int ipt) {

    for (int k=0; k<nsvar; k++) {                                   // if any moments are zero, return with zero wts and absc
        if (M[k] <= 0.0)
            return;
    }

    int N = nsvar;                                                  // local nsvar
    bool negs = false;                                              // flag for downselecting if there are negative wts/abs

    vector<double> w_temp(N/2,0.0);
    vector<double> a_temp(N/2,0.0);

    do {                                                            // downselection loop

        negs = false;                                               // reset flag

        for (int k=0; k<N/2; k++) {                                 // reinitialize wts and abs with zeros
            w_temp[k] = 0.0;
            a_temp[k] = 0.0;
        }

        if (N == 2) {                                               // in 2 moment case, return monodisperse output
           wts[0]  = M[0];
           absc[0] = M[1]/M[0];
           return;
        }

        //pdAlg(N,N/2,ymom,wts,absc);                               // PD algorithm
        wheeler(M, N/2, w_temp, a_temp);                           // wheeler algorithm

        for (int k=0; k<N/2; k++) {
            if (w_temp[k] < 0.0 || a_temp[k] < 0.0 || a_temp[k] > 1.0)
                negs = true;
        }

        if (negs == true) {                                         // if we found negative values
            N = N - 2;                                              // downselect to two fewer moments and try again
            w_temp.resize(N/2);
            a_temp.resize(N/2);
        }

    } while (negs == true);                                         // end of downselection loop

    for (int k = 0; k < w_temp.size(); k++) {                       // assign temporary variables to output
        wts[k]  = w_temp[k];
        absc[k] = a_temp[k];
    }

}
