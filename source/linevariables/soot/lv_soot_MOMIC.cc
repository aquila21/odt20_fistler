/**
 * @file lv_soot_MOMIC.cc
 * Header file for class lv_soot_MOMIC
 * @author Victoria B. Lansinger
 */

#include "lv_soot_MOMIC.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! lv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 *
 *    copied from lv_ygas.cc
 */

void lv_soot_MOMIC::getRhsSrc(const int ipt) {

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

    if(kMe==0) {                                                                // to save cost, compute needed terms for all lv_soot_MOMIC objects using this one.

        for(int k=0; k<nsvar; k++)                                              // for each mom/section,
            rrSvar.at(k).resize(odtl->ngrd);                                    // resize to match grid size

        for(int i=iS; i<=iE; i++) {                                             // loop over grid points

            // odtl->odtc->setGasStateAtPt(i);                                  // set gas state; uncomment to couple gas/soot
            odtl->odtc->enforceSootMom();                                       // make sure moments are positive or zero

            // get moments from odtline

            vector <double> M(nsvar,0.0);
            for (int k=0; k<nsvar; k++) {
                M[k] = odtl->svar[k]->d[i] * odtl->rho->d[i];
            }

            // determine how many moments to use

            int N = nsvar;                                                      // local number of moments
            downselectIfNeeded(i, M, N);                                        // downselect() will not change anything if all moment values >0

            // get chemical soot rates

            double Jnuc = chem->getNucleationRate(i);                           // #/m3*s
            double Kgrw = chem->getGrowthRate(i);                               // kg/m2*s
            double Koxi = chem->getOxidationRate(i);                            // kg/m2*s
            // NOTE: MOMIC has its own coagulation rate, implemented in getCoag

            // calculate MOMIC source terms; see Frenklach 2002 MOMIC paper

            double Acoef = M_PI*pow(abs(6.0/M_PI/chem->rho_soot),2.0/3.0);      // Acoef = kmol^2/3 / kg^2/3

            // nucleation terms
            vector<double> Mnuc(nsvar,0.0);                                     // Mnuc initialized with nsvar to be consistent with rrSvar sizing

            for (int k=0; k<N; k++) {                                           // values assigned with N (downselected number of moments)
                Mnuc[k] = pow(chem->Cmin*chem->MWc/chem->Na,k) * Jnuc;          // Nr = m_min^r * Jnuc
            }

            // growth terms
            vector<double> Mgrw(nsvar,0.0);

            for (int k=0; k<N; k++) {
                Mgrw[k] = k * Acoef * MOMIC(k-1.0/3.0,M) * Kgrw;                // Gk = r * M_(k-1/3)*Ap*ks
            }

            Mgrw[0] = 0.0;                                                      // by definition

            // oxidation terms
            vector<double> Moxi(nsvar,0.0);

            for (int k=0; k<N; k++) {
                Moxi[k] = k * Acoef * MOMIC(k-1.0/3.0,M) * Koxi;                // Xk = r * M_(k-1/3)*Ap*ks
            }

            Moxi[0] = 0.0;                                                      // by definition

            // coagulation terms
            vector<double> Mcoa(nsvar,0.0);

            if (odtl->odtp->coagulation_mech != "NONE") {
                for (int k=0; k<N; k++) {
                    Mcoa[k] = getCoag(odtl->temp->d[i],odtl->odtp->pres,odtl->dvisc->d[i],M,k);
                }
            }

            Mcoa[1] = 0.0;                                                      // by definition; redundant

            // combinine to make full source terms

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
/*! lagrangeInterp function
 *
 *      Calculates the Lagrange interpolated value from whole order moments.
 *
 *      @param x_i  \input      x value of desired interpolation
 *      @param x    \input      vector of x values to interpolate amongst
 *      @param y    \input      vector of y values to interpolate amongst
 *      @param y_i  \output     interpolated y value
 *
 */

double lv_soot_MOMIC::lagrangeInterp(double x_i, vector<double> x, vector<double> y) {

    double y_i = 0.0;

    for (int j = 0; j < x.size(); j++) {
        double L = 1.0;
        for (int m = 0; m < x.size(); m++) {
            if (m != j) {
                L = L * (x_i - x[m])/(x[j] - x[m]);
            }
        }
        y_i = y_i + y[j] * L;
    }

    return y_i;

}

////////////////////////////////////////////////////////////////////////////////
/*! MOMIC function
 *
 *      Calculates the desired fractional moment by lagrange interpolation
 *      between whole order moments. Because it uses log moments, it will crash
 *      if any moment is less than or equal to zero.
 *
 *      @param p     \input     desired interpolation value
 *      @param M     \input     vector of whole order moments
 *
 */

double lv_soot_MOMIC::MOMIC(double p, vector<double> M) {

    if (p == 0) {
        return M[0];
    }

    int size = M.size();
    if (p < 0 && M.size() != 2) {
        size = 3;
    }

    vector<double> log_mu(size,0.0);
    vector<double> x(size,0.0);

    for (int i = 0; i < size; i++) {
        log_mu[i] = log10(M[i] / M[0]);     // reduced moments
        x[i] = i;
    }

    double log_mu_p = lagrangeInterp(p,x,log_mu);

    return pow(10.0, log_mu_p) * M[0];

}

////////////////////////////////////////////////////////////////////////////////
/*! f_grid function
 *
 *      Calculates the grid function described in Frenklach 2002 MOMIC paper
 *      using lagrange interpolation between whole order moments
 *
 *      @param x     \input x grid point
 *      @param y     \input y grid point
 *      @param M     \input vector of whole order moments
 *
 */

double lv_soot_MOMIC::f_grid(int x, int y, vector<double> M) {

    double f1_0 = MOMIC(x-1.0/2.0,M)*MOMIC(y+1.0/6.0,M) + 2.0*MOMIC(x-1.0/6.0,M)*MOMIC(y-1.0/6.0,M) + MOMIC(x+1.0/6.0,M)*MOMIC(y-1.0/2.0,M);

    double f1_1 = MOMIC(x-1.0/2.0,M)*MOMIC(y+7.0/6.0,M) + 2.0*MOMIC(x-1.0/6.0,M)*MOMIC(y+5.0/6.0,M) + MOMIC(x+1.0/6.0,M)*MOMIC(y+1.0/2.0,M) +
                  MOMIC(x+1.0/2.0,M)*MOMIC(y+1.0/6.0,M) + 2.0*MOMIC(x+5.0/6.0,M)*MOMIC(y-1.0/6.0,M) + MOMIC(x+7.0/6.0,M)*MOMIC(y-1.0/2.0,M);

    if (y >= 4) {

        vector<double> temp_x(2,0.0);
        temp_x[0] = 0.0;
        temp_x[1] = 1.0;

        vector<double> temp_y(2,0.0);
        temp_y[0] = log10(f1_0);
        temp_y[1] = log10(f1_1);

        double value = lagrangeInterp(1.0/2.0, temp_x, temp_y);

        return pow(10.0, value);
    }

    double f1_2 =     MOMIC(x-1.0/2.0,M)*MOMIC(y+13.0/6.0,M) + 2.0*MOMIC(x-1.0 /6.0,M)*MOMIC(y+11.0/6.0,M) +     MOMIC(x+1.0 /6.0,M)*MOMIC(y+3.0/2.0,M) +
                  2.0*MOMIC(x+1.0/2.0,M)*MOMIC(y+7.0 /6.0,M) + 4.0*MOMIC(x+5.0 /6.0,M)*MOMIC(y+5.0 /6.0,M) + 2.0*MOMIC(x+7.0 /6.0,M)*MOMIC(y+1.0/2.0,M) +
                      MOMIC(x+3.0/2.0,M)*MOMIC(y+1.0 /6.0,M) + 2.0*MOMIC(x+11.0/6.0,M)*MOMIC(y-1.0 /6.0,M) +     MOMIC(x+13.0/6.0,M)*MOMIC(y-1.0/2.0,M);

    if (y >= 3) {

        vector<double> temp_x(3,0.0);
        temp_x[0] = 0.0;
        temp_x[1] = 1.0;
        temp_x[2] = 2.0;

        vector<double> temp_y(3,0.0);
        temp_y[0] = log10(f1_0);
        temp_y[1] = log10(f1_1);
        temp_y[2] = log10(f1_2);

        double value = lagrangeInterp(1.0/2.0, temp_x, temp_y);

        return pow(10.0, value);
    }

    double f1_3 =     MOMIC(x-1.0/2.0,M)*MOMIC(y+19.0/6.0,M) + 2.0*MOMIC(x-1.0 /6.0,M)*MOMIC(y+17.0/6.0,M) +     MOMIC(x+1.0 /6.0,M)*MOMIC(y+5.0/2.0,M) +
                  3.0*MOMIC(x+1.0/2.0,M)*MOMIC(y+13.0/6.0,M) + 6.0*MOMIC(x+5.0 /6.0,M)*MOMIC(y+11.0/6.0,M) + 3.0*MOMIC(x+7.0 /6.0,M)*MOMIC(y+3.0/2.0,M) +
                  3.0*MOMIC(x+3.0/2.0,M)*MOMIC(y+7.0 /6.0,M) + 6.0*MOMIC(x+11.0/6.0,M)*MOMIC(y+5.0 /6.0,M) + 3.0*MOMIC(x+13.0/6.0,M)*MOMIC(y+1.0/2.0,M) +
                      MOMIC(x+5.0/2.0,M)*MOMIC(y+1.0 /6.0,M) + 2.0*MOMIC(x+17.0/6.0,M)*MOMIC(y-1.0 /6.0,M) +     MOMIC(x+19.0/6.0,M)*MOMIC(y-1.0/2.0,M);

    vector<double> temp_x(4,0.0);
    temp_x[0] = 0.0;
    temp_x[1] = 1.0;
    temp_x[2] = 2.0;
    temp_x[3] = 3.0;

    vector<double> temp_y(4,0.0);
    temp_y[0] = log10(f1_0);
    temp_y[1] = log10(f1_1);
    temp_y[2] = log10(f1_2);
    temp_y[3] = log10(f1_3);

    double value = lagrangeInterp(1.0/2.0, temp_x, temp_y);

    return pow(10.0, value);

}

////////////////////////////////////////////////////////////////////////////////
/*! getCoag function
 *
 *      Calculates coagulation rate for MOMIC based on a weighted average of
 *      continuum and free-molecular values. See Frenklach's 2002 MOMIC paper.
 *      Adapted from python code by Alex Josephson.
 *
 *      @param T    \input  local gas temperature (K)
 *      @param P    \input  local gas pressure (Pa)
 *      @param mu   \input  local dynamic/absolute gas viscosity (kg/s/m)
 *      @param M    \input  vector of previous moment values
 *      @param r    \input  number of the moment to be calculated
 *
 */

double lv_soot_MOMIC::getCoag(double T, double P, double mu, vector<double> M, int r) {

    if (r == 1) {                                               // coagulation does not affect M1
        return 0;                                               // this is a shortcut to save a little computation
    }
    double kb    = 1.38064852E-23;
    double rho_soot = 1850.0;

    // Calculate Knudsen number to determine regime

    double mu_1     = M[1]/M[0];                                // average particle mass (kg)
    double d_g      = pow(6.0*kb*T/P/M_PI, 1.0/3.0);            // average gas molecular diameter (m)
    double d_p      = pow(6.0*mu_1/rho_soot/M_PI, 1.0/3.0);     // average particle diameter (m)
    double lambda_g = kb*T/(pow(2.0,0.5)*M_PI*pow(d_g,2.0)*P);  // gas mean free path (m)
    double Kn       = lambda_g/d_p;                             // Knudsen number

    // Continuum regime

    double Rate_C = 0.0;

    double K_C = 2.0*kb*T/(3.0*mu);
    double K_Cprime = 1.257*lambda_g*pow(M_PI*rho_soot/6.0,1.0/3.0);

    if (r == 0) {
        Rate_C = -K_C*(pow(M[0],2.0) + MOMIC(1.0/3.0,M)*MOMIC(-1.0/3.0,M) +
                  K_Cprime*(3.0*MOMIC(-1.0/3.0,M)*M[0] + MOMIC(2.0/3.0,M)*MOMIC(1.0/3.0,M)));
    }
    else {
        Rate_C = 0.0;
        for (int k=1; k<r; k++) {
            if (k > r-k) {
                Rate_C = Rate_C;
            }
            else {
                Rate_C = Rate_C + binomCoeff(r,k)*
                         (2.0*M[k]*M[r-k] + MOMIC(k+1.0/3.0,M)*MOMIC(r-k-1.0/3.0,M) + MOMIC(k-1.0/3.0,M)*MOMIC(r-k+1.0/3.0,M) +
                          2.0*K_Cprime* (2.0*MOMIC(k-1.0/3.0,M)*M[r-k] + M[k]*MOMIC(r-k-1.0/3.0,M) + MOMIC(k-2.0/3.0,M)*MOMIC(r-k+1.0/3.0,M)));
            }
        }
        Rate_C = 0.5*K_C*Rate_C;
    }

    // Free-molecular regime

    double Rate_F = 0.0;

    double K_f = 2.2*pow(3.0/(4.0*M_PI*rho_soot),2.0/3.0) * pow(8.0*M_PI*kb*T,1.0/2.0);

    if (r == 0) {
        Rate_F = -0.5*K_f*f_grid(0,0,M);
    }
    else {
        Rate_F = 0.0;
        for (int k=1; k<r; k++) {
            if (k > r-k) {
                Rate_F = Rate_F;
            }
            else {
                Rate_F = Rate_F + binomCoeff(r,k)*f_grid(k,r-k,M);
            }
        }
        Rate_F = 0.5*K_f*Rate_F;
    }

    return Rate_F/(1+1/Kn) + Rate_C/(1+Kn);

}

////////////////////////////////////////////////////////////////////////////////
/*! factorial function
 *
 *      Calculates the factorial of an integer n. Recursive.
 *
 *      @param n    \input      factorial input
 *      @param      \output     n!
 */

int lv_soot_MOMIC::factorial(int n) {

    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;

}

////////////////////////////////////////////////////////////////////////////////
/*! binomCoeff function
 *
 *      Calculates binomial coefficient
 *
 *      (r)      r!
 *      ( ) = ---------
 *      (k)   (r-k)!*k!
 */

double lv_soot_MOMIC::binomCoeff(int r, int k) {

    return 1.0 * factorial(r) / factorial(r-k) / factorial(k);      // multiply by 1.0 to force float division

}

////////////////////////////////////////////////////////////////////////////////
/*! downselectIfNeeded function
 *
 *      Checks the moment set can be used with MOMIC. Downselects if needed.
 *      Zero and negative moment values will trip the flag for downselection.
 *      The floor value is N = 2.
 *
 *      For the case where M0 <= 0, we set N = 0, which means that the source
 *      terms will be initalized but not calculated such that the rates are all
 *      zero. If M1 <= 0, we assign it a value based on M0 and a lognormal
 *      distribution that matches the initialized profile in the odtcase. This
 *      is a workaround since the lagrange interpolation can't handle M1 = 0.
 *
 *      @param ipt  \input  grid point to evaluate at
 *      @param &M   \input  moment set
 *      @param &N   \input  number of downselected moments
 *
 */

void lv_soot_MOMIC::downselectIfNeeded(int ipt, vector<double> &M, int &N) {

    // CHECK: M0 <= 0.0

    if (M[0] <= 0.0) {
        N = 0;
        return;
    }

    // CHECK: M1 <= 0.0

    if (M[1] <= 0.0) {
        double M0 = 1.0;
        double sigL = 3.0;                              // sigL and mavg should be same as in odtcase
        double mavg = 1.0E-21;
        M[1] = M0 * mavg * exp(0.5 * pow(sigL,2.0));    // give M1 a value based on M0 and lognormal dist.
        odtl->svar[1]->d[ipt] = M[1];
    }

    // CHECK: all remaining moments

    bool zeros = false;                                 // flag for breaking do/while loop

    do {

        if (N <= 2) {                                   // will downselect to 2 moments, but not further
            break;
        }

        zeros = false;                                  // reset flag

        for (int i=0; i<N; i++) {
           if (M[i] <= 0.0) {
               zeros = true;
           }
        }

        if (zeros == true) {                            // if flagged, downselect by two moments
            N = N - 1;
        }

    } while (zeros != 0);                               // keep going until we get no flag

    M.resize(N);                                        // resize M based on downselected N

}
