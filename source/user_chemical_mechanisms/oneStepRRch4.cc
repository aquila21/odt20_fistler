/**
 * @file oneStepRRch4.cc
 */

using namespace std;
#include <vector>
#include <cmath>
#include <iostream>
#include<stdlib.h>


void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rrsp) {
/** inputs:
 * rho [=] kg/m3
 * temp [=] K
 * yi are species mass fractions
 * outputs:
 * rrsp [=] kmol/m3*s
 */

//    CH4 + 2 O2 => 2 H2O + CO2

    int iO2   = 0;
    int iCH4  = 1;
    int iCO2  = 2;
    int iH2O  = 3;
    int iN2   = 4;

    double rr_f = 1.3E8  * exp(-24358.0/temp);

    double cO2   = rho*yi[0]/31998.8;         // mol/cm3
    double cCH4  = rho*yi[1]/16042.6;         // mol/cm3

    double rr = rr_f *
                pow(fabs(cO2),   1.30) *
                pow(fabs(cCH4), -0.3+1.3*exp(-800*yi[iCH4]));

    rrsp[iO2]   = -2000.0*rr;                        // kmol/m3*s
    rrsp[iCH4]  = -1000.0*rr;
    rrsp[iH2O]  = 2000.0*rr;
    rrsp[iCO2]  = 1000.0*rr;
    rrsp[iN2]   = 0.0;
  }
