/**
 * @file lv_hr.cc
 * Header file for class lv_hr
 */


#include "lv_hr.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
// Declare the global function prototype so it can be used in this source file

void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rr);


////////////////////////////////////////////////////////////////////////////////
/*! lv_hr  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_hr::lv_hr(odtline  *line,
             const      string s,
             const bool Lt,
             const bool Lo) {

    odtl          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(odtl->ngrd, 0.0);

}

////////////////////////////////////////////////////////////////////////////////
/*! Set heat release rate (J/m3*s) from the gas state
 *  @param ipt \input optional point to compute at
 */

void lv_hr::setVar(const int ipt){

    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be -1" << endl;
        exit(0);
    }

    d.resize(odtl->ngrd);

    int nsp = odtl->gas->nSpecies();
    vector<double> rr(nsp);
    vector<double> yi(nsp);
    vector<double> hsp(nsp);
    double GasConstant = 8314.47215;             // J/kmol*K

    for(int i=0; i<odtl->ngrd; i++){
        odtl->odtc->setGasStateAtPt(i);
#ifdef PROBLEMSPECIFICRR
        odtl->gas->getMassFractions( &yi[0] );
        getProblemSpecificRR(odtl->gas->density(), odtl->gas->temperature(), odtl->odtp->pres, &yi.at(0), &rr.at(0));
#else
        odtl->gas->getNetProductionRates(&rr.at(0));
#endif
        odtl->gas->getEnthalpy_RT(&hsp.at(0));               // non-dimensional enthalpy
        d.at(i) = 0.0;
        for(int k=0; k<nsp; k++)
            d.at(i) -= rr.at(k)/odtl->gas->molecularWeight(k)*hsp.at(k)*odtl->gas->temperature()*GasConstant;
    }
}

