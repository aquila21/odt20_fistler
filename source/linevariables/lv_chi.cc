/**
 * @file lv_chi.cc
 * Header file for class lv_chi
 */


#include "lv_chi.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! lv_chi  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_chi::lv_chi(odtline  *line,
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
/** Set scalar dissipation rate (chi) \cond
 *  @param ipt \input optional point to compute at
 *  Chi = 2D(df/dx)^2
 *  Compute the derivative as follows:
 *        *       *                 *
 *       i-1      i                i+1
 *            a            b
 *  Linearly interpolate the derivatives at a and b to point i
 *  f' = ( d2*fa' + d1*fb' ) / (d1+d2) where d1 is dist between i and a
 *  and d2 is distance between b and i \endcond
 */

void lv_chi::setVar(const int ipt){

    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be -1" << endl;
        exit(0);
    }

    d.resize(odtl->ngrd);

    odtl->mixf->setVar();     // this may be redundant

    vector<double> gradZ(odtl->ngrd);
    vector<double> Dthm(odtl->ngrd);

    //-------------- Get thermal diffusivity

    for(int i=0; i<odtl->ngrd; i++) {
        odtl->odtc->setGasStateAtPt(i);
        double tcond = odtl->tran->thermalConductivity();
        double cp    = odtl->gas->cp_mass();
        Dthm.at(i) = tcond/odtl->rho->d.at(i)/cp;
    }

    //------------- Compute chi

    gradZ.at(0) = (odtl->mixf->d.at(1)-odtl->mixf->d.at(0))/(odtl->pos->d.at(1)-odtl->pos->d.at(0));
    double d1, d2;
    for(int i=1; i<odtl->ngrd-1; i++) {
        d1 = 0.5*(odtl->pos->d.at(i)  -odtl->pos->d.at(i-1));
        d2 = 0.5*(odtl->pos->d.at(i+1)-odtl->pos->d.at(i));
        gradZ.at(i) = (d2*d2*(odtl->mixf->d.at(i)  -odtl->mixf->d.at(i-1)) +
                    d1*d1*(odtl->mixf->d.at(i+1)-odtl->mixf->d.at(i)))/(2*d1*d2*(d1+d2));
    }
    gradZ.at(odtl->ngrd-1) = (odtl->mixf->d.at(odtl->ngrd-1)-odtl->mixf->d.at(odtl->ngrd-2)) /
                          (odtl->pos->d.at(odtl->ngrd-1)-odtl->pos->d.at(odtl->ngrd-2));

    for(int i=0; i<odtl->ngrd; i++)
        d.at(i) = 2.0*Dthm.at(i)*gradZ.at(i)*gradZ.at(i);

}

