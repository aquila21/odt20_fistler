/**
 * @file lv_rho.cc
 * Header file for class lv_rho
 */


#include "lv_rho.h"
#include "odtline.h"

////////////////////////////////////////////////////////////////////////////////
/*! lv_rho  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_rho::lv_rho(odtline    *line,
               const      string s,
               const bool Lt,
               const bool Lo) {

    odtl          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(odtl->ngrd, odtl->odtp->rho0);

}

////////////////////////////////////////////////////////////////////////////////
/*! lv_rho merger2cells function
 *
 * Function presumes that the variable being merged is a quantity per unit mass.
 * Merging conservatively: (rho*phi*dx)_merged = (rho*phi*dx)_1 + (rho*phi*dx)_2
 * Then solve for phi_merged.
 *
 * @param imrg    \input merge cells imrg and imrg+1
 * @param imrg \input merge cells imrg and imrg+1
 * @param m1   \input mass in cell imrg
 * @param m2   \input mass in cell imrg
 * @param LconstVolume \input (for posf, default is false)
 */

void lv_rho::merge2cells(const int    imrg,
                         const double m1,
                         const double m2,
                         const bool   LconstVolume) {

    odtl->odtc->setGasStateAtPt(imrg);
    d.at(imrg) = odtl->gas->density();
    d.erase(d.begin() + imrg+1);

}

////////////////////////////////////////////////////////////////////////////////
/*! lv_rho setVar function
 *  @param ipt \input optional point to compute at
 */

void lv_rho::setVar(const int ipt){

    d.resize(odtl->ngrd, odtl->odtp->rho0);
    if(ipt == -1)
        for(int i=0; i<odtl->ngrd; i++) {
            odtl->odtc->setGasStateAtPt(i);
            d.at(i) = odtl->gas->density();
        }
    else {
        odtl->odtc->setGasStateAtPt(ipt);
        d.at(ipt) = odtl->gas->density();
    }
}

