/**
 * @file lv_rho_const.cc
 * Header file for class lv_rho_const
 */


#include "lv_rho_const.h"
#include "odtline.h"

////////////////////////////////////////////////////////////////////////////////
/*! lv_rho_const  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_rho_const::lv_rho_const(odtline    *line,
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
/*! lv_rho_const merger2cells function
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

void lv_rho_const::merge2cells(const int    imrg,
                               const double m1,
                               const double m2,
                               const bool   LconstVolume) {

    d.erase(d.begin() + imrg+1);

}

////////////////////////////////////////////////////////////////////////////////
/*! lv_rho_const setVar function
 *  @param ipt \input optional point to compute at
 */

void lv_rho_const::setVar(const int ipt){

    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be -1" << endl;
        exit(0);
    }
 
    d.resize(odtl->ngrd, odtl->odtp->rho0);
    for(int i=0; i < odtl->ngrd; i++) d.at(i) = odtl->odtp->rho0;
}

