/**
 * @file lv_dvisc_const.cc
 * Header file for class lv_dvisc_const
 */


#include "lv_dvisc_const.h"
#include "odtline.h"

////////////////////////////////////////////////////////////////////////////////
/*! lv_dvisc_const  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_dvisc_const::lv_dvisc_const(odtline    *line,
                               const      string s,
                               const bool Lt,
                               const bool Lo) {

    odtl          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(odtl->ngrd, odtl->odtp->kvisc0 * odtl->odtp->rho0);

    if(Lt){
        *odtl->io->ostrm << endl << "ERROR, you set dvisc to be transported. Resetting L_transported to false" << endl;
        L_transported = false;
    }

}

////////////////////////////////////////////////////////////////////////////////
/*! merger2cells function
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

void lv_dvisc_const::merge2cells(const int    imrg,
                                 const double m1,
                                 const double m2,
                                 const bool   LconstVolume) {

    d.erase(d.begin() + imrg+1);

}

////////////////////////////////////////////////////////////////////////////////
/*! lv_dvisc_const setVar function
 *  @param ipt \input optional point to compute at
 */

void lv_dvisc_const::setVar(const int ipt){
    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be -1" << endl;
        exit(0);
    }
    d.resize(odtl->ngrd, odtl->odtp->kvisc0 * odtl->odtp->rho0);
    for(int i=0;i<odtl->ngrd; i++ )
        d.at(i) = odtl->odtp->kvisc0 * odtl->odtp->rho0;
}

