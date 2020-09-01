/**
 * @file lv_sdiff_const.cc
 * Header file for class lv_sdiff_const
 */

#include "lv_sdiff_const.h"
#include "odtline.h"

////////////////////////////////////////////////////////////////////////////////
/*! lv_sdiff_const  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_sdiff_const::lv_sdiff_const(odtline    *line,
                               const      string s,
                               const bool Lt,
                               const bool Lo) {

    odtl          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(odtl->ngrd, odtl->odtp->sdiff0);

    if(Lt){
        *odtl->io->ostrm << endl << "WARNING, you set sdiff to be transported. Resetting L_transported to false" << endl;
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

void lv_sdiff_const::merge2cells(const int    imrg,
                                 const double m1,
                                 const double m2,
                                 const bool   LconstVolume) {

    d.erase(d.begin() + imrg+1);
   
}

////////////////////////////////////////////////////////////////////////////////
/*! lv_sdiff_const setVar function
 *  @param ipt \input optional point to compute at
 */

void lv_sdiff_const::setVar(const int ipt){
    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be = -1" << endl; 
        exit(0);
    }
    d.resize(odtl->ngrd, odtl->odtp->sdiff0);
}

