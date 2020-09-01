/**
 * @file lv_dvisc.cc
 * Header file for class lv_dvisc
 */


#include "lv_dvisc.h"
#include "odtline.h"

////////////////////////////////////////////////////////////////////////////////
/*! lv_dvisc  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_dvisc::lv_dvisc(odtline    *line,
                   const      string s,
                   const bool Lt,
                   const bool Lo) {

    odtl          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(odtl->ngrd, odtl->odtp->kvisc0 * odtl->odtp->rho0);

    if(Lt){
        *odtl->io->ostrm << endl << "WARNING, you set dvisc to be transported. Resetting L_transported to false" << endl;
        L_transported = false;
    }

}

////////////////////////////////////////////////////////////////////////////////
/*! lv_dvisc setVar function
 *  @param ipt \input optional point to compute at
 */

void lv_dvisc::setVar(const int ipt){

    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be -1" << endl;
        exit(0);
    }

    d.resize(odtl->ngrd, odtl->odtp->kvisc0 * odtl->odtp->rho0);
    for(int i=0; i<odtl->ngrd; i++) {
        odtl->odtc->setGasStateAtPt(i);
        d.at(i) = odtl->tran->viscosity();
    }
}

