/**
 * @file lv_sdiff.cc
 * Header file for class lv_sdiff
 */


#include "lv_sdiff.h"
#include "odtline.h"

////////////////////////////////////////////////////////////////////////////////
/*! lv_sdiff  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_sdiff::lv_sdiff(odtline    *line,
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
/*! lv_sdiff setVar function
 *  @param ipt \input optional point to compute at
 */

void lv_sdiff::setVar(const int ipt){

    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be -1" << endl; 
        exit(0);
    }

    d.resize(odtl->ngrd, odtl->odtp->sdiff0);
    for(int i=0; i<odtl->ngrd; i++) {
        odtl->odtc->setGasStateAtPt(i);
        d.at(i) = odtl->tran->thermalConductivity(); //mk - WORK STATUS: container for scalar diffusivity with c_p = 1
    }
}

