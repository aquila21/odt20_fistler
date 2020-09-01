/**
 * @file lv_ygas_noRxn.cc
 * Header file for class lv_ygas_noRxn
 */


#include "lv_ygas_noRxn.h"
#include "odtline.h"


////////////////////////////////////////////////////////////////////////////////
/*! lv source term part of the rhs function.
 *  @param ipt \input optional point to compute source at.
 */

void lv_ygas_noRxn::getRhsSrc(const int ipt){

    if(!L_transported)
        return;

    rhsSrc = vector<double>(odtl->ngrd, 0.0);

}

