/**
 * @file lv_temp.cc
 * Header file for class lv_temp
 */


#include "lv_temp.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! lv_temp  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_temp::lv_temp(odtline  *line,
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
/*! Set temperature from the gas state
 *  @param ipt \input optional point to compute at
 */

void lv_temp::setVar(const int ipt){

    d.resize(odtl->ngrd);
    if(ipt == -1)
        for(int i=0; i<odtl->ngrd; i++) {
            odtl->odtc->setGasStateAtPt(i);
            d.at(i) = odtl->gas->temperature();
        }
    else {
        odtl->odtc->setGasStateAtPt(ipt);
        d.at(ipt) = odtl->gas->temperature();
    }
}

