/**
 * @file lv_aDL.cc
 * Header file for class lv_aDL
 */


#include "lv_aDL.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! lv_aDL  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_aDL::lv_aDL(odtline  *line,
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
/*! Set aDL
 *  @param ipt \input optional point to compute at
 */

void lv_aDL::setVar(const int ipt){

    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be -1" << endl;
        exit(0);
    }

    d.resize(odtl->ngrd,0.0);
}

