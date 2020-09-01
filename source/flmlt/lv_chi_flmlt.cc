/**
 * @file lv_chi_flmlt.cc
 * Header file for class lv_chi_flmlt
 */


#include "lv_chi_flmlt.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! lv_chi_flmlt  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_chi_flmlt::lv_chi_flmlt(odtline  *line,
                 const               string s,
                 const bool          Lt,
                 const bool          Lo) {

    odtl          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(odtl->ngrd, 0.0);

    if(odtl->io->inputFile["flmlt"]["chiStoic"])
        chiStoic = odtl->io->inputFile["flmlt"]["chiStoic"].as<double>();
    else {
        cout << endl << "ERROR: missing chiStoic parameter in input file" << endl;
        exit(0);
    }
    chi0 = chiStoic / pow(1.0-pow(2.0*odtl->strm->mixfStoic-1.0, 2.0), 2.0);
    //chi0 = chiStoic; //doldb

}

////////////////////////////////////////////////////////////////////////////////
/** Set scalar dissipation rate (chi) \cond
 *  @param ipt \input optional point to compute at
 *  Chi = chi0 * (1-(2*mixf - 1)^2)^2
 */

void lv_chi_flmlt::setVar(const int ipt){

    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be < 1" << endl;
        exit(0);
    }

    d.resize(odtl->ngrd);

    for(int i=0; i<odtl->ngrd; i++)
        d.at(i) = chi0 * pow(1.0-pow(2.0*odtl->pos->d.at(i)-1.0, 2.0), 2.0);

}

