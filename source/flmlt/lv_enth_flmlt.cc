/**
 * @file lv_enth_flmlt.cc
 * Header file for class lv_enth_flmlt
 */

#include "lv_enth_flmlt.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! lv_enth_flmlt  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_enth_flmlt::lv_enth_flmlt(odtline    *line,
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
/*! lv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 *  todo: add in pressure term: unsteady, and nonuniform.
 */

void lv_enth_flmlt::getRhsSrc(const int ipt){

    if(!L_transported)
        return;

    if(LagSrc)
        return;

    rhsSrc.resize(odtl->ngrd, 0.0);

}

////////////////////////////////////////////////////////////////////////////////
/*! lv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 * Note: temperature variable should be set already.
 * Note: species mass fluxes should be set already.
 */

void lv_enth_flmlt::getRhsMix(const vector<double> &gf,
                        const vector<double> &dxc){

    if(!L_transported)
        return;

    rhsMix.resize(odtl->ngrd, 0.0);
}

////////////////////////////////////////////////////////////////////////////////
/*! lv_enth_flmlt setVar function
 *  @param ipt \input optional point to compute at
 */

void lv_enth_flmlt::setVar(const int ipt){

    d.resize(odtl->ngrd);
    for(int i=0; i<odtl->ngrd; i++)
        d.at(i) = odtl->strm->h0 * (1.0-odtl->pos->d.at(i)) +
            odtl->strm->h1 * (odtl->pos->d.at(i));
}

