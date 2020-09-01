/**
 * @file lv_enth_hips.cc
 * Header file for class lv_enth_hips
 */

#include "lv_enth_hips.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/** lv_enth_hips  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_enth_hips::lv_enth_hips(odtline      *line,
                           const string s,
                           const bool   Lt,
                           const bool   Lo) {

    odtl          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(odtl->ngrd, 0.0);

    rhsSrc.resize(odtl->ngrd, 0.0);
    rhsMix.resize(odtl->ngrd, 0.0);

}

////////////////////////////////////////////////////////////////////////////////
/** lv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 *      This parameter is for the inerhited interface.
 *  This function just returns since the enthalpy source is zero in hips
 *     (well, as currently defined).
 */

void lv_enth_hips::getRhsSrc(const int ipt){
        return;
}

////////////////////////////////////////////////////////////////////////////////
/** lv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 *      These parameters are for the inherited interface (not used here).
 * Solving: dA/dt = (B-A)/τ,
 *          dB/dt = (A-B)/τ,
 * where A is enthalpy of one parcel, and B is enthalpy of its neighbor.
 * The solution is A = A_0 + (B_0-A_0)/2*(1-exp(-2*t/τ)),
 *                 B = B_0 + (A_0-B_0)/2*(1-exp(-2*t/τ)).
 * Note, we only need to solve for A at every point as long as we reference B.
 * Since we know the analytic solution, form the rate for use with Explicit Euler:
 * Rate_A = (A(t+Δt) - A(t)) / Δt = (B_0-A_0)/(2*Δt)*(1-exp(-2Δt/τ)).
 * Hence, Explicit Euler will give the exact solution.
 */

void lv_enth_hips::getRhsMix(const vector<double> &gf,
                             const vector<double> &dxc){

    if(!L_transported)
        return;

    int ime;
    int inb;

    for(int i=0; i<odtl->ngrd; i++) {

        ime = odtl->solv->pLoc[i];
        inb = (odtl->solv->pLoc[i]%2 == 0) ? odtl->solv->pLoc[i+1] : odtl->solv->pLoc[i-1];

        rhsMix[ime] = (d.at(inb)-d.at(ime))/(2.0*odtl->mimx->dt) *
                      (1.0 - exp(-2.0*odtl->mimx->dt/odtl->solv->tMix));
    }

}

