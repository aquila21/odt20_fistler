/**
 * @file lv_ygas_cold_hips.cc
 * Header file for class lv_ygas_cold_hips
 */


#include "lv_ygas_cold_hips.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
// Setup static members

int lv_ygas_cold_hips::nspc;

////////////////////////////////////////////////////////////////////////////////
/** lv_ygas  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_ygas_cold_hips::lv_ygas_cold_hips(odtline    *line,
                           const      string s,
                           const bool Lt,
                           const bool Lo) {

    odtl               = line;
    var_name           = s;
    L_transported      = Lt;
    L_output           = Lo;
    d                  = vector<double>(odtl->ngrd, 0.0);

    rhsSrc.resize(odtl->ngrd, 0.0);
    rhsMix.resize(odtl->ngrd, 0.0);

    nspc = 4;
    string spName(var_name.begin()+2, var_name.end());        // var_name is like y_O2. Just want the O2 part.

    if(spName=="A") kMe = 0;
    if(spName=="B") kMe = 1;
    if(spName=="R") kMe = 2;
    if(spName=="P") kMe = 3;

    Da1 = odtl->io->inputFile["rxnParams"]["Da1"].as<double>();
    Da2 = odtl->io->inputFile["rxnParams"]["Da2"].as<double>();

}

////////////////////////////////////////////////////////////////////////////////
/** Source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 *  Gas temperature/rho needs to be set to use problem specific RR
 */

void lv_ygas_cold_hips::getRhsSrc(const int ipt) {

    if(!L_transported)
        return;

    static vector<vector<double> > rrSpc(nspc, vector<double>(odtl->ngrd));    // [nspc][ngrd]
    static vector<double>          yi(nspc);       // [nspc]
    static vector<double>          rr(nspc);       // [nspc]

    int iS, iE;
    if(ipt==-1) {
        iS = 0;
        iE = odtl->ngrd-1;
    }
    else {
        iS = ipt;
        iE = ipt;
    }

    if(kMe==0) {                 // to save cost, compute needed terms for all lv_ygas_cold_hips objects using this one.

        for(int i=iS; i<=iE; i++) {

            double A = odtl->ysp[0]->d.at(i);
            double B = odtl->ysp[1]->d.at(i);
            double R = odtl->ysp[2]->d.at(i);
            double P = odtl->ysp[3]->d.at(i);

            rrSpc.at(0).at(i) = -Da1*A*B    - 0.5*Da2*A*R;
            rrSpc.at(1).at(i) = -Da1*A*B;
            rrSpc.at(2).at(i) = 2.0*Da1*A*B - Da2*A*R;
            rrSpc.at(3).at(i) = 1.5*Da2*A*R;
        }
    }

    for(int i=iS; i<=iE; i++)    
        rhsSrc.at(i) = rrSpc.at(kMe).at(i);
}

////////////////////////////////////////////////////////////////////////////////
/** lv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 *      These parameters are for the inherited interface (not used here).
 * Solving: dA/dt = (B-A)/τ,
 *          dB/dt = (A-B)/τ,
 * where A is Y_k of one parcel, and B is Y_k of its neighbor.
 * The solution is A = A_0 + (B_0-A_0)/2*(1-exp(-2*t/τ)),
 *                 B = B_0 + (A_0-B_0)/2*(1-exp(-2*t/τ)).
 * Note, we only need to solve for A at every point as long as we reference B.
 * Since we know the analytic solution, form the rate for use with Explicit Euler:
 * Rate_A = (A(t+Δt) - A(t)) / Δt = (B_0-A_0)/(2*Δt)*(1-exp(-2Δt/τ)).
 * Hence, Explicit Euler will give the exact solution.
 */

void lv_ygas_cold_hips::getRhsMix(const vector<double> &gf,
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

