/**
 * @file lv_ygas_hips.cc
 * Header file for class lv_ygas_hips
 */


#include "lv_ygas_hips.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
// Setup static members

int lv_ygas_hips::nspc;

///////////////////////////////////////////////////////////////////////////////
// Declare the global function prototype so it can be used in this source file

void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rr);

////////////////////////////////////////////////////////////////////////////////
/** lv_ygas  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_ygas_hips::lv_ygas_hips(odtline    *line,
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

    nspc = odtl->gas->nSpecies();
    string spName(var_name.begin()+2, var_name.end());        // var_name is like y_O2. Just want the O2 part.
    kMe = odtl->gas->speciesIndex(spName);

}

////////////////////////////////////////////////////////////////////////////////
/** Source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 *  Gas temperature/rho needs to be set to use problem specific RR
 */

void lv_ygas_hips::getRhsSrc(const int ipt) {

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

    if(kMe==0) {                 // to save cost, compute needed terms for all lv_ygas_hips objects using this one.

        for(int i=iS; i<=iE; i++) {
#ifdef PROBLEMSPECIFICRR
            // rho and T should already be set (make sure).
            for(int k=0; k<nspc; k++)
                yi.at(k) = odtl->ysp[k]->d.at(i);
            getProblemSpecificRR(odtl->rho->d.at(i), odtl->temp->d.at(i), odtl->odtp->pres, &yi.at(0), &rr.at(0));
#else
            odtl->odtc->setGasStateAtPt(i);
            odtl->gas->getNetProductionRates(&rr.at(0));
#endif
            for (int k=0; k<nspc; k++)
                rrSpc.at(k).at(i) = rr.at(k) * odtl->gas->molecularWeight(k) / odtl->rho->d.at(i);   // kmol/(m³ s)*(kg/kmol)*(kg/m3) = 1/s
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

void lv_ygas_hips::getRhsMix(const vector<double> &gf,
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

