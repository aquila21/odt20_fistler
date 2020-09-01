/**
 * @file lv_ygas_flmlt.cc
 * Header file for class lv_ygas_flmlt
 */


#include "lv_ygas_flmlt.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
// Setup static members

int lv_ygas_flmlt::nspc;

///////////////////////////////////////////////////////////////////////////////
// Declare the global function prototype so it can be used in this source file

void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rr);

////////////////////////////////////////////////////////////////////////////////
/*! lv_ygas  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_ygas_flmlt::lv_ygas_flmlt(odtline  *line,
                 const      string s,
                 const bool Lt,
                 const bool Lo) {

    odtl               = line;
    var_name           = s;
    L_transported      = Lt;
    L_output           = Lo;
    d                  = vector<double>(odtl->ngrd, 0.0);

    nspc = odtl->gas->nSpecies();

    string spName(var_name.begin()+2, var_name.end());        // var_name is like y_O2. Just want the O2 part.
    kMe                = odtl->gas->speciesIndex(spName);

}

////////////////////////////////////////////////////////////////////////////////
/*! lv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 * Gas temperature needs to be set to use problem specific RR
 */

void lv_ygas_flmlt::getRhsSrc(const int ipt) {

    if(!L_transported)
        return;

    rhsSrc.resize(odtl->ngrd, 0.0);

    static vector<vector<double> > rrSpc(nspc);    // [nspc][ngrd]
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

    if(kMe==0) {                 // to save cost, compute needed terms for all lv_ygas_flmlt objects using this one.

        for(int k=0; k<nspc; k++)
            rrSpc.at(k).resize(odtl->ngrd);

        for(int i=iS; i<=iE; i++) {
#ifdef PROBLEMSPECIFICRR
            // make sure rho and T are set first (though it should be for the diffuser at least).
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
/*! lv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 * Note: temperature variable should be set already.
 *
 * We are using a finite difference formulation here (different that the "ODT" formulation).
 */

void lv_ygas_flmlt::getRhsMix(const vector<double> &gf,
                        const vector<double> &dxc){

    if(!L_transported) return;

    rhsMix.resize(odtl->ngrd, 0.0);

    //------------------ Compute the mixing term

    double dp, dm;
    double d2YdZ2;
    int i;

    //-------- left point

    i = 0;
    dp = odtl->pos->d.at(i+1) - odtl->pos->d.at(i);
    dm = odtl->pos->d.at(i)   - 0.0;
    d2YdZ2 = 2.0/(dp+dm)*( (d.at(i+1)-d.at(i)             )/dp -
                           (d.at(i)  -odtl->strm->y0[kMe] )/dm );
    rhsMix.at(i) = odtl->chi->d.at(i)/2.0 * d2YdZ2;

    //-------- interior points

    for(i=1; i<odtl->ngrd-1; i++) {
        dp = odtl->pos->d.at(i+1) - odtl->pos->d.at(i);
        dm = odtl->pos->d.at(i)   - odtl->pos->d.at(i-1);
        d2YdZ2 = 2.0/(dp+dm)*( (d.at(i+1) - d.at(i)  )/dp -
                               (d.at(i)   - d.at(i-1))/dm );
        rhsMix.at(i) = odtl->chi->d.at(i)/2.0 * d2YdZ2;
    }

    //-------- right point

    i = odtl->ngrd-1;
    dp = 1.0                  - odtl->pos->d.at(i);
    dm = odtl->pos->d.at(i)   - odtl->pos->d.at(i-1);
    d2YdZ2 = 2.0/(dp+dm)*( (odtl->strm->y1[kMe] - d.at(i)  )/dp -
                           (d.at(i)             - d.at(i-1))/dm );
    rhsMix.at(i) = odtl->chi->d.at(i)/2.0 * d2YdZ2;

}

