/**
 * @file lv_soot_flmlt.cc
 * Header file for class lv_soot_flmlt
 * @author Victoria B. Lansinger and David Lignell
 */

#include "lv_soot_flmlt.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

/*//////////////////////////////////////////////////////////////////////////////
 *! lv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 *
 * See page 191 of Lignell_thesis, equation 7.9.
 * Equation is dm/dt = C*dm/dZ + S.
 * This is hyperbolic and C is a "convective wind".
 * We'll upwind the discretization of dm/dZ depending on the sign of C.
 * Note, here S does not include the m_r factor and also the M_r/rho term
 *        that appears in the thesis. The M_r/rho is the reaction source term
 *        computed elsewhere. S without the m_r factor keeps S and C independent of m_r
 *        so that we can compute them once for kMe = 0. The m_r factor is then put
 *        in at the bottom.
 * As for ygas, we are just using a finite difference formulation here.
 */

void lv_soot_flmlt::getRhsMix(const vector<double> &gf,
                              const vector<double> &dxc){

    if(!L_transported) return;

    rhsMix.resize(odtl->ngrd, 0.0);

    //------------------

    static vector<double> D;
    static vector<double> mu;
    static vector<double> beta;
    static vector<double> C;
    static vector<double> S;
    double dp, dm;
    int i;

    //------------------

    if(kMe==0) {

        D.resize(odtl->ngrd);
        mu.resize(odtl->ngrd);
        beta.resize(odtl->ngrd);
        C.resize(odtl->ngrd);
        S.resize(odtl->ngrd);

        //------------------

        for(i=0; i<odtl->ngrd; i++){
            //odtl->odtc->setGasStateAtPt(i);   // keep commented for decoupled soot
            D[i] = odtl->tran->thermalConductivity() / odtl->rho->d[i] / odtl->gas->cp_mass();
            mu[i] = odtl->tran->viscosity();
            beta[i] = sqrt(odtl->chi->d[i]*0.5/D[i]);
        }

        //------------------ Set C and S

        double dTdZ, drDbdZ, dmbTdZ, d2TdZ2;

        //-------- left point

        i = 0;
        dp = odtl->pos->d.at(i+1) - odtl->pos->d.at(i);
        dm = odtl->pos->d.at(i)   - 0.0;
        dTdZ   = (odtl->temp->d[i] - odtl->strm->T0)/dm;
        drDbdZ = (odtl->rho->d[i]*D[i]*beta[i] - 0.0)/dm;
        dmbTdZ = (mu[i]*beta[i]/odtl->temp->d[i] - 0.0)/dm;
        d2TdZ2 = 2.0/(dp+dm)*( (odtl->temp->d.at(i+1) - odtl->temp->d.at(i))/dp -
                               (odtl->temp->d.at(i)   - odtl->strm->T0     )/dm );

        C[i] = 0.554*mu[i]/odtl->rho->d[i]/odtl->temp->d[i]*beta[i]*beta[i]*dTdZ - beta[i]/odtl->rho->d[i]*drDbdZ;
        S[i] = 0.554/odtl->rho->d[i]*( beta[i]*dTdZ*dmbTdZ + beta[i]*beta[i]*mu[i]/odtl->temp->d[i]*d2TdZ2 );

        //-------- interior points

        for(i=1; i<odtl->ngrd-1; i++) {
            dp = odtl->pos->d.at(i+1) - odtl->pos->d.at(i);
            dm = odtl->pos->d.at(i)   - odtl->pos->d.at(i-1);
            dTdZ   = (odtl->temp->d[i+1] - odtl->temp->d[i-1])/(dm+dp);
            drDbdZ = (odtl->rho->d[i+1]*D[i+1]*beta[i+1] - odtl->rho->d[i-1]*D[i-1]*beta[i-1])/(dm+dp);
            dmbTdZ = (mu[i+1]*beta[i+1]/odtl->temp->d[i+1] - mu[i-1]*beta[i-1]/odtl->temp->d[i-1])/(dm+dp);
            d2TdZ2 = 2.0/(dp+dm)*( (odtl->temp->d[i+1] - odtl->temp->d[i]  )/dp -
                                   (odtl->temp->d[i]   - odtl->temp->d[i-1])/dm );

            C[i] = 0.554*mu[i]/odtl->rho->d[i]/odtl->temp->d[i]*beta[i]*beta[i]*dTdZ - beta[i]/odtl->rho->d[i]*drDbdZ;
            S[i] = 0.554/odtl->rho->d[i]*( beta[i]*dTdZ*dmbTdZ + beta[i]*beta[i]*mu[i]/odtl->temp->d[i]*d2TdZ2 );
        }

        //-------- right point

        i = odtl->ngrd-1;
        dp = 1.0                  - odtl->pos->d.at(i);
        dm = odtl->pos->d.at(i)   - odtl->pos->d.at(i-1);
        dTdZ   = (odtl->strm->T1 - odtl->temp->d[i])/dp;
        drDbdZ = (0.0 - odtl->rho->d[i]*D[i]*beta[i])/dp;
        dmbTdZ = (0.0 - mu[i]*beta[i]/odtl->temp->d[i])/dp;
        d2TdZ2 = 2.0/(dp+dm)*( (odtl->strm->T1      - odtl->temp->d.at(i)  )/dp -
                               (odtl->temp->d.at(i) - odtl->temp->d.at(i-1))/dm );

        C[i] = 0.554*mu[i]/odtl->rho->d[i]/odtl->temp->d[i]*beta[i]*beta[i]*dTdZ - beta[i]/odtl->rho->d[i]*drDbdZ;
        S[i] = 0.554/odtl->rho->d[i]*( beta[i]*dTdZ*dmbTdZ + beta[i]*beta[i]*mu[i]/odtl->temp->d[i]*d2TdZ2 );

    }

    //------------------ Set mixing term
    // dm/dt = ( C*dm/dZ + m*S ) + RxnTerm
    // Upwind the C*dm/dZ part of the mixing term.

    double dmdZ;

    //-------- left point

    i = 0;
    dp = odtl->pos->d.at(i+1) - odtl->pos->d.at(i);
    dm = odtl->pos->d.at(i)   - 0.0;
    dmdZ = (C[i] < 0.0) ? (d[i]-0.0)/dm : (d[i+1]-d[i])/dp;
    rhsMix[i] = C[i]*dmdZ + d[i]*S[i];              // d*S puts in the m_r factor on S.

    //-------- interior points

    for(i=1; i<odtl->ngrd-1; i++) {
        dp = odtl->pos->d.at(i+1) - odtl->pos->d.at(i);
        dm = odtl->pos->d.at(i)   - odtl->pos->d.at(i-1);
        dmdZ = (C[i] < 0.0) ? (d[i]-d[i-1])/dm : (d[i+1]-d[i])/dp;
        rhsMix[i] = C[i]*dmdZ + d[i]*S[i];
    }

    //-------- right point

    i = odtl->ngrd-1;
    dp = 1.0                  - odtl->pos->d.at(i);
    dm = odtl->pos->d.at(i)   - odtl->pos->d.at(i-1);
    dmdZ = (C[i] < 0.0) ? (d[i]-d[i-1])/dm : (0.0-d[i])/dp;
    rhsMix[i] = C[i]*dmdZ + d[i]*S[i];

}












//void lv_soot_flmlt::getRhsMix(const vector<double> &gf,
//                        const vector<double> &dxc){
//
//    if(!L_transported) return;
//
//    rhsMix.resize(odtl->ngrd, 0.0);
//
//    //------------------ Compute the mixing term
//
//    double dp, dm;
//    double d2mdZ2;
//    int i;
//
//    //-------- left point
//
//    i = 0;
//    dp = odtl->pos->d.at(i+1) - odtl->pos->d.at(i);
//    dm = odtl->pos->d.at(i)   - 0.0;
//    d2mdZ2 = 2.0/(dp+dm)*( (d.at(i+1) - d.at(i) )/dp -
//                           (d.at(i)   - 0.0     )/dm );
//    rhsMix.at(i) = odtl->chi->d.at(i)/2.0 * d2mdZ2;
//
//    //-------- interior points
//
//    for(i=1; i<odtl->ngrd-1; i++) {
//        dp = odtl->pos->d.at(i+1) - odtl->pos->d.at(i);
//        dm = odtl->pos->d.at(i)   - odtl->pos->d.at(i-1);
//        d2mdZ2 = 2.0/(dp+dm)*( (d.at(i+1) - d.at(i)  )/dp -
//                               (d.at(i)   - d.at(i-1))/dm );
//        rhsMix.at(i) = odtl->chi->d.at(i)/2.0 * d2mdZ2;
//    }
//
//    //-------- right point
//
//    i = odtl->ngrd-1;
//    dp = 1.0                  - odtl->pos->d.at(i);
//    dm = odtl->pos->d.at(i)   - odtl->pos->d.at(i-1);
//    d2mdZ2 = 2.0/(dp+dm)*( (0.0                 - d.at(i)  )/dp -
//                           (d.at(i)             - d.at(i-1))/dm );
//    rhsMix.at(i) = odtl->chi->d.at(i)/2.0 * d2mdZ2;
//
//}

