/**
 * @file lv_ygas.cc
 * Header file for class lv_ygas
 */


#include "lv_ygas.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
// Setup static members

int lv_ygas::nspc;

///////////////////////////////////////////////////////////////////////////////
// Declare the global function prototype so it can be used in this source file

void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rr);

////////////////////////////////////////////////////////////////////////////////
/*! lv_ygas  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_ygas::lv_ygas(odtline  *line,
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

    aP   = odtl->io->params["aP"]   ? odtl->io->params["aP"].as<double>()       : 1.0;
    //aP_x = odtl->io->params["aP_x"] ? odtl->io->params["aP_x"].as<double>()/4.6 : 1.0E-10;
    aP_x = odtl->io->params["aP_x"] ? odtl->io->params["aP_x"].as<double>() : 1.0E-10;

}

////////////////////////////////////////////////////////////////////////////////
/*! lv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 * Gas temperature needs to be set to use problem specific RR
 */

void lv_ygas::getRhsSrc(const int ipt) {

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

    if(kMe==0) {                 // to save cost, compute needed terms for all lv_ygas objects using this one.

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
        rhsSrc.at(i) = rrSpc.at(kMe).at(i) *
             (odtl->mimx->time < aP_x ? aP : 1.0 ); //doldb this line
             //(1.0+(aP-1.0)*exp(-odtl->mimx->time/aP_x)); //doldb this line

    if(odtl->odtp->Lspatial)
        for(int i=iS; i<=iE; i++)
            rhsSrc.at(i) /= odtl->uvel->d.at(i);
}

////////////////////////////////////////////////////////////////////////////////
/*! lv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 * Note: temperature variable should be set already.
 * Note: species mass fluxes should be set already.
 */

void lv_ygas::getRhsMix(const vector<double> &gf,
                        const vector<double> &dxc){

    if(!L_transported) return;

    rhsMix.resize(odtl->ngrd, 0.0);

    setFlux(gf, dxc);

    //------------------ Compute the mixing term

    for(int i=0,ip=1; i<odtl->ngrd; i++, ip++)
       rhsMix.at(i) = -odtl->odtp->cCoord / (odtl->rho->d.at(i) * dxc.at(i)) *
                    (flux.at(ip) * pow(abs(odtl->posf->d.at(ip)), odtl->odtp->cCoord-1) -
                     flux.at(i)  * pow(abs(odtl->posf->d.at(i) ), odtl->odtp->cCoord-1));

    if(odtl->odtp->Lspatial)
        for(int i=0; i<odtl->ngrd; i++)
            rhsMix.at(i) /= odtl->uvel->d.at(i);

}

////////////////////////////////////////////////////////////////////////////////
/*! lv set face fluxes
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 */

void lv_ygas::setFlux(const vector<double> &gf,
                      const vector<double> &dxc){

    static vector<vector<double> > rhoD(nspc);
    static vector<vector<double> > rhoD_f(nspc);
    static vector<vector<double> > rhoDYinvM(nspc);
    static vector<vector<double> > rhoDYinvM_f(nspc);
    static vector<vector<double> > ysp_f(nspc);
    static vector<double>          Di(nspc);
    static vector<double>          dMdx;
    static vector<double>          MMw;

    if(kMe==0) {                 // to save cost, compute needed terms for all lv_ygas objects using this one.

        MMw.resize(odtl->ngrd);
        dMdx.resize(odtl->ngrdf);
        for(int k=0; k<nspc; k++) {
            rhoD.at(k).resize(odtl->ngrd);
            rhoD_f.at(k).resize(odtl->ngrdf);
            rhoDYinvM.at(k).resize(odtl->ngrd);
            rhoDYinvM_f.at(k).resize(odtl->ngrdf);
            ysp_f.at(k).resize(odtl->ngrdf);
            odtl->ysp[k]->flux.resize(odtl->ngrdf);
        }
        for(int i=0; i<odtl->ngrd; i++) {
            odtl->odtc->setGasStateAtPt(i);
            MMw.at(i) = odtl->gas->meanMolecularWeight();
            odtl->tran->getMixDiffCoeffs(&Di.at(0));
            for (int k=0; k<nspc; k++) {
                rhoD.at(k).at(i)      = odtl->rho->d.at(i)*Di.at(k);
                rhoDYinvM.at(k).at(i) = rhoD.at(k).at(i)*odtl->ysp[k]->d.at(i)/MMw.at(i);
            }
        }
        //for (int k=0; k<nspc; k++) {      // this breaks sometimes due to division issues, use linear instead (below)
        //    interpVarToFacesHarmonic(rhoD.at(k),      rhoD_f.at(k));
        //    interpVarToFacesHarmonic(rhoDYinvM.at(k), rhoDYinvM_f.at(k));
        //}

        for (int k=0; k<nspc; k++) {         // linear interpolation.
            for(int i=0; i<odtl->ngrdf; i++){
                rhoD_f.at(k).at(i)      = linearInterpToFace(i, rhoD.at(k));
                rhoDYinvM_f.at(k).at(i) = linearInterpToFace(i, rhoDYinvM.at(k));
                ysp_f.at(k).at(i)       = linearInterpToFace(i, odtl->ysp[k]->d);
            }
        }

        dMdx.at(0)          = 0.0;
        dMdx.at(odtl->ngrd) = 0.0;
        for (int i=1, im=0; i < odtl->ngrd; i++, im++)
            dMdx.at(i) = gf.at(i) * (MMw.at(i) - MMw.at(im));

        //==========================

        //---------- Interior faces

        double jstar;             // correction flux so that all fluxes sum to zero. This is equal to using a correction velocity
                                  // j_i_corrected = j_i - Yi*jstar; jstar = sum(j_i).
                                  // the previous approch just made j_last = -sum(j_all_but_last), where N2 was last.
        for (int i=1, im=0; i < odtl->ngrd; i++, im++) {
            jstar = 0.0;
            for(int k=0; k<nspc; k++) {
                odtl->ysp[k]->flux.at(i) = -rhoD_f.at(k).at(i)*gf.at(i)*(odtl->ysp[k]->d.at(i) - odtl->ysp[k]->d.at(im))
                    -rhoDYinvM_f.at(k).at(i)*dMdx.at(i);
                jstar += odtl->ysp[k]->flux.at(i);
            }
            for(int k=0; k<nspc; k++)
                odtl->ysp[k]->flux.at(i) -= ysp_f.at(k).at(i) * jstar;
        }

        //---------- Boundary faces

        for(int k=0; k<nspc; k++) {
            odtl->ysp[k]->flux.at(0)          = 0.0;                          // for wall or outflow; todo: wall flame specifics
            odtl->ysp[k]->flux.at(odtl->ngrd) = 0.0;
        }
    }
}

