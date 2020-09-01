/**
 * @file lv_soot.cc
 * Header file for class lv_soot
 *
 * @author Victoria B. Lansinger
 */

#include "lv_soot.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
// Static members

int         lv_soot::N = 0;
soot_chem*     lv_soot::chem = new soot_chem();

////////////////////////////////////////////////////////////////////////////////
/*! lv_soot  constructor function
 *
 * @param
 * @param
 */

lv_soot::lv_soot(odtline   *line,
                 const      string s,
                 const bool Lt,
                 const bool Lo) {

    odtl               = line;
    var_name           = s;
    L_transported      = Lt;
    L_output           = Lo;
    d                  = vector<double>(odtl->ngrd, 0.0);

    kMe         = N++;                    ///< kMe

    if (kMe == 0) {
        chem->init(odtl);
    }

    nsvar = odtl->odtp->nsvar;

}

////////////////////////////////////////////////////////////////////////////////
/*! lv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 * Note: temperature variable should be set already.
 * Note: species mass fluxes should be set already.
 */

void lv_soot::getRhsMix(const vector<double> &gf,
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

void lv_soot::setFlux(const vector<double> &gf,
                      const vector<double> &dxc){

    int nVar = odtl->io->sootParams["nsvar"].as<int>();

    static vector<vector<double> > dvisc(nVar);
    static vector<vector<double> > dvisc_f(nVar);

    if(kMe==0) {                 // to save cost, compute needed terms for all lv_soot objects using this one.

        for(int k=0; k<nVar; k++) {
            dvisc.at(k).resize(odtl->ngrd);
            dvisc_f.at(k).resize(odtl->ngrdf);
            odtl->svar[k]->flux.resize(odtl->ngrdf);
        }
        for(int i=0; i<odtl->ngrd; i++) {
            // odtl->odtc->setGasStateAtPt(i);   // keep this commented for decoupled soot
            for (int k=0; k<nVar; k++) {
                dvisc.at(k).at(i) = odtl->dvisc->d.at(i);
            }
        }

        //for (int k=0; k<nVar; k++) {      // this breaks sometimes due to division issues, use linear instead (below)
        //    interpVarToFacesHarmonic(dvisc.at(k),      dvisc_f.at(k));
        //}

        for (int k=0; k<nVar; k++) {         // linear interpolation.
            for(int i=0; i<odtl->ngrdf; i++){
                dvisc_f.at(k).at(i)     = linearInterpToFace(i, dvisc.at(k));
            }
        }

        //==========================

        //---------- Interior faces

        for (int i=1, im=0; i < odtl->ngrd; i++, im++) {

            //-------------- Thermophoretic transport is "hyperbolic", flux looks like f=v*phi, so upwind it for stability

            for(int k=0; k<nVar; k++){
                double vel = -0.55415 * dvisc_f.at(k).at(i)*(log(odtl->temp->d.at(i))-log(odtl->temp->d.at(im)))*gf.at(i);
                double ii = (vel > 0) ? im : i;
                odtl->svar[k]->flux.at(i) = vel * odtl->svar[k]->d.at(ii);
            }
        }

        //---------- Boundary faces

        for(int k=0; k<nVar; k++) {
            odtl->svar[k]->flux.at(0)          = 0.0;                          // for wall or outflow; todo: wall flame specifics
            odtl->svar[k]->flux.at(odtl->ngrd) = 0.0;
        }
    }
}
