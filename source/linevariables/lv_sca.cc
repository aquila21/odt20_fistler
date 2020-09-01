/**
 * @file lv_sca.cc
 * Header file for class lv_sca
 */


#include "lv_sca.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! lv_sca  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_sca::lv_sca(odtline  *line,
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
/*! Set the passive scalar
 *  @param ipt \input optional point to set at
 *
 * NOTE: The passive scalar evolves by mixing and stirring through the
 *  velocity field. The scalar distribution is a result of the flow and
 *  therefore does not require a setVar() method. Future changes, e.g.
 *  for the point-wise array manipulation, may be inserted here.
 */

//lv_sca::setVar(const int ipt) {}

////////////////////////////////////////////////////////////////////////////////
/*! lv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 */

void lv_sca::getRhsSrc(const int ipt){

    if(!L_transported)
        return;

    rhsSrc.resize(odtl->ngrd, 0.0);

    //-------------------------

    // Volume forcing like in Kawamura's heated channel would go in here.

}

////////////////////////////////////////////////////////////////////////////////
/*! lv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 */

void lv_sca::getRhsMix(const vector<double> &gf,
                       const vector<double> &dxc){

    rhsMix.resize(odtl->ngrd, 0.0);

    //------------------ Set fluxes

    flux.resize(odtl->ngrdf);
    vector<double> rho_f(odtl->ngrdf);

    interpVarToFacesHarmonic(odtl->rho->d, rho_f);

    //---------- Interior faces

    for (int i=1, im=0; i<odtl->ngrd; i++, im++)
        flux.at(i) = -gf.at(i) * rho_f.at(i)*odtl->odtp->sdiff0 * (d.at(i) - d.at(im)); // rho*sdiff0 is the "dynamic" diffusion coefficient

    //---------- Boundary faces

    if(odtl->odtp->bcType=="OUTFLOW") {
        flux.at(0) = 0.0;
        flux.at(odtl->ngrd) = 0.0;
    }
    else if(odtl->odtp->bcType=="WALL") {
        double bclo = odtl->odtp->sBClo;
        double bchi = odtl->odtp->sBChi;
        flux.at(0)          = -gf.at(0)          * rho_f.at(0)          * odtl->odtp->sdiff0          * (d.at(0) - bclo);
        flux.at(odtl->ngrd) = -gf.at(odtl->ngrd) * rho_f.at(odtl->ngrd) * odtl->odtp->sdiff0          * (bchi    - d.at(odtl->ngrd-1));
    }
    else if(odtl->odtp->bcType=="WALL_OUT") {
        double bclo = odtl->odtp->sBClo;
        flux.at(0)          = -gf.at(0) * rho_f.at(0) * odtl->odtp->sdiff0 * (d.at(0) - bclo);
        flux.at(odtl->ngrd) = 0.0;
    }
    else if(odtl->odtp->bcType=="PERIODIC") {
        flux.at(0)          = -gf.at(0) * rho_f.at(0) * odtl->odtp->sdiff0 * (d.at(0) - d.at(odtl->ngrd-1));
        flux.at(odtl->ngrd) = flux.at(0);
    }
    else {
        *odtl->io->ostrm << endl << "ERROR: bcType not recognized in lv_sca::getRhsMix" << endl;
        exit(0);
    }

    //------------------ Compute the mixing term

    for(int i=0,ip=1; i<odtl->ngrd; i++, ip++)
       rhsMix.at(i) = -odtl->odtp->cCoord / (odtl->rho->d.at(i) * dxc.at(i)) *
                    (flux.at(ip) * pow(abs(odtl->posf->d.at(ip)), odtl->odtp->cCoord-1) -
                     flux.at(i)  * pow(abs(odtl->posf->d.at(i) ), odtl->odtp->cCoord-1));

    if(odtl->odtp->Lspatial)
        for(int i=0; i<odtl->ngrd; i++)
            rhsMix.at(i) /= odtl->uvel->d.at(i);

}

