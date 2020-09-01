/**
 * @file lv_mixf.cc
 * Header file for class lv_mixf
 */


#include "lv_mixf.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! lv_mixf  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_mixf::lv_mixf(odtline  *line,
                 const      string s,
                 const bool Lt,
                 const bool Lo) {

    odtl          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(odtl->ngrd, 0.0);

    Dmf           = odtl->io->streamProps["Dmf"] ? odtl->io->streamProps["Dmf"].as<double>() : 0.0;
    if(L_transported && Dmf == 0.0) {
        cout << endl << "ERROR: if you are transporting mixture fraction, you need to set Dmf";
        exit(0);
    }

}

////////////////////////////////////////////////////////////////////////////////
/*! Set mixture fraction from the gas state
 *  @param ipt \input optional point to compute at
 */

void lv_mixf::setVar(const int ipt){

    if(L_transported)
        return;        // don't set mixf from other quantities if its transported

    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be -1" << endl;
        exit(0);
    }

    d.resize(odtl->ngrd);

    int nsp = odtl->gas->nSpecies();
    vector<double> yi(nsp);
    for(int i=0; i<odtl->ngrd; i++) {
        for(int k=0; k<nsp; k++)
            yi.at(k) = odtl->ysp[k]->d.at(i);
        d.at(i) = odtl->strm->getMixtureFraction(&yi.at(0));
    }
}

////////////////////////////////////////////////////////////////////////////////
/*! lv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 *  todo: add in pressure term: unsteady, and nonuniform.
 */

void lv_mixf::getRhsSrc(const int ipt){

    if(!L_transported)
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

void lv_mixf::getRhsMix(const vector<double> &gf,
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

void lv_mixf::setFlux(const vector<double> &gf,
                      const vector<double> &dxc){

    flux.resize(odtl->ngrdf);
    vector<double> rhoD(odtl->ngrd, Dmf);
    vector<double> rhoD_f(odtl->ngrdf);

    for(int i=0; i<odtl->ngrd; i++)
        rhoD.at(i) *= odtl->rho->d.at(i);
    for(int i=0; i<odtl->ngrdf; i++)
        rhoD_f.at(i) = linearInterpToFace(i, rhoD);

    //==========================

    //---------- Interior faces

    for (int i=1, im=0, sum=0.0; i < odtl->ngrd; i++, im++)
        flux.at(i) = -rhoD_f.at(i)*gf.at(i)*(d.at(i) - d.at(im));

    //---------- Boundary faces

    flux.at(0)          = 0.0;
    flux.at(odtl->ngrd) = 0.0;
}

