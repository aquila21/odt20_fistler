/**
 * @file lv_uvw.cc
 * Header file for class lv_uvw
 */


#include "lv_uvw.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! lv_uvw  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_uvw::lv_uvw(odtline  *line,
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
 */

void lv_uvw::getRhsSrc(const int ipt){

    if(!L_transported)
        return;

    rhsSrc.resize(odtl->ngrd, 0.0);

    //-------------------------

    if(ipt==-1) {
        if(var_name == "uvel" && odtl->odtp->cCoord != 3.0) {
            rhsSrc = vector<double>(odtl->ngrd, -odtl->odtp->dPdx);
            for(int i=0; i<odtl->ngrd; i++) {
                if(odtl->odtp->Lbuoyant)
                    rhsSrc.at(i) += (odtl->rho->d.at(i) - odtl->rho->d.at(odtl->ngrd-1))*odtl->odtp->g;
                rhsSrc.at(i) /= odtl->rho->d.at(i);
            }
        }

        if(odtl->odtp->Lspatial)
            for(int i=0; i<odtl->ngrd; i++)
                rhsSrc.at(i) /= odtl->uvel->d.at(i);
    }
    else {
        if(var_name == "uvel" && odtl->odtp->cCoord != 3.0) {
            rhsSrc.at(ipt) = -odtl->odtp->dPdx;
            if(odtl->odtp->Lbuoyant)
                rhsSrc.at(ipt) += (odtl->rho->d.at(ipt) - odtl->rho->d.at(odtl->ngrd-1))*odtl->odtp->g;
            rhsSrc.at(ipt) /= odtl->rho->d.at(ipt);
        }

        if(odtl->odtp->Lspatial)
            rhsSrc.at(ipt) /= odtl->uvel->d.at(ipt);
    }
}

////////////////////////////////////////////////////////////////////////////////
/*! lv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 */

void lv_uvw::getRhsMix(const vector<double> &gf,
                       const vector<double> &dxc){

    rhsMix.resize(odtl->ngrd, 0.0);

    //------------------ Set fluxes

    flux.resize(odtl->ngrdf);
    vector<double> dvisc_f(odtl->ngrdf);

    interpVarToFacesHarmonic(odtl->dvisc->d, dvisc_f);

    //---------- Interior faces

    for (int i=1, im=0; i < odtl->ngrd; i++, im++)
        flux.at(i) = -gf.at(i) * dvisc_f.at(i)*(d.at(i) - d.at(im));

    //---------- Boundary faces

    if(odtl->odtp->bcType=="OUTFLOW") {
        flux.at(0) = 0.0;
        flux.at(odtl->ngrd) = 0.0;
    }
    else if(odtl->odtp->bcType=="WALL") {
        double bclo = var_name=="uvel" ? odtl->odtp->uBClo : var_name=="vvel" ? odtl->odtp->vBClo : odtl->odtp->wBClo;
        double bchi = var_name=="uvel" ? odtl->odtp->uBChi : var_name=="vvel" ? odtl->odtp->vBChi : odtl->odtp->wBChi;
        flux.at(0)          = -gf.at(0)          * dvisc_f.at(0)          * (d.at(0) - bclo);
        flux.at(odtl->ngrd) = -gf.at(odtl->ngrd) * dvisc_f.at(odtl->ngrd) * (bchi - d.at(odtl->ngrd-1));
    }
    else if(odtl->odtp->bcType=="WALL_OUT") {
        double bclo = var_name=="uvel" ? odtl->odtp->uBClo : var_name=="vvel" ? odtl->odtp->vBClo : odtl->odtp->wBClo;
        flux.at(0)          = -gf.at(0) * dvisc_f.at(0) * (d.at(0) - bclo);
        flux.at(odtl->ngrd) = 0.0;
    }
    else if(odtl->odtp->bcType=="PERIODIC") {
		if(odtl->odtp->probType=="SHEARFLOW"){
			double bclo = var_name=="uvel" ? odtl->odtp->uBClo : var_name=="vvel" ? odtl->odtp->vBClo : odtl->odtp->wBClo;
        	double bchi = var_name=="uvel" ? odtl->odtp->uBChi : var_name=="vvel" ? odtl->odtp->vBChi : odtl->odtp->wBChi;
        	flux.at(0)          = -gf.at(0)          * dvisc_f.at(0)          * (d.at(0) - bclo);
        	flux.at(odtl->ngrd) = -gf.at(odtl->ngrd) * dvisc_f.at(odtl->ngrd) * (bchi - d.at(odtl->ngrd-1));
		}
		else{
			int im = odtl->ngrd - 1;
            flux.at(0)          = -gf.at(0) * dvisc_f.at(0) * (d.at(0) - d.at(im));
            flux.at(odtl->ngrd) = flux.at(0);
		}
    }
    else {
        *odtl->io->ostrm << endl << "ERROR: bcType not recognized in lv_uvw::getRhsMix" << endl;
        exit(0);
    }

    //------------------ Compute the mixing term

    for(int i=0,ip=1; i<odtl->ngrd; i++, ip++){
       rhsMix.at(i) = -odtl->odtp->cCoord / (odtl->rho->d.at(i) * dxc.at(i)) *
                    (flux.at(ip) * pow(abs(odtl->posf->d.at(ip)), odtl->odtp->cCoord-1) -
                     flux.at(i)  * pow(abs(odtl->posf->d.at(i) ), odtl->odtp->cCoord-1));
       if(rhsMix.at(i) != rhsMix.at(i)){
                    odtl->io->writeDataFile("odt_error.dat", 0.0);
                    exit(0);
       }
    }
    if(odtl->odtp->Lspatial)
        for(int i=0; i<odtl->ngrd; i++)
            rhsMix.at(i) /= odtl->uvel->d.at(i);

}


