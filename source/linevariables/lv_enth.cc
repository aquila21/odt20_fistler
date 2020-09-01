/**
 * @file lv_enth.cc
 * Header file for class lv_enth
 */


#include "lv_enth.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! lv_enth  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_enth::lv_enth(odtline    *line,
                 const      string s,
                 const bool Lt,
                 const bool Lo) {

    odtl          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(odtl->ngrd, 0.0);
    nspc          = odtl->gas->nSpecies();

    rad = new radiation();
    rad->init(odtl);

}

////////////////////////////////////////////////////////////////////////////////
/*! lv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 *  todo: add in pressure term: unsteady, and nonuniform.
 */

void lv_enth::getRhsSrc(const int ipt){

    if(!L_transported)
        return;

    if(LagSrc)
        return;

    rhsSrc.resize(odtl->ngrd, 0.0);

    int iS, iE;
    if(ipt==-1) {
        iS = 0;
        iE = odtl->ngrd-1;
    }
    else {
        iS = ipt;
        iE = ipt;
    }

    //----------- Compute the radiative source

    if(odtl->odtp->radType!="NONE") {

        vector<vector<double> > xMoleSp(odtl->ngrd, vector<double>(nspc));

        for(int i=0; i<odtl->ngrd; i++){
            odtl->odtc->setGasStateAtPt(i);
            odtl->gas->getMoleFractions(&xMoleSp.at(i).at(0));
        }
        rad->getRadHeatSource(xMoleSp, odtl->temp->d, rhsSrc);
    }

    if(odtl->odtp->Lspatial)
        for(int i=iS; i<=iE; i++)
            rhsSrc.at(i) /= odtl->uvel->d.at(i);

    LagSrc = true;      // this is reset to false in setCaseSpecificVars


}

////////////////////////////////////////////////////////////////////////////////
/*! lv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 * Note: temperature variable should be set already.
 * Note: species mass fluxes should be set already.
 */

void lv_enth::getRhsMix(const vector<double> &gf,
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

void lv_enth::setFlux(const vector<double> &gf,
                      const vector<double> &dxc){


    flux.resize(odtl->ngrdf, 0.0);

    vector<double>          tcond(odtl->ngrd, 0.0);
    vector<double>          tcond_f(odtl->ngrdf, 0.0);
    vector<vector<double> > hsp(nspc, vector<double>(odtl->ngrd));
    vector<double>          hh(nspc);
    double                  GasConstant = 8314.47215;             // J/kmol*K

    for(int i=0; i<odtl->ngrd; i++) {
        odtl->odtc->setGasStateAtPt(i);
        tcond.at(i) = odtl->tran->thermalConductivity();    // W/m*K
        odtl->gas->getEnthalpy_RT(&hh.at(0));               // non-dimensional enthalpy
        for (int k=0; k<nspc; k++)
            hsp.at(k).at(i) = hh.at(k) * odtl->temp->d.at(i) * GasConstant / odtl->gas->molecularWeight(k);        // J/kg
    }

    interpVarToFacesHarmonic(tcond, tcond_f);

    //========== Do the thermal conduction portion of the flux

    //---------- Interior faces

    for (int i=1, im=0; i < odtl->ngrd; i++, im++)
        flux.at(i) = -gf.at(i) * tcond_f.at(i)*(odtl->temp->d.at(i) - odtl->temp->d.at(im));

    //---------- Boundary faces

    if(odtl->odtp->bcType=="OUTFLOW") {
        flux.at(0)          = 0.0;
        flux.at(odtl->ngrd) = 0.0;
    }
    else if(odtl->odtp->bcType=="WALL") {
        if(odtl->odtp->hWallBCtype=="ADIABATIC") {
            flux.at(0)          = 0.0;
            flux.at(odtl->ngrd) = 0.0;
        }
        else if(odtl->odtp->hWallBCtype=="ISOTHERMAL") {
            flux.at(0)          = -gf.at(0)          * tcond_f.at(0)          * (odtl->temp->d.at(0)  - odtl->odtp->TBClo);
            flux.at(odtl->ngrd) = -gf.at(odtl->ngrd) * tcond_f.at(odtl->ngrd) * (odtl->odtp->TBChi - odtl->temp->d.at(odtl->ngrd-1));
        }
        else {
            cout << endl << "ERROR: hWallBCtype unknown" << endl;
            exit(0);
        }
    }
    else if(odtl->odtp->bcType=="WALL_OUT") {
        if(odtl->odtp->hWallBCtype=="ADIABATIC") {
            flux.at(0)          = 0.0;
            flux.at(odtl->ngrd) = 0.0;
        }
        else if(odtl->odtp->hWallBCtype=="ISOTHERMAL") {
            flux.at(0)          = -gf.at(0)          * tcond_f.at(0)          * (odtl->temp->d.at(0) - odtl->odtp->TBClo);
            flux.at(odtl->ngrd) = 0.0;
        }
        else {
            cout << endl << "ERROR: hWallBCtype unknown" << endl;
            exit(0);
        }
    }
    else if(odtl->odtp->bcType=="PERIODIC") {
        int im = odtl->ngrd - 1;
        flux.at(0)          = -gf.at(0)          * tcond_f.at(0)          * (odtl->temp->d.at(0) - odtl->temp->d.at(im));
        flux.at(odtl->ngrd) = flux.at(0);
    }
    else {
        *odtl->io->ostrm << endl << "ERROR: bcType not recognized in lv_enth::setFlux" << endl;
        exit(0);
    }

    //========== Add in the species flux portion.


    double h_f;
    double hjsum;
    for(int i=0; i<odtl->ngrdf; i++){
        hjsum = 0.0;
        for(int k=0; k<nspc; k++) {
            h_f   =  linearInterpToFace(i, hsp.at(k));
            hjsum += h_f * odtl->ysp[k]->flux.at(i);
        }
        flux.at(i) += hjsum;
    }

}

