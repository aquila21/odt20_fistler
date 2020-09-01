/**
 * @file odtcase_jetMixlRxn.cc
 * Header file for class odtcase_jetMixlRxn
 */


#include "odtcase_jetMixlRxn.h"
#include "odtline.h"
#include "lv.h"
#include "lv_pos.h"
#include "lv_posf.h"
#include "lv_rho.h"
#include "lv_dvisc.h"
#include "lv_uvw.h"
#include "lv_enth.h"
#include "lv_temp.h"
#include "lv_ygas.h"
#include "lv_mixf.h"
#include "lv_chi.h"
#include "lv_hr.h"

#include <cmath>
#include <string>

////////////////////////////////////////////////////////////////////////////////
/** odtcase_jetMixlRxn initialization function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void odtcase_jetMixlRxn::init(odtline *p_odtl) {

    odtl = p_odtl;

    vector<double> gammas(4,0.0);
    gammas[0] = 2.0/odtl->gas->atomicWeight(odtl->gas->elementIndex("C"));
    gammas[1] = 0.5/odtl->gas->atomicWeight(odtl->gas->elementIndex("H"));
    gammas[2] = -1.0/odtl->gas->atomicWeight(odtl->gas->elementIndex("O"));
    gammas[3] = 0.0;
    odtl->strm->init(odtl, gammas);

    odtl->v.push_back(new lv_pos(   odtl, "pos",     false, true ));   // last are: L_transported, L_output
    odtl->v.push_back(new lv_posf(  odtl, "posf",    false, true ));
    odtl->v.push_back(new lv_rho(   odtl, "rho",     false, true ));
    odtl->v.push_back(new lv_dvisc( odtl, "dvisc",   false, true ));
    odtl->v.push_back(new lv_uvw(   odtl, "uvel",    true,  true ));
    odtl->v.push_back(new lv_uvw(   odtl, "vvel",    true,  true ));
    odtl->v.push_back(new lv_uvw(   odtl, "wvel",    true,  true ));
    odtl->v.push_back(new lv_temp(  odtl, "temp",    false, true ));
    odtl->v.push_back(new lv_mixf(  odtl, "mixf",    false, true ));
    odtl->v.push_back(new lv_chi(   odtl, "chi",     false, true ));
    odtl->v.push_back(new lv_hr(    odtl, "hr",      false, true ));
    for(int k=0; k<odtl->gas->nSpecies(); k++)
        odtl->v.push_back(new lv_ygas(odtl, "y_"+odtl->gas->speciesName(k), true, true ));
    odtl->v.push_back(new lv_enth(  odtl, "enth",    true,  true ));

    int ii = 0;
    odtl->pos    = odtl->v.at(ii++);
    odtl->posf   = odtl->v.at(ii++);
    odtl->rho    = odtl->v.at(ii++);
    odtl->dvisc  = odtl->v.at(ii++);
    odtl->uvel   = odtl->v.at(ii++);
    odtl->vvel   = odtl->v.at(ii++);
    odtl->wvel   = odtl->v.at(ii++);
    odtl->temp   = odtl->v.at(ii++);
    odtl->mixf   = odtl->v.at(ii++);
    odtl->chi    = odtl->v.at(ii++);
    odtl->hr     = odtl->v.at(ii++);
    odtl->ysp = odtl->v.begin()+ii;       // access as odtl->ysp[k]->d[i], etc. where k is the species starting from 0.
    ii += odtl->gas->nSpecies();
    odtl->enth   = odtl->v.at(ii++);

    //------------------- set variables used for mesh adaption

    vector<lv*> phi;
    phi.push_back(odtl->uvel);
    phi.push_back(odtl->temp);
    odtl->mesher->init(odtl, phi);

    //------------------- set profiles

    string jetOrMixl  = odtl->io->initParams["jetOrMixl"].as<string>();
    double delta_mixf = odtl->io->initParams["delta_mixf"].as<double>();
    double fyc1       = odtl->io->initParams["fyc1"].as<double>();
    double delta_vel  = odtl->io->initParams["delta_vel"].as<double>();
    double vyc1       = odtl->io->initParams["vyc1"].as<double>();
    double vel_min    = odtl->io->initParams["vel_min"].as<double>();
    double vel_max    = odtl->io->initParams["vel_max"].as<double>();
    double vel_diff   = vel_max - vel_min;

    //--------------------

    if(jetOrMixl == "MIXL") {    // mixing layer

        fyc1 += 0.5*(odtl->posf->d.at(0)+odtl->posf->d.at(odtl->ngrd)); // fyc1=dist from center --> dist from left
        for(int i=0; i<odtl->ngrd; i++)
            odtl->mixf->d.at(i) = 0.5*(1.0+tanh(2.0/delta_mixf*(odtl->pos->d.at(i)-fyc1)));

        vyc1 += 0.5*(odtl->posf->d.at(0)+odtl->posf->d.at(odtl->ngrd));  // vyc1=dist from center --> dist from left
        for(int i=0; i<odtl->ngrd; i++)
            odtl->uvel->d.at(i) = vel_diff*0.5*(1.0+tanh(2.0/delta_vel*(odtl->pos->d.at(i)-vyc1))) + vel_min;
    }
    //--------------------

    else if(jetOrMixl == "JET") { // jet

        double fyc2       = odtl->io->initParams["fyc2"].as<double>();
        double vyc2       = odtl->io->initParams["vyc2"].as<double>();
        fyc1 += 0.5*(odtl->posf->d.at(0)+odtl->posf->d.at(odtl->ngrd)); // fyc1=dist from center --> dist from left
        fyc2 += 0.5*(odtl->posf->d.at(0)+odtl->posf->d.at(odtl->ngrd)); // fyc2=dist from center --> dist from left
        for(int i=0; i<odtl->ngrd; i++){
            odtl->mixf->d.at(i) = 0.5*(1.0+tanh(2.0/delta_mixf*(odtl->pos->d.at(i)-fyc1))) *
                               0.5*(1.0+tanh(2.0/delta_mixf*(fyc2-odtl->pos->d.at(i))));
        }

        vyc1 += 0.5*(odtl->posf->d.at(0)+odtl->posf->d.at(odtl->ngrd));  // vyc1=dist from center --> dist from left
        vyc2 += 0.5*(odtl->posf->d.at(0)+odtl->posf->d.at(odtl->ngrd));  // vyc2=dist from center --> dist from left
        for(int i=0; i<odtl->ngrd; i++){
            odtl->uvel->d.at(i) = vel_diff * 0.5 * (1.0 + tanh(2.0 / delta_vel * (odtl->pos->d.at(i) - vyc1))) *
                0.5 * (1.0 + tanh(2.0 / delta_vel * (vyc2 - odtl->pos->d.at(i)))) + vel_min;
        }
    }
    else {
        *odtl->io->ostrm << endl << "Error setting the odtline: option " << jetOrMixl
                         << " not recognized " << endl;
        exit(0);
    }

    //--------------------

    int nsp = odtl->gas->nSpecies();
    vector<double> ysp(nsp);               // dummy storage
    for(int i=0; i<odtl->ngrd; i++) {
        odtl->strm->getProdOfCompleteComb(odtl->mixf->d.at(i), ysp, odtl->enth->d.at(i), odtl->temp->d.at(i));
        for(int k=0; k<nsp; k++)
            odtl->ysp[k]->d.at(i) = ysp.at(k);
    }

    enforceMassFractions();

    odtl->rho->setVar();
    odtl->dvisc->setVar();

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the soltuion.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */

void odtcase_jetMixlRxn::setGasStateAtPt(const int &ipt) {

    int nsp = odtl->gas->nSpecies();
    vector<double> yi(nsp);
    for(int k=0; k<nsp; k++)
        yi.at(k) = odtl->ysp[k]->d.at(ipt);
    odtl->gas->setState_PY(odtl->odtp->pres, &yi.at(0));
    odtl->gas->setState_HP(odtl->enth->d.at(ipt), odtl->odtp->pres, 1.E-10);

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */
void odtcase_jetMixlRxn::setCaseSpecificVars() {

    enforceMassFractions();
    odtl->rho->setVar();
    odtl->dvisc->setVar();
    odtl->temp->setVar();

    odtl->enth->LagSrc = false;     // reset to false; in enth source the src is computed if false, then set to true on subsequent calls.
}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion. See cvodeDriver.cc
 *  These should not be transported. Those are already set.
 */
void odtcase_jetMixlRxn::setCaseSpecificVars_cvode(const int &ipt) {
    odtl->rho->setVar(ipt);
    odtl->temp->setVar(ipt);
}
