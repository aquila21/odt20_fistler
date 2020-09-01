/**
 * @file odtcase_coldPropaneJet.cc
 * Header file for class odtcase_coldPropaneJet
 */


#include "odtcase_coldPropaneJet.h"
#include "odtline.h"
#include "lv.h"
#include "lv_pos.h"
#include "lv_posf.h"
#include "lv_rho.h"
#include "lv_dvisc.h"
#include "lv_uvw.h"
#include "lv_ygas_noRxn.h"
#include "lv_mixf.h"
#include "lv_chi.h"
#include "lv_aDL.h"

#include "interp_linear.h"

#include <cmath>
#include <string>

////////////////////////////////////////////////////////////////////////////////
/** odtcase_coldPropaneJet initialization function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void odtcase_coldPropaneJet::init(odtline *p_odtl) {

    odtl = p_odtl;

    vector<double> gammas(4,0.0);
    gammas[0] = 1.0;
    odtl->strm->init(odtl,gammas);

    odtl->v.push_back(new lv_pos(   odtl, "pos",     false, true ));   // last are: L_transported, L_output
    odtl->v.push_back(new lv_posf(  odtl, "posf",    false, true ));
    odtl->v.push_back(new lv_rho(   odtl, "rho",     false, true ));
    odtl->v.push_back(new lv_dvisc( odtl, "dvisc",   false, true ));
    odtl->v.push_back(new lv_uvw(   odtl, "uvel",    true,  true ));
    odtl->v.push_back(new lv_uvw(   odtl, "vvel",    true,  true ));
    odtl->v.push_back(new lv_uvw(   odtl, "wvel",    true,  true ));
    odtl->v.push_back(new lv_mixf(  odtl, "mixf",    false, true ));
    odtl->v.push_back(new lv_chi(   odtl, "chi",     false, true ));
    odtl->v.push_back(new lv_aDL(   odtl, "aDL",     false, false ));
    for(int k=0; k<odtl->gas->nSpecies(); k++)
        odtl->v.push_back(new lv_ygas_noRxn(odtl, "y_"+odtl->gas->speciesName(k), true, true ));

    int ii = 0;
    odtl->pos    = odtl->v.at(ii++);
    odtl->posf   = odtl->v.at(ii++);
    odtl->rho    = odtl->v.at(ii++);
    odtl->dvisc  = odtl->v.at(ii++);
    odtl->uvel   = odtl->v.at(ii++);
    odtl->vvel   = odtl->v.at(ii++);
    odtl->wvel   = odtl->v.at(ii++);
    odtl->mixf   = odtl->v.at(ii++);
    odtl->chi    = odtl->v.at(ii++);
    odtl->aDL    = odtl->v.at(ii++);
    odtl->ysp = odtl->v.begin()+ii;       // access as odtl->ysp[k]->d[i], etc. where k is the species starting from 0.
    ii += odtl->gas->nSpecies();

    //------------------- set variables used for mesh adaption

    vector<lv*> phi;
    phi.push_back(odtl->uvel);
    odtl->mesher->init(odtl, phi);

    //------------------- set profiles

    double djeti = odtl->io->initParams["djeti"].as<double>();
    double djeto = odtl->io->initParams["djeto"].as<double>();

    vector<double> xprof, uprof;
    for(int i=0; i<odtl->io->initParams["vprof"].size(); i++){
        xprof.push_back(odtl->io->initParams["vprof"][i][0].as<double>());
        uprof.push_back(odtl->io->initParams["vprof"][i][1].as<double>());
    }

    Linear_interp Linterp(xprof, uprof);
    for(int i=0; i<odtl->ngrd; i++)
        odtl->uvel->d.at(i) = Linterp.interp(odtl->pos->d.at(i));

    for(int i=0; i<odtl->ngrd; i++)
        odtl->mixf->d.at(i) = (odtl->pos->d.at(i) >= -djeti/2.0 && odtl->posf->d.at(i) <= djeti/2.0) ? 1.0 : 0.0;

    //--------------------

    for(int i=0; i<odtl->ngrd; i++)
        for(int k=0; k<odtl->gas->nSpecies(); k++)
            odtl->ysp[k]->d.at(i) = odtl->mixf->d.at(i) * odtl->strm->y1[k] + (1.0-odtl->mixf->d.at(i))*odtl->strm->y0[k];

    enforceMassFractions();

    odtl->rho->setVar();
    odtl->dvisc->setVar();
    odtl->aDL->setVar();

    //-------------------

    if(odtl->odtp->Lsolver!="EXPLICIT") {
            cout << "\nError Lsolver needs to be EXPLICIT for this case (no rxn sources)" << endl;
            exit(0);
    }

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the soltuion.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */

void odtcase_coldPropaneJet::setGasStateAtPt(const int &ipt) {

    int nsp = odtl->gas->nSpecies();
    vector<double> yi(nsp);
    for(int k=0; k<nsp; k++)
        yi.at(k) = odtl->ysp[k]->d.at(ipt);
    odtl->gas->setState_TPY( odtl->strm->T0, odtl->odtp->pres, &yi[0] );

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */
void odtcase_coldPropaneJet::setCaseSpecificVars() {

    enforceMassFractions();
    odtl->rho->setVar();
    odtl->dvisc->setVar();
}

