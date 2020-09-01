/**
 * @file odtcase_RT.cc
 * Header file for class odtcase_RT
 */

#include "odtcase_RT.h"
#include "odtline.h"
#include "lv.h"
#include "lv_pos.h"
#include "lv_posf.h"
#include "lv_rho_mf.h"
#include "lv_dvisc_const.h"
#include "lv_uvw.h"
#include "lv_mixf.h"

#include <cmath>
#include <string>

////////////////////////////////////////////////////////////////////////////////
/** odtcase_RT initialization function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void odtcase_RT::init(odtline *p_odtl) {

    odtl = p_odtl;

    odtl->v.push_back(new lv_pos(         odtl, "pos",     false, true  ));   // last are: L_transported, L_output
    odtl->v.push_back(new lv_posf(        odtl, "posf",    false, true  ));
    odtl->v.push_back(new lv_uvw(         odtl, "uvel",    true,  true  ));
    odtl->v.push_back(new lv_uvw(         odtl, "vvel",    true,  true  ));
    odtl->v.push_back(new lv_uvw(         odtl, "wvel",    true,  true  ));
    odtl->v.push_back(new lv_mixf(        odtl, "mixf",    true,  true  ));
    odtl->v.push_back(new lv_rho_mf(      odtl, "rho",     false, true  ));
    odtl->v.push_back(new lv_dvisc_const( odtl, "dvisc",   false, false ));

    int ii = 0;
    odtl->pos    = odtl->v.at(ii++);
    odtl->posf   = odtl->v.at(ii++);
    odtl->uvel   = odtl->v.at(ii++);
    odtl->vvel   = odtl->v.at(ii++);
    odtl->wvel   = odtl->v.at(ii++);
    odtl->mixf   = odtl->v.at(ii++);
    odtl->rho    = odtl->v.at(ii++);
    odtl->dvisc  = odtl->v.at(ii++);

    //------------------- set variables used for mesh adaption

    vector<lv*> phi;
    phi.push_back(odtl->uvel);
    phi.push_back(odtl->mixf);
    odtl->mesher->init(odtl, phi);

    //------------------- set profiles

    double Zcent = odtl->io->initParams["Zcent"].as<double>();
    double Ztran = odtl->io->initParams["Ztran"].as<double>();

    Zcent += 0.5*(odtl->posf->d.at(0)+odtl->posf->d.at(odtl->ngrd)); // Zcent=dist from center --> dist from left

    for(int i=0; i<odtl->ngrd; i++)
        odtl->mixf->d.at(i) = 0.5*(1.0+tanh(2.0/Ztran*(odtl->pos->d.at(i)-Zcent)));
        // odtl->mixf->d.at(i) = (odtl->pos->d.at(i) >= 0.5*(odtl->pos->d[0]+odtl->pos->d[odtl->ngrd-1])) ? 1.0 : 0.0;

    //--------------------

    odtl->rho->setVar();
    odtl->dvisc->setVar();

    //-------------------

    if(odtl->odtp->Lsolver!="EXPLICIT") {
            cout << "\nError Lsolver should be EXPLICIT for this case (no rxn sources)" << endl;
            exit(0);
    }

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */
void odtcase_RT::setCaseSpecificVars() {
    odtl->rho->setVar();
    odtl->dvisc->setVar();
}


