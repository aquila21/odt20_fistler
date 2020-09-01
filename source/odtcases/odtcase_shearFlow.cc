/**
 * @file odtcase_shearFlow.cc
 * Header file for class odtcase_shearFlow
 */


#include "odtcase_shearFlow.h"
#include "odtline.h"
#include "lv.h"
#include "lv_pos.h"
#include "lv_posf.h"
#include "lv_rho_const.h"
#include "lv_dvisc_const.h"
#include "lv_uvw.h"
#include "lv_ygas_noRxn.h"
#include "lv_mixf.h"
#include "lv_chi.h"
#include "lv_aDL.h"

#include "interp_linear.h"

#include <cmath>
#include <string>

////////////////////////////////////////////////////////////////////////////////
/** odtcase_coldJet initialization function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void odtcase_shearFlow::init(odtline *p_odtl) {

    odtl = p_odtl;

    odtl->v.push_back(new lv_pos(         odtl, "pos",     false, true ));   // last are: L_transported, L_output
    odtl->v.push_back(new lv_posf(        odtl, "posf",    false, true ));
    odtl->v.push_back(new lv_rho_const(   odtl, "rho",     false, true ));
    odtl->v.push_back(new lv_dvisc_const( odtl, "dvisc",   false, true ));
    odtl->v.push_back(new lv_uvw(         odtl, "uvel",    true,  true ));
    odtl->v.push_back(new lv_uvw(         odtl, "vvel",    true,  true ));
    odtl->v.push_back(new lv_uvw(         odtl, "wvel",    true,  true ));

    int ii = 0;
    odtl->pos    = odtl->v.at(ii++);
    odtl->posf   = odtl->v.at(ii++);
    odtl->rho    = odtl->v.at(ii++);
    odtl->dvisc  = odtl->v.at(ii++);
    odtl->uvel   = odtl->v.at(ii++);
    odtl->vvel   = odtl->v.at(ii++);
    odtl->wvel   = odtl->v.at(ii++);

    //------------------- set variables used for mesh adaption

    vector<lv*> phi;
    phi.push_back(odtl->uvel);
    odtl->mesher->init(odtl, phi);

    //------------------- set profiles

    double domainLength  = odtl->io->params["domainLength"].as<double>();
	double shearRate 	 = odtl->io->initParams["Srate"].as<double>();
	double pi            = 3.141592654;

    for(int i=0; i<odtl->ngrd; i++){
        odtl->uvel->d.at(i) = shearRate*(odtl->pos->d.at(i)+0.5*odtl->odtp->domainLength)+odtl->odtp->uBClo;
    }
		
	odtl->odtp->uBChi = shearRate*domainLength+odtl->odtp->uBClo;

    //--------------------

    odtl->rho->setVar();
    odtl->dvisc->setVar();

    //-------------------

    if(odtl->odtp->Lsolver!="EXPLICIT") {
            cout << "\nError Lsolver needs to be EXPLICIT for this case (no rxn sources)" << endl;
            exit(0);
    }

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */
void odtcase_shearFlow::setCaseSpecificVars() {

    odtl->rho->setVar();
    odtl->dvisc->setVar();
}

