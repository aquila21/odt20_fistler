/**
 * @file odtcase_coldJet.cc
 * Header file for class odtcase_coldJet
 */


#include "odtcase_coldJet.h"
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

void odtcase_coldJet::init(odtline *p_odtl) {

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

    double djeti      = odtl->io->initParams["djeti"].as<double>();
    double vel_min    = odtl->io->initParams["vel_min"].as<double>();
    double vel_max    = odtl->io->initParams["vel_max"].as<double>();
    double delta_vel  = odtl->io->initParams["delta_vel"].as<double>();
    double vel_diff   = vel_max - vel_min;

	double ti = 0.05;
    double vyc1 = -0.5*djeti;
    double vyc2 =  0.5*djeti;

	double r = odtl->rand->getRand();
    int sign = 1;
    if(r >= 0 && r < 0.5)
        sign = -1;
    double b = 2*vyc2/(3.14*5);

    for(int i=0; i<odtl->ngrd; i++){
        odtl->uvel->d.at(i) = vel_diff * 0.5 * (1.0 + tanh(2.0 / delta_vel * (odtl->pos->d.at(i) - vyc1))) *
            0.5 * (1.0 + tanh(2.0 / delta_vel * (vyc2 - odtl->pos->d.at(i)))) + vel_min;
		if(odtl->pos->d.at(i) > vyc1 && odtl->pos->d.at(i) < vyc2){
            odtl->uvel->d.at(i) +=sign*ti*vel_max*cos(odtl->pos->d.at(i)/b);
        }	
    }

    //--------- hack in a vprof profile

    if(djeti < 0) {
        vector<double> xprof, uprof;
        for(int i=0; i<odtl->io->initParams["vprof"].size(); i++){
			double r = 4*odtl->rand->getRand();
			if(r >= 2.0)
        		r -= 4;
            xprof.push_back(odtl->io->initParams["vprof"][i][0].as<double>()*-1*djeti);
            uprof.push_back((odtl->io->initParams["vprof"][i][1].as<double>()+odtl->io->initParams["vRMSprof"][i][1].as<double>()*r)*odtl->io->initParams["vel_max"].as<double>());
        }
        Linear_interp Linterp(xprof, uprof);
        for(int i=0; i<odtl->ngrd; i++)
            odtl->uvel->d.at(i) = Linterp.interp(odtl->pos->d.at(i))+vel_min;
    }

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
void odtcase_coldJet::setCaseSpecificVars() {

    odtl->rho->setVar();
    odtl->dvisc->setVar();
}

