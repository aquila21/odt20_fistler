/**
 * @file odtcase_HIT.cc
 * Header file for class odtcase_HIT
 */


#include "odtcase_HIT.h"
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
/** odtcase_HIT initialization function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void odtcase_HIT::init(odtline *p_odtl) {

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
    for(int i=0; i<odtl->ngrd; i++){
        odtl->uvel->d.at(i) = 0.0;
        odtl->vvel->d.at(i) = 0.0;
        odtl->wvel->d.at(i) = 0.0;
    }
    
/*
    double u0 = odtl->io->initParams["turbInt"].as<double>();
    double pi = 3.141592654;  

    for(int j=1; j<7; j++){
        double r1 = odtl->rand->getRand();
        double r2 = odtl->rand->getRand();
        double r3 = odtl->rand->getRand();
        if(j < 3){
            for(int i=0; i<odtl->ngrd; i++){
                odtl->uvel->d.at(i) += 0.5*u0*sin(j*odtl->pos->d.at(i)+r1*2*pi);
                odtl->vvel->d.at(i) += 0.5*u0*sin(j*odtl->pos->d.at(i)+r2*2*pi);
                odtl->wvel->d.at(i) += 0.5*u0*sin(j*odtl->pos->d.at(i)+r3*2*pi);
            }
        }
        else{
            for(int i=0; i<odtl->ngrd; i++){
	        odtl->uvel->d.at(i) += u0/j*sin(j*odtl->pos->d.at(i)+r1*2*pi);
	        odtl->vvel->d.at(i) += u0/j*sin(j*odtl->pos->d.at(i)+r1*2*pi);
	        odtl->wvel->d.at(i) += u0/j*sin(j*odtl->pos->d.at(i)+r1*2*pi);
	    }
        }
    }
    double uAv = 0.0;
    double vAv = 0.0;
    double wAv = 0.0;

    for(int i=0; i<odtl->ngrd; i++){
        uAv += odtl->uvel->d.at(i);
        vAv += odtl->vvel->d.at(i);
        wAv += odtl->wvel->d.at(i);
    }
    uAv = uAv/odtl->ngrd;
    vAv = vAv/odtl->ngrd;
    wAv = wAv/odtl->ngrd;
    for(int i=0; i<odtl->ngrd; i++){
        odtl->uvel->d.at(i) -= uAv;
        odtl->vvel->d.at(i) -= vAv;
        odtl->wvel->d.at(i) -= wAv;
    }
*/
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
void odtcase_HIT::setCaseSpecificVars() {

    odtl->rho->setVar();
    odtl->dvisc->setVar();
}

