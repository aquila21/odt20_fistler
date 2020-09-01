/**
 * @file odtcase_channelScalar.cc
 * Header file for class odtcase_channelScalar
 */


#include "odtcase_channelScalar.h"
#include "odtline.h"
#include "lv.h"
#include "lv_pos.h"
#include "lv_posf.h"
#include "lv_rho_const.h"
#include "lv_dvisc_const.h"
#include "lv_uvw.h"
#include "lv_sca.h"

////////////////////////////////////////////////////////////////////////////////
/** Initialization
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void odtcase_channelScalar::init(odtline *p_odtl){

    odtl = p_odtl;

    odtl->v.push_back(new lv_pos(        odtl, "pos",   false, true));   // last are: L_transported, L_output
    odtl->v.push_back(new lv_posf(       odtl, "posf",  false, true));
    odtl->v.push_back(new lv_rho_const(  odtl, "rho",   false, false));
    odtl->v.push_back(new lv_dvisc_const(odtl, "dvisc", false, false));
    odtl->v.push_back(new lv_uvw(        odtl, "uvel",  true,  true));
    odtl->v.push_back(new lv_uvw(        odtl, "vvel",  true,  true));
    odtl->v.push_back(new lv_uvw(        odtl, "wvel",  true,  true));
    odtl->v.push_back(new lv_sca(        odtl, "sca",   true,  true));

    int k=0;
    odtl->pos   = odtl->v.at(k++);
    odtl->posf  = odtl->v.at(k++);
    odtl->rho   = odtl->v.at(k++);
    odtl->dvisc = odtl->v.at(k++);
    odtl->uvel  = odtl->v.at(k++);
    odtl->vvel  = odtl->v.at(k++);
    odtl->wvel  = odtl->v.at(k++);
    odtl->sca   = odtl->v.at(k++);

    //------------------- set variables used for mesh adaption

    vector<lv*> phi;
    phi.push_back(odtl->uvel);
    phi.push_back(odtl->sca);
    odtl->mesher->init(odtl, phi);

    //------------------- default vel. & sca. values (0.0) are fine, along with rho, dvisc.

    for(int i=0; i<odtl->ngrd; i++) {
        //---------- nonzero velocity profile (mind the BCs)
        //odtl->uvel->d[i] = 10*odtl->pos->d.at(i); //doldbg
        //odtl->uvel->d[i] = 10*0.016/4.0/0.002*(1.0-odtl->pos->d.at(i)*odtl->pos->d.at(i)); //doldbg
        //---------- non-zero scalar profile (mind the BCs)
        double xm = (odtl->pos->d.at(i) - odtl->posf->d.at(0))/odtl->Ldomain();
        odtl->sca->d.at(i) = odtl->odtp->sBClo * (1.0-xm) + odtl->odtp->sBChi * xm;
    }

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the soltuion.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */

void odtcase_channelScalar::setCaseSpecificVars() {

    odtl->rho->setVar();
    odtl->dvisc->setVar();

}
