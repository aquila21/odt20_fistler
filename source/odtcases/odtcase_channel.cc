/**
 * @file odtcase_channel.cc
 * Header file for class odtcase_channel
 */


#include "odtcase_channel.h"
#include "odtline.h"
#include "lv.h"
#include "lv_pos.h"
#include "lv_posf.h"
#include "lv_rho_const.h"
#include "lv_dvisc_const.h"
#include "lv_uvw.h"

////////////////////////////////////////////////////////////////////////////////
/** Initialization
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void odtcase_channel::init(odtline *p_odtl){

    odtl = p_odtl;

    odtl->v.push_back(new lv_pos(        odtl, "pos",   false, true));   // last are: L_transported, L_output
    odtl->v.push_back(new lv_posf(       odtl, "posf",  false, true));
    odtl->v.push_back(new lv_rho_const(  odtl, "rho",   false, false));
    odtl->v.push_back(new lv_dvisc_const(odtl, "dvisc", false, false));
    odtl->v.push_back(new lv_uvw(        odtl, "uvel",  true,  true));
    odtl->v.push_back(new lv_uvw(        odtl, "vvel",  true,  true));
    odtl->v.push_back(new lv_uvw(        odtl, "wvel",  true,  true));

    odtl->pos   = odtl->v.at(0);
    odtl->posf  = odtl->v.at(1);
    odtl->rho   = odtl->v.at(2);
    odtl->dvisc = odtl->v.at(3);
    odtl->uvel  = odtl->v.at(4);
    odtl->vvel  = odtl->v.at(5);
    odtl->wvel  = odtl->v.at(6);

    //------------------- set variables used for mesh adaption

    vector<lv*> phi;
    phi.push_back(odtl->uvel);
    odtl->mesher->init(odtl, phi);

    //------------------- default velocity values (0.0) are fine, along with rho, dvisc.

    //for(int i=0; i<odtl->uvel->d.size(); i++)
    //  odtl->uvel->d[i] = 10*odtl->pos->d.at(i);
      //odtl->uvel->d[i] = 10*0.016/4.0/0.002*(1.0-odtl->pos->d.at(i)*odtl->pos->d.at(i));

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the soltuion.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */

void odtcase_channel::setCaseSpecificVars() {

    odtl->rho->setVar();
    odtl->dvisc->setVar();

}
