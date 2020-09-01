/**
 * @file odtcase_jetFlame.cc
 * Header file for class odtcase_jetFlame
 */


#include "odtcase_jetFlame.h"
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
#include "lv_aDL.h"
#include "soot/lv_soot.h"
#include "soot/lv_soot_MONO.h"
//#include "soot/lv_soot_LOGN.h"
#include "soot/lv_soot_QMOM.h"
#include "soot/lv_soot_CQMOM.h"
#include "soot/lv_soot_MOMIC.h"

#include "interp_linear.h"

#include <cmath>
#include <string>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
/** odtcase_jetFlame initialization function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void odtcase_jetFlame::init(odtline *p_odtl) {

    odtl = p_odtl;

    LisFlameD = odtl->io->initParams["LisFlameD"] ? odtl->io->initParams["LisFlameD"].as<bool>() : false;

    vector<double> gammas(4,0.0);
    gammas[0] = 2.0/odtl->gas->atomicWeight(odtl->gas->elementIndex("C"));
    gammas[1] = 0.5/odtl->gas->atomicWeight(odtl->gas->elementIndex("H"));
    gammas[2] = -1.0/odtl->gas->atomicWeight(odtl->gas->elementIndex("O"));
    gammas[3] = 0.0;
    if(LisFlameD) gammas[2] = 0.0;
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
    //odtl->v.push_back(new lv_aDL(   odtl, "aDL",     false, false ));
    for(int k=0; k<odtl->gas->nSpecies(); k++)
        odtl->v.push_back(new lv_ygas(odtl, "y_"+odtl->gas->speciesName(k), true, true ));
    odtl->v.push_back(new lv_enth(  odtl, "enth",    true,  true ));

    // Add soot moments to variable list
    if (odtl->odtp->Lsoot) {

        string PSD_method = odtl->io->sootParams["PSD_method"].as<string>();
        stringstream ss;

        if (PSD_method == "MONO") {
            odtl->odtp->nsvar = 2;
            for(int k=0; k<2; k++) {
                ss.str(""); ss.clear(); ss << k;
                odtl->v.push_back(new lv_soot_MONO(odtl, "M"+ss.str(), true, true ));
            }
        }
        //else if (PSD_method == "LOGN") {
        //    odtl->odtp->nsvar = 3;
        //    for(int k=0; k<3; k++) {
        //        odtl->v.push_back(new lv_soot_LOGN(odtl, "M"+to_string(k), true, true ));
        //    }
        //}
        else if (PSD_method == "QMOM") {
            for(int k=0; k<odtl->odtp->nsvar; k++) {
                ss.str(""); ss.clear(); ss << k;
                odtl->v.push_back(new lv_soot_QMOM(odtl, "M"+ss.str(), true, true ));
            }
        }   // end QMOM
//        else if (PSD_method == "CQMOM") {
//            odtl->odtp->nsvar = 2*odtl->odtp->nsvar_v/2*odtl->odtp->nsvar_s/2 + odtl->odtp->nsvar_v/2;  // stored and accessed by column
//            for (int k=0; k<odtl->odtp->nsvar_v; k++) {
//                odtl->v.push_back(new lv_soot_CQMOM(odtl, "M"+to_string(k)+","+'0', true, true ));          // first s column: M00, M10, M20, etc.
//            }
//            for (int k=1; k<odtl->odtp->nsvar_s; k++) {
//                for (int j=0; j<odtl->odtp->nsvar_v/2; j++) {
//                    odtl->v.push_back(new lv_soot_CQMOM(odtl, "M"+to_string(j)+","+to_string(k), true, true ));       // other s columns: M01, M11, M02, M12, etc.
//                }
//            }
//        }   // end CQMOM
        else if (PSD_method == "MOMIC") {
            for(int k=0; k<odtl->odtp->nsvar; k++) {
                ss.str(""); ss.clear(); ss << k;
                odtl->v.push_back(new lv_soot_MOMIC(odtl, "M"+ss.str(), true, true ));
            }
        }   // end MOMIC


    }   // end if Lsoot

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
    //odtl->aDL    = odtl->v.at(ii++);
    odtl->ysp = odtl->v.begin()+ii;       // access as odtl->ysp[k]->d[i], etc. where k is the species starting from 0.
    ii += odtl->gas->nSpecies();
    odtl->enth   = odtl->v.at(ii++);
    if (odtl->odtp->Lsoot == true) {
        odtl->svar     = odtl->v.begin()+ii;    // access as odtl->svar[k]->d[i], etc. where k is the species starting from 0.
        ii += odtl->odtp->nsvar;
    }

    //------------------- set variables used for mesh adaption

    vector<lv*> phi;
    phi.push_back(odtl->uvel);
    phi.push_back(odtl->temp);
    odtl->mesher->init(odtl, phi);

    //------------------- set profiles

    double d_f    = odtl->io->initParams["d_f"].as<double>();
    double d_p    = odtl->io->initParams["d_p"].as<double>();
    double Z_p    = odtl->io->initParams["Z_p"].as<double>();
    double dTrans = odtl->io->initParams["dTrans"].as<double>();
    double U_f    = odtl->io->initParams["U_f"].as<double>();
    double U_p    = odtl->io->initParams["U_p"].as<double>();
    double U_a    = odtl->io->initParams["U_a"].as<double>();

    //--------------------

    double r_f = d_f/2.0;
    double r_p = d_p/2.0;
    double p1;                 // dummy for setting profiles
    double p2;                 // dummy for setting profiles

    for(int i=0; i<odtl->ngrd; i++){
        p1 = 0.5*(1.0+tanh(2.0/dTrans*(odtl->pos->d.at(i)+r_f))) *
             0.5*(1.0+tanh(2.0/dTrans*(r_f-odtl->pos->d.at(i))));
        p2 = 0.5*(1.0+tanh(2.0/dTrans*(odtl->pos->d.at(i)+r_p))) *
             0.5*(1.0+tanh(2.0/dTrans*(r_p-odtl->pos->d.at(i))));
        odtl->mixf->d.at(i) = p1*(1.0-Z_p) + p2*Z_p;
        if(U_f > 0)
            odtl->uvel->d.at(i) = p1*(U_f-U_p-U_a) + p2*(U_p-U_a) + U_a;
    }

    if(U_f < 0) {
        vector<double> xprof, uprof;
        for(int i=0; i<odtl->io->initParams["vprof"].size(); i++){
            xprof.push_back(odtl->io->initParams["vprof"][i][0].as<double>());
            uprof.push_back(odtl->io->initParams["vprof"][i][1].as<double>());
        }
        Linear_interp Linterp(xprof, uprof);
        for(int i=0; i<odtl->ngrd; i++)
            odtl->uvel->d.at(i) = Linterp.interp(odtl->pos->d.at(i));
    }

    //------------------- set initial soot profile for QMOM
    if (odtl->odtp->Lsoot) {
        if (odtl->odtp->PSD_method == "QMOM" || odtl->odtp->PSD_method == "MOMIC") {
            double M0 = 1.0E0;
            double sigL = 3.0;
            double mavg = 1.0E-21;
            for (int k=0; k<odtl->odtp->nsvar; k++) {
                for(int j=0; j<odtl->ngrd; j++) {
                    odtl->svar[k]->d[j] = M0 * pow(mavg, k) * exp(0.5 * pow(k,2) * pow(sigL,2));
                }
            }
        }
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

void odtcase_jetFlame::setGasStateAtPt(const int &ipt) {

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
void odtcase_jetFlame::setCaseSpecificVars() {

    enforceSootMom();
    enforceMassFractions();
    odtl->enth->setVar();
    odtl->rho->setVar();
    odtl->dvisc->setVar();
    odtl->temp->setVar();
    odtl->chi->setVar();

    odtl->enth->LagSrc = false;     // reset to false; in enth source the src is computed if false, then set to true on subsequent calls.
}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion. See cvodeDriver.cc
 *  These should not be transported. Those are already set.
 */
void odtcase_jetFlame::setCaseSpecificVars_cvode(const int &ipt) {
    odtl->rho->setVar(ipt);
    odtl->temp->setVar(ipt);
}
