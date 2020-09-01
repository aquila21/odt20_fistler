/**
 * @file odtcase_flmlt.cc
 * Header file for class odtcase_flmlt
 */

#include "odtcase_flmlt.h"
#include "odtline.h"
#include "lv.h"
#include "lv_pos.h"
#include "lv_posf.h"
#include "lv_rho.h"
#include "lv_dvisc.h"
#include "lv_enth_flmlt.h"
#include "lv_temp.h"
#include "lv_ygas_flmlt.h"
#include "lv_chi_flmlt.h"
#include "lv_soot_flmlt_MONO.h"
//#include "lv_soot_flmlt_LOGN.h"
#include "lv_soot_flmlt_QMOM.h"
#include "lv_soot_flmlt_CQMOM.h"
#include "lv_soot_flmlt_MOMIC.h"

#include "interp_linear.h"

#include <cmath>
#include <string>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
/** odtcase_flmlt initialization function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void odtcase_flmlt::init(odtline *p_odtl) {

    odtl = p_odtl;

    vector<double> gammas(4,0.0);
    gammas[0] = 2.0/odtl->gas->atomicWeight(odtl->gas->elementIndex("C"));
    gammas[1] = 0.5/odtl->gas->atomicWeight(odtl->gas->elementIndex("H"));
    gammas[2] = -1.0/odtl->gas->atomicWeight(odtl->gas->elementIndex("O"));
    gammas[3] = 0.0;
    odtl->strm->init(odtl, gammas);

    odtl->v.push_back(new lv_pos(           odtl, "pos",     false, true ));   // last are: L_transported, L_output
    odtl->v.push_back(new lv_posf(          odtl, "posf",    false, true ));
    odtl->v.push_back(new lv_rho(           odtl, "rho",     false, true ));
    odtl->v.push_back(new lv_dvisc(         odtl, "dvisc",   false, true ));
    odtl->v.push_back(new lv_temp(          odtl, "temp",    false, true ));
    odtl->v.push_back(new lv_chi_flmlt(     odtl, "chi",     false, true ));
    for(int k=0; k<odtl->gas->nSpecies(); k++)
        odtl->v.push_back(new lv_ygas_flmlt(odtl, "y_"+odtl->gas->speciesName(k), true, true ));
    odtl->v.push_back(new lv_enth_flmlt(    odtl, "enth",    false,  true ));

    // Add soot moments to variable list
    if (odtl->odtp->Lsoot) {

        string PSD_method = odtl->io->sootParams["PSD_method"].as<string>();
        stringstream ss;

        if (PSD_method == "MONO") {
            odtl->odtp->nsvar = 2;
            for(int k=0; k<2; k++) {
                ss.str(""); ss.clear(); ss << k;
                odtl->v.push_back(new lv_soot_flmlt_MONO(odtl, "M"+ss.str(), true, true ));
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
                odtl->v.push_back(new lv_soot_flmlt_QMOM(odtl, "M"+ss.str(), true, true ));
            }
        }   // end QMOM
        //else if (PSD_method == "CQMOM") {
        //    odtl->odtp->nsvar = 2*odtl->odtp->nsvar_v/2*odtl->odtp->nsvar_s/2 + odtl->odtp->nsvar_v/2;  // stored and accessed by column
        //    for (int k=0; k<odtl->odtp->nsvar_v; k++) {
        //        odtl->v.push_back(new lv_soot_flmlt_CQMOM(odtl, "M"+to_string(k)+","+'0', true, true ));          // first s column: M00, M10, M20, etc.
        //    }
        //    for (int k=1; k<odtl->odtp->nsvar_s; k++) {
        //        for (int j=0; j<odtl->odtp->nsvar_v/2; j++) {
        //            odtl->v.push_back(new lv_soot_flmlt_CQMOM(odtl, "M"+to_string(j)+","+to_string(k), true, true ));       // other s columns: M01, M11, M02, M12, etc.
        //        }
        //    }
        //}   // end CQMOM
        else if (PSD_method == "MOMIC") {
            for(int k=0; k<odtl->odtp->nsvar; k++) {
                ss.str(""); ss.clear(); ss << k;
                odtl->v.push_back(new lv_soot_flmlt_MOMIC(odtl, "M"+ss.str(), true, true ));
            }
        }   // end MOMIC

    }

    int ii = 0;
    odtl->pos    = odtl->v.at(ii++);
    odtl->posf   = odtl->v.at(ii++);
    odtl->rho    = odtl->v.at(ii++);
    odtl->dvisc  = odtl->v.at(ii++);
    odtl->temp   = odtl->v.at(ii++);
    odtl->chi    = odtl->v.at(ii++);
    odtl->ysp    = odtl->v.begin()+ii;          // access as odtl->ysp[k]->d[i], etc. where k is the species starting from 0.
    ii += odtl->gas->nSpecies();
    odtl->enth   = odtl->v.at(ii++);
    if (odtl->odtp->Lsoot == true) {
        odtl->svar     = odtl->v.begin()+ii;    // access as odtl->svar[k]->d[i], etc. where k is the species starting from 0.
        ii += odtl->odtp->nsvar;
    }

    //------------------- set variables used for mesh adaption

    vector<lv*> phi;
    phi.push_back(odtl->temp);
    odtl->mesher->init(odtl, phi);

    //------------------- set profiles

    if(odtl->odtp->domainLength != 1.0) {
        cout << endl << "ERROR: for flmlt, domainLength is mixf and must be 1" << endl;
        exit(0);
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

    double dx = odtl->odtp->domainLength / odtl->ngrd;
    odtl->posf->d.at(0) = 0.0;
    for(int i=1; i<odtl->ngrdf; i++)
        odtl->posf->d.at(i) = odtl->posf->d.at(i-1) + dx;
    odtl->posf->d.at(odtl->ngrd) = odtl->odtp->domainLength;
    odtl->pos->setVar();

    int nsp = odtl->gas->nSpecies();
    vector<double> ysp(nsp);               // dummy storage
    for(int i=0; i<odtl->ngrd; i++) {
        odtl->strm->getProdOfCompleteComb(odtl->pos->d.at(i), ysp, odtl->enth->d.at(i), odtl->temp->d.at(i));
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

void odtcase_flmlt::setGasStateAtPt(const int &ipt) {

    int nsp = odtl->gas->nSpecies();
    vector<double> yi(nsp);
    for(int k=0; k<nsp; k++) {
        yi.at(k) = odtl->ysp[k]->d.at(ipt);
    }

    odtl->gas->setState_PY(odtl->odtp->pres, &yi.at(0));
    odtl->gas->setState_HP(odtl->enth->d.at(ipt), odtl->odtp->pres, 1.E-10);

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */
void odtcase_flmlt::setCaseSpecificVars() {

    enforceSootMom();
    enforceMassFractions();
    odtl->enth->setVar();
    odtl->rho->setVar();
    odtl->dvisc->setVar();
    odtl->temp->setVar();
    odtl->chi->setVar();
}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion. See cvodeDriver.cc
 *  These should not be transported. Those are already set.
 */
void odtcase_flmlt::setCaseSpecificVars_cvode(const int &ipt) {

    odtl->rho->setVar(ipt);
    odtl->temp->setVar(ipt);
}
