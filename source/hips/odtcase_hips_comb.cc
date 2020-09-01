/**
 * @file odtcase_hips_comb.cc
 * Header file for class odtcase_hips_comb
 */

#include "odtcase_hips_comb.h"
#include "odtline.h"
#include "lv.h"
#include "lv_rho.h"
#include "lv_temp.h"
#include "lv_mixf_hips.h"
#include "lv_enth_hips.h"
#include "lv_ygas_hips.h"

////////////////////////////////////////////////////////////////////////////////
/** odtcase_hips_comb initialization function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void odtcase_hips_comb::init(odtline *p_odtl) {

    odtl = p_odtl;

    vector<double> gammas(4,0.0);
    gammas[0] = 2.0/odtl->gas->atomicWeight(odtl->gas->elementIndex("C"));
    gammas[1] = 0.5/odtl->gas->atomicWeight(odtl->gas->elementIndex("H"));
    gammas[2] = -1.0/odtl->gas->atomicWeight(odtl->gas->elementIndex("O"));
    gammas[3] = 0.0;
    odtl->strm->init(odtl, gammas);

    odtl->v.push_back(new lv_mixf_hips(    odtl, "mixf",    true,  true ));
    odtl->v.push_back(new lv_rho(          odtl, "rho",     false, true ));
    odtl->v.push_back(new lv_temp(         odtl, "temp",    false, true ));
    odtl->v.push_back(new lv_enth_hips(    odtl, "enth",    true,  true ));
    for(int k=0; k<odtl->gas->nSpecies(); k++)
        odtl->v.push_back(new lv_ygas_hips(odtl, "y_"+odtl->gas->speciesName(k), true, true ));

    int ii = 0;
    odtl->mixf   = odtl->v.at(ii++);
    odtl->rho    = odtl->v.at(ii++);
    odtl->temp   = odtl->v.at(ii++);
    odtl->enth   = odtl->v.at(ii++);
    odtl->ysp    = odtl->v.begin()+ii;          // access as odtl->ysp[k]->d[i], etc. where k is the species starting from 0.
    ii += odtl->gas->nSpecies();

    //-------------------- initialize profiles

    for(int i=0; i<odtl->ngrd; i++)
        odtl->mixf->d.at(i) = i<odtl->ngrd/2 ? 0.0 : 0.0637;

    int nsp = odtl->gas->nSpecies();
    vector<double> ysp(nsp);               // dummy storage
    for(int i=0; i<odtl->ngrd; i++) {
        odtl->strm->getProdOfCompleteComb(odtl->mixf->d.at(i), ysp, odtl->enth->d.at(i), odtl->temp->d.at(i));
        for(int k=0; k<nsp; k++)
            odtl->ysp[k]->d.at(i) = ysp.at(k);
    }

    enforceMassFractions();

    odtl->rho->setVar();

    //------------------- set minimial mesher

    vector<lv*> phi;
    odtl->mesher->init(odtl, phi);

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the soltuion.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */

void odtcase_hips_comb::setGasStateAtPt(const int &ipt) {

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
void odtcase_hips_comb::setCaseSpecificVars() {

    enforceMassFractions();
    odtl->rho->setVar();
    odtl->temp->setVar();
}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion. See cvodeDriver.cc
 *  These should not be transported. Those are already set.
 */
void odtcase_hips_comb::setCaseSpecificVars_cvode(const int &ipt) {

    odtl->rho->setVar(ipt);
    odtl->temp->setVar(ipt);
}
