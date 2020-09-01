/**
 * @file odtcase_hips_simpleRxn.cc
 * Header file for class odtcase_hips_simpleRxn
 */

#include "odtcase_hips_simpleRxn.h"
#include "odtline.h"
#include "lv.h"
#include "lv_mixf_hips.h"
#include "lv_ygas_cold_hips.h"

////////////////////////////////////////////////////////////////////////////////
/** odtcase_hips_simpleRxn initialization function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void odtcase_hips_simpleRxn::init(odtline *p_odtl) {

    odtl = p_odtl;

    odtl->v.push_back(new lv_mixf_hips(     odtl, "mixf", true, true ));
    odtl->v.push_back(new lv_ygas_cold_hips(odtl, "y_A",  true, true ));
    odtl->v.push_back(new lv_ygas_cold_hips(odtl, "y_B",  true, true ));
    odtl->v.push_back(new lv_ygas_cold_hips(odtl, "y_R",  true, true ));
    odtl->v.push_back(new lv_ygas_cold_hips(odtl, "y_P",  true, true ));

    int ii = 0;
    odtl->mixf   = odtl->v.at(ii++);
    odtl->ysp    = odtl->v.begin()+ii;          // access as odtl->ysp[k]->d[i], etc. where k is the species starting from 0.
    ii += 4;

    //-------------------- initialize profiles

    for(int i=0; i<odtl->ngrd; i++){
        odtl->mixf->d.at(i) = i<odtl->ngrd/2  ? 0.0 : 1.0;
        odtl->ysp[0]->d[i]  = i<odtl->ngrd/2  ? 1.0 : 0.0;
        odtl->ysp[1]->d[i]  = i<odtl->ngrd/2  ? 0.0 : 1.0;
        odtl->ysp[2]->d[i]  = 0.0;
        odtl->ysp[3]->d[i]  = 0.0;
    }

    enforceMassFractions();

    //------------------- set minimial mesher

    vector<lv*> phi;
    odtl->mesher->init(odtl, phi);

}

