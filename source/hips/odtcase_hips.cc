/**
 * @file odtcase_hips.cc
 * Header file for class odtcase_hips
 */

#include "odtcase_hips.h"
#include "odtline.h"
#include "lv.h"
#include "lv_mixf_hips.h"

////////////////////////////////////////////////////////////////////////////////
/** odtcase_hips initialization function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void odtcase_hips::init(odtline *p_odtl) {

    odtl = p_odtl;

    odtl->v.push_back(new lv_mixf_hips(    odtl, "mixf",    true,  true ));

    int ii = 0;
    odtl->mixf   = odtl->v.at(ii++);

    //-------------------- initialize profiles

    for(int i=0; i<odtl->ngrd; i++)
        odtl->mixf->d.at(i) = i<odtl->ngrd/2 ? 0.0 : 1.0;

    //------------------- set minimial mesher

    vector<lv*> phi;
    odtl->mesher->init(odtl, phi);

}

