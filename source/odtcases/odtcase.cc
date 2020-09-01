/**
 * @file odtcase.cc
 * Header file for class odtcase
 */

#include "odtcase.h"
#include "odtline.h"

////////////////////////////////////////////////////////////////////////////////
/** Make sure mass fractions are normalized and bounded between 0 and 1
 */
void odtcase::enforceMassFractions() {

    double sum;
    for(int i=0; i<odtl->ngrd; i++) {
        sum = 0.0;
        for(int k=0; k<odtl->gas->nSpecies(); k++) {
            if(odtl->ysp[k]->d.at(i) < 0.0) odtl->ysp[k]->d.at(i) = 0.0;
            if(odtl->ysp[k]->d.at(i) > 1.0) odtl->ysp[k]->d.at(i) = 1.0;
            sum += odtl->ysp[k]->d.at(i);
        }
        for(int k=0; k<odtl->gas->nSpecies(); k++)
            odtl->ysp[k]->d.at(i) /= sum;
    }
}

////////////////////////////////////////////////////////////////////////////////
/** Make sure mass fractions are normalized and bounded between 0 and 1
 */
void odtcase::enforceSootMom() {

    if (!odtl->odtp->Lsoot)
        return;

    double sum;
    for(int i=0; i<odtl->ngrd; i++) {
        for(int k=0; k<odtl->odtp->nsvar; k++)
            if(odtl->svar[k]->d[i] < 0.0)
                odtl->svar[k]->d[i] = 0.0;
    }
}
