
#include "micromixer_flmlt.h"
#include "odtline.h"

///////////////////////////////////////////////////////////////////////////////
/** Set time step size.
 *  This is based on a diffusive (or other) timescale.
 *  This is a uniform step size for the given integration period.
 *  The actual step size will be rest based on a data dump or tend.
 */

void micromixer_flmlt::setNominalStepSize() {

    odtl->mesher->setGridDx(odtl, dx);

    double coef = 0.0;
    double dmb;
    for (int i=0; i < odtl->ngrd; i++) {
        dmb = 0.5*odtl->chi->d.at(i)/dx[i]/dx[i];
        if (dmb > coef)
            coef = dmb;
    }
    dtStepNominal = odtl->odtp->diffCFL * 0.5 / coef;
}


///////////////////////////////////////////////////////////////////////////////

/** Adapt during diffusion for spatial cases for which grid contraction
 *  results in small grid cells
 */

bool micromixer_flmlt::adaptGridsIfNeeded() {

    if(tNextAdapt == -1.0)   // just to initialize it
        tNextAdapt = 1.0/odtl->chi->chi0 * 10;
    if(time > tNextAdapt) {
        odtl->mesher->adaptGrid(0, odtl->ngrd-1);
        tNextAdapt += 1.0/odtl->chi->chi0 * 10;
        return true;
    }
    return false;
}

