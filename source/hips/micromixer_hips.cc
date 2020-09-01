
#include "micromixer_hips.h"
#include "odtline.h"

///////////////////////////////////////////////////////////////////////////////
/** Set the nominal time step size.
 *  The actual step size will be rest based on a data dump or tend.
 */

void micromixer_hips::setNominalStepSize() {
    dtStepNominal = odtl->odtp->diffCFL * odtl->solv->tMix;
}

///////////////////////////////////////////////////////////////////////////////

