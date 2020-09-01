
#include "solver_flmlt.h"
#include "odtline.h"

///////////////////////////////////////////////////////////////////////////////
/**
 */

void solver_flmlt::init(odtline *p_odtl) {
    odtl = p_odtl;
}

///////////////////////////////////////////////////////////////////////////////
/** The actual solver
 */

void solver_flmlt::calculateSolution() {

    odtl->io->writeDataFile("odt_init.dat", 0.0);

    odtl->mesher->adaptGrid(0, odtl->ngrd-1);

    odtl->io->writeDataFile("odt_init_adpt.dat", 0.0);

    odtl->mimx->advanceOdt(0.0, odtl->odtp->tEnd);

}
