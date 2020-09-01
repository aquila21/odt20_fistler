/**
 * @file micromixer_flmlt.h
 * Header file for class micromixer_flmlt
 */

#ifndef MICROMIXER_FLMLT_H
#define MICROMIXER_FLMLT_H

#include "micromixer.h"

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing micromixer_flmlt object
 *
 *  @author David O. Lignell
 */

class micromixer_flmlt : public micromixer {

    //////////////////// DATA MEMBERS //////////////////////

    private:

        double tNextAdapt;


    //////////////////// MEMBER FUNCTIONS /////////////////

    protected:

        virtual void setGf()              {return;}
        virtual void setGridDxcDx()       {return;}
        virtual void set_oldrho_or_rhov() {return;}
        virtual void setNominalStepSize();   ///< sets a nominal dt for the whole period

        virtual bool adaptGridsIfNeeded();

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        micromixer_flmlt() : micromixer() { tNextAdapt = -1.0; }
        virtual ~micromixer_flmlt(){ if(cvode) delete cvode; }

};


////////////////////////////////////////////////////////////////////////////////

#endif

