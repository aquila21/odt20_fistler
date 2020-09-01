/**
 * @file micromixer_hips.h
 * Header file for class micromixer_hips
 */

#ifndef MICROMIXER_HIPS_H
#define MICROMIXER_HIPS_H

#include "micromixer.h"

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing micromixer_hips object
 *
 *  @author David O. Lignell
 */

class micromixer_hips : public micromixer {

    //////////////////// DATA MEMBERS //////////////////////

    private:


    //////////////////// MEMBER FUNCTIONS /////////////////

    protected:

        virtual void setGf()              {return;}
        virtual void setGridDxcDx()       {return;}
        virtual void set_oldrho_or_rhov() {return;}
        virtual void setNominalStepSize();     ///< sets a nominal dt for the whole period

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        micromixer_hips() : micromixer() { }
        virtual ~micromixer_hips(){ if(cvode) delete cvode; }

};


////////////////////////////////////////////////////////////////////////////////

#endif

