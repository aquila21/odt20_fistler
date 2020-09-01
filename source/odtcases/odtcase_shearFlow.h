/**
 * @file odtcase_shearFlow.h
 * Header file for class odtcase_shearFlow
 */

#ifndef ODTCASE_SHEARFLOW_H
#define ODTCASE_SHEARFLOW_H

#include "odtcase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child odtcase_shearFlow of parent odtcase object.
 *
 *  @author Marco Fistler
 */

class odtcase_shearFlow : public odtcase {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

    void setCaseSpecificVars();

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        odtcase_shearFlow(){}
        virtual void init(odtline *p_odtl);
        ~odtcase_shearFlow(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
