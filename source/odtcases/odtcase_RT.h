/**
 * @file odtcase_RT.h
 * Header file for class odtcase_RT
 */

#ifndef ODTCASE_RT_H
#define ODTCASE_RT_H

#include "odtcase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child odtcase_RT of parent odtcase object.
 *
 *  @author David O. Lignell
 */

class odtcase_RT : public odtcase {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setCaseSpecificVars();

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        odtcase_RT(){}
        virtual void init(odtline *p_odtl);
        ~odtcase_RT(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
