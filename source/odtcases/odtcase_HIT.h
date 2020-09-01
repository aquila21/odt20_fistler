/**
 * @file odtcase_HIT.h
 * Header file for class odtcase_HIT
 */

#ifndef ODTCASE_HIT_H
#define ODTCASE_HIT_H

#include "odtcase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child odtcase_HIT of parent odtcase object.
 *
 *  @author Marco Fistler
 */

class odtcase_HIT : public odtcase {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

    void setCaseSpecificVars();

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        odtcase_HIT(){}
        virtual void init(odtline *p_odtl);
        ~odtcase_HIT(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
