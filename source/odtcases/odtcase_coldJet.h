/**
 * @file odtcase_coldJet.h
 * Header file for class odtcase_coldJet
 */

#ifndef ODTCASE_COLDJET_H
#define ODTCASE_COLDJET_H

#include "odtcase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child odtcase_coldJet of parent odtcase object.
 *
 *  @author David O. Lignell
 */

class odtcase_coldJet : public odtcase {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

    void setCaseSpecificVars();

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        odtcase_coldJet(){}
        virtual void init(odtline *p_odtl);
        ~odtcase_coldJet(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
