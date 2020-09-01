/**
 t @file odtcase_coldPropaneJet.h
 * Header file for class odtcase_coldPropaneJet
 */

#ifndef ODTCASE_COLDPROPANEJET_H
#define ODTCASE_COLDPROPANEJET_H

#include "odtcase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child odtcase_coldPropaneJet of parent odtcase object.
 *
 *  @author David O. Lignell
 */

class odtcase_coldPropaneJet : public odtcase {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setGasStateAtPt(const int &ipt);
        virtual void setCaseSpecificVars();
        virtual void setCaseSpecificVars_cvode(const int &ipt) { }

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        odtcase_coldPropaneJet(){}
        virtual void init(odtline *p_odtl);
        ~odtcase_coldPropaneJet(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
