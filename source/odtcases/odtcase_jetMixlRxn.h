/**
 * @file odtcase_jetMixlRxn.h
 * Header file for class odtcase_jetMixlRxn
 */

#ifndef ODTCASE_JETMIXLRXN_H
#define ODTCASE_JETMIXLRXN_H

#include "odtcase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child odtcase_jetMixlRxn of parent odtcase object.
 *
 *  @author David O. Lignell
 */

class odtcase_jetMixlRxn : public odtcase {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setGasStateAtPt(const int &ipt);
        virtual void setCaseSpecificVars();
        virtual void setCaseSpecificVars_cvode(const int &ipt);

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        odtcase_jetMixlRxn(){}
        virtual void init(odtline *p_odtl);
        ~odtcase_jetMixlRxn(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
