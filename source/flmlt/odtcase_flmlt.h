/**
 * @file odtcase_flmlt.h
 * Header file for class odtcase_flmlt
 */

#ifndef ODTCASE_FLMLT_H
#define ODTCASE_FLMLT_H

#include "odtcase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child odtcase_flmlt of parent odtcase object.
 *
 *  @author David O. Lignell
 */

class odtcase_flmlt : public odtcase {


    //////////////////// DATA MEMBERS //////////////////////

    private:
        bool LisFlameD;


    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void setGasStateAtPt(const int &ipt);
        virtual void setCaseSpecificVars();
        virtual void setCaseSpecificVars_cvode(const int &ipt);


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        odtcase_flmlt(){}
        virtual void init(odtline *p_odtl);
        ~odtcase_flmlt(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
