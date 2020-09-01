/**
 * @file odtcase_hips_comb.h
 * Header file for class odtcase_hips_comb
 */

#ifndef ODTCASE_HIPS_COMB_H
#define ODTCASE_HIPS_COMB_H

#include "odtcase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child odtcase_hips_comb of parent odtcase object.
 *
 *  @author David O. Lignell
 */

class odtcase_hips_comb : public odtcase {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void setGasStateAtPt(const int &ipt);
        virtual void setCaseSpecificVars();
        virtual void setCaseSpecificVars_cvode(const int &ipt);


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        odtcase_hips_comb(){}
        virtual void init(odtline *p_odtl);
        ~odtcase_hips_comb(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
