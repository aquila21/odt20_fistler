/**
 * @file odtcase_hips_simpleRxn.h
 * Header file for class odtcase_hips_simpleRxn
 */

#ifndef ODTCASE_HIPS_SIMPLERXN_H
#define ODTCASE_HIPS_SIMPLERXN_H

#include "odtcase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child odtcase_hips_comb of parent odtcase object.
 *
 *  @author David O. Lignell
 */

class odtcase_hips_simpleRxn : public odtcase {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        odtcase_hips_simpleRxn(){}
        virtual void init(odtline *p_odtl);
        ~odtcase_hips_simpleRxn(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
