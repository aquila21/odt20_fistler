/**
 * @file odtcase_hips.h
 * Header file for class odtcase_hips
 */

#ifndef ODTCASE_HIPS_H
#define ODTCASE_HIPS_H

#include "odtcase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child odtcase_hips of parent odtcase object.
 *
 *  @author David O. Lignell
 */

class odtcase_hips : public odtcase {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        odtcase_hips(){}
        virtual void init(odtline *p_odtl);
        ~odtcase_hips(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
