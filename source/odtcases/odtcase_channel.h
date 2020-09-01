/**
 * @file odtcase_channel.h
 * Header file for class odtcase_channel
 */

#ifndef ODTCASE_CHANNEL_H
#define ODTCASE_CHANNEL_H

#include "odtcase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child odtcase_channel of parent odtcase object.
 *
 *  @author David O. Lignell
 */

class odtcase_channel : public odtcase {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setCaseSpecificVars();

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        odtcase_channel(){}
        virtual void init(odtline *p_odtl);
        ~odtcase_channel(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
