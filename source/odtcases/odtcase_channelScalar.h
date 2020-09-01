/**
 * @file odtcase_channelScalar.h
 * Header file for class odtcase_channelScalar
 */

#ifndef ODTCASE_CHANNEL_SCALAR_H
#define ODTCASE_CHANNEL_SCALAR_H

#include "odtcase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child odtcase_channelScalar of parent odtcase object.
 *
 *  @author David O. Lignell
 *  @author Marten Klein
 */

class odtcase_channelScalar : public odtcase {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setCaseSpecificVars();

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        odtcase_channelScalar(){}
        virtual void init(odtline *p_odtl);
        ~odtcase_channelScalar(){}

};

////////////////////////////////////////////////////////////////////////////////

#endif
