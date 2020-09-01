/**
 * @file odtcase_jetFlame.h
 * Header file for class odtcase_jetFlame
 */

#ifndef ODTCASE_JETFLAME_H
#define ODTCASE_JETFLAME_H

#include "odtcase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child odtcase_jetFlame of parent odtcase object.
 *
 *  @author David O. Lignell
 */

class odtcase_jetFlame : public odtcase {


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

        odtcase_jetFlame(){}
        virtual void init(odtline *p_odtl);
        ~odtcase_jetFlame(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
