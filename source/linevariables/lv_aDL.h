/**
 * @file lv_aDL.h
 * Header file for class lv_aDL
 */

#ifndef LV_ADL_H
#define LV_ADL_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_mixf of parent lv object.
 *
 *  @author David O. Lignell
 */

class lv_aDL : public lv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:
        virtual void setVar(const int ipt=-1);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_aDL(){}
        lv_aDL(odtline     *line,
               const string s,
               const bool   Lt,
               const bool   Lo=true);

        virtual ~lv_aDL(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
