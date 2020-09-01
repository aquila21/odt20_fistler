/**
 * @file lv_dvisc.h
 * Header file for class lv_dvisc
 */

#ifndef LV_DVISC_H
#define LV_DVISC_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_dvisc of parent lv object.
 *  Dynamic viscosity
 *
 *  @author David O. Lignell
 */

class lv_dvisc : public lv {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setVar(const int ipt=-1);

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_dvisc(){}
        lv_dvisc(odtline      *line,
                 const string s,
                 const bool   Lt,
                 const bool   Lo=true);

        virtual ~lv_dvisc(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
