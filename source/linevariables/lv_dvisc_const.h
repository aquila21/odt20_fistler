/**
 * @file lv_dvisc_const.h
 * Header file for class lv_dvisc_const
 */

#ifndef LV_DVISC_CONST_H
#define LV_DVISC_CONST_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_dvisc_const of parent lv object.
 *  Dynamic viscosity
 *
 *  @author David O. Lignell
 */

class lv_dvisc_const : public lv {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setVar(const int ipt=-1);

        virtual void   merge2cells(const int    imrg,
                                   const double m2,
                                   const double m1,
                                   const bool   LconstVolume=false);

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_dvisc_const(){}
        lv_dvisc_const(odtline      *line,
                       const string s,
                       const bool   Lt,
                       const bool   Lo=true);

        virtual ~lv_dvisc_const(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
