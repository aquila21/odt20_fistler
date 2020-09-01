/**
 * @file lv_rho_const.h
 * Header file for class lv_rho_const
 */

#ifndef LV_RHO_CONST_H
#define LV_RHO_CONST_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_rho_const of parent lv object.
 *
 *  @author David O. Lignell
 */

class lv_rho_const : public lv {

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

        lv_rho_const(){}
        lv_rho_const(odtline      *line,
                     const string s,
                     const bool   Lt,
                     const bool   Lo=true);

        virtual ~lv_rho_const(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
