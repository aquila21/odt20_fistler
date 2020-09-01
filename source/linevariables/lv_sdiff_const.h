/**
 * @file lv_sdiff_const.h
 * Header file for class lv_sdiff_const
 */

#ifndef LV_SDIFF_CONST_H
#define LV_SDIFF_CONST_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_sdiff_const of parent lv object.
 *  Passive scalar diffusivity 
 *  
 *  @author David O. Lignell
 *  @author Marten Klein
 */

class lv_sdiff_const : public lv {

    public:

    //////////////////// DATA MEMBERS //////////////////////
        

    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setVar(const int ipt=-1);

        virtual void merge2cells(const int    imrg,
                                 const double m2,
                                 const double m1,
                                 const bool   LconstVolume=false);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_sdiff_const(){}
        lv_sdiff_const(odtline      *line,
                       const string s,
                       const bool   Lt=false,
                       const bool   Lo=false);

        virtual ~lv_sdiff_const(){}

};

////////////////////////////////////////////////////////////////////////////////

#endif
