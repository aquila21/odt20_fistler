/**
 * @file lv_sdiff.h
 * Header file for class lv_sdiff
 */

#ifndef LV_SDIFF_H
#define LV_SDIFF_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_sdiff of parent lv object.
 *  Passive scalar diffusivity 
 *  
 *  @author David O. Lignell
 *  @author Marten Klein
 */

class lv_sdiff : public lv {

    public:

    //////////////////// DATA MEMBERS //////////////////////
        

    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setVar(const int ipt=-1);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_sdiff(){}
        lv_sdiff(odtline      *line,
                 const string s,
                 const bool   Lt,
                 const bool   Lo=true);

        virtual ~lv_sdiff(){}

};

////////////////////////////////////////////////////////////////////////////////

#endif
