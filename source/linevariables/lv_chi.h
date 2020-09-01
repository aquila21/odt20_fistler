/**
 * @file lv_chi.h
 * Header file for class lv_chi
 */

#ifndef LV_CHI_H
#define LV_CHI_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_chi of parent lv object.
 *
 *  @author David O. Lignell
 */

class lv_chi : public lv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:
        virtual void setVar(const int ipt=-1);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_chi(){}
        lv_chi(odtline     *line,
               const string s,
               const bool   Lt,
               const bool   Lo=true);

        virtual ~lv_chi(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
