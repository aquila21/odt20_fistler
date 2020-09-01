/**
 * @file lv_hr.h
 * Header file for class lv_hr
 */

#ifndef LV_HR_H
#define LV_HR_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_hr of parent lv object.
 *
 *  @author David O. Lignell
 */

class lv_hr : public lv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:
        virtual void setVar(const int ipt=-1);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_hr(){}
        lv_hr(odtline     *line,
              const string s,
              const bool   Lt,
              const bool   Lo=true);

        virtual ~lv_hr(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
