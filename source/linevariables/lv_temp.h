/**
 * @file lv_temp.h
 * Header file for class lv_temp
 */

#ifndef LV_TEMP_H
#define LV_TEMP_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_temp of parent lv object.
 *
 *  @author David O. Lignell
 */

class lv_temp : public lv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:
        virtual void setVar(const int ipt=-1);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_temp(){}
        lv_temp(odtline     *line,
                const string s,
                const bool   Lt,
                const bool   Lo=true);

        virtual ~lv_temp(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
