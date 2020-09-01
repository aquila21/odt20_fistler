/**
 * @file lv_soot_LOGN.h
 * Header file for class lv_soot_LOGN
 */

#ifndef lv_soot_LOGN_H
#define lv_soot_LOGN_H

#include "lv_soot.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_soot_LOGN of parent lv object.
 *
 *  @author Victoria B. Lansinger
 */

class lv_soot_LOGN : virtual public lv_soot {

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsSrc(const int ipt=-1);

    private:

        double Mk(double k, const int i);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_soot_LOGN(odtline      *line,
                     const string  s,
                     const bool    Lt,
                     const bool    Lo=true) : lv_soot(line, s, Lt, Lo) {}

        virtual ~lv_soot_LOGN(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif


