/**
 * @file lv_soot_flmlt.h
 * Header file for class lv_soot_flmlt
 */

#ifndef lv_soot_FLMLT_H
#define lv_soot_FLMLT_H

#include "lv_soot.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_soot_flmlt of parent lv object.
 *
 *  @author Victoria B. Lansinger
 */

class lv_soot_flmlt : virtual public lv_soot {

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_soot_flmlt(odtline      *line,
                      const string  s,
                      const bool    Lt,
                      const bool    Lo=true) : lv_soot(line, s, Lt, Lo) {}

        virtual ~lv_soot_flmlt(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif


