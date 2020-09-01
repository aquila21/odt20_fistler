/**
 * @file lv_soot_flmlt_MOMIC.h
 * Header file for class lv_soot_flmlt_MOMIC
 *
 * https://en.wikipedia.org/wiki/Dominance_(C%2B%2B)
 */

#ifndef lv_soot_flmlt_MOMIC_H
#define lv_soot_flmlt_MOMIC_H

#include "lv_soot_MOMIC.h"
#include "lv_soot_flmlt.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_soot_MOMIC of parent lv object.
 *
 *  @author Victoria B. Lansinger
 *
 *  NOTE: getRhsMix will come from lv_soot_flmlt, and getRhsSrc will come from lv_soot_MOMIC
 */

class lv_soot_flmlt_MOMIC : public lv_soot_flmlt, public lv_soot_MOMIC {

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_soot_flmlt_MOMIC(odtline      *line,
                           const string  s,
                           const bool    Lt,
                           const bool    Lo=true) :
                 lv_soot(line, s, Lt, Lo),
                 lv_soot_flmlt(line, s, Lt, Lo),
                 lv_soot_MOMIC(line, s, Lt, Lo) {}

        virtual ~lv_soot_flmlt_MOMIC(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif


