/**
 * @file lv_soot_flmlt_CQMOM.h
 * Header file for class lv_soot_flmlt_CQMOM
 *
 * https://en.wikipedia.org/wiki/Dominance_(C%2B%2B)
 */

#ifndef lv_soot_flmlt_CQMOM_H
#define lv_soot_flmlt_CQMOM_H

#include "lv_soot_CQMOM.h"
#include "lv_soot_flmlt.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_soot_CQMOM of parent lv object.
 *
 *  @author Victoria B. Lansinger
 *
 *  NOTE: getRhsMix will come from lv_soot_flmlt, and getRhsSrc will come from lv_soot_CQMOM
 */

class lv_soot_flmlt_CQMOM : public lv_soot_flmlt, public lv_soot_CQMOM {

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_soot_flmlt_CQMOM(odtline      *line,
                           const string  s,
                           const bool    Lt,
                           const bool    Lo=true) :
                 lv_soot(line, s, Lt, Lo),
                 lv_soot_flmlt(line, s, Lt, Lo),
                 lv_soot_CQMOM(line, s, Lt, Lo) {}

        virtual ~lv_soot_flmlt_CQMOM(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif


