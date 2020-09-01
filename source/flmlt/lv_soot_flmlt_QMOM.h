/**
 * @file lv_soot_flmlt_QMOM.h
 * Header file for class lv_soot_flmlt_QMOM
 *
 * https://en.wikipedia.org/wiki/Dominance_(C%2B%2B)
 */

#ifndef lv_soot_flmlt_QMOM_H
#define lv_soot_flmlt_QMOM_H

#include "lv_soot_QMOM.h"
#include "lv_soot_flmlt.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_soot_QMOM of parent lv object.
 *
 *  @author Victoria B. Lansinger
 *
 *  NOTE: getRhsMix will come from lv_soot_flmlt, and getRhsSrc will come from lv_soot_QMOM
 */

class lv_soot_flmlt_QMOM : public lv_soot_flmlt, public lv_soot_QMOM {

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_soot_flmlt_QMOM(odtline      *line,
                           const string  s,
                           const bool    Lt,
                           const bool    Lo=true) :
                 lv_soot(line, s, Lt, Lo),
                 lv_soot_flmlt(line, s, Lt, Lo),
                 lv_soot_QMOM(line, s, Lt, Lo) {}

        virtual ~lv_soot_flmlt_QMOM(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif


