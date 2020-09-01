/**
 * @file lv_soot_MONO.h
 * Header file for class lv_soot_MONO
 */

#ifndef lv_soot_MONO_H
#define lv_soot_MONO_H

#include "lv_soot.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_soot_MONO of parent lv object.
 *
 *  @author Victoria B. Lansinger
 */

class lv_soot_MONO : virtual public lv_soot {

    //////////////////// DATA MEMBERS //////////////////////

    private:

        vector<double>        wts;        ///< weights of the particle size distribution
        vector<double>        absc;       ///< abscissas of the particle size distribution

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsSrc(const int ipt=-1);


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_soot_MONO(odtline      *line,
                     const string  s,
                     const bool    Lt,
                     const bool    Lo=true) : lv_soot(line, s, Lt, Lo) {}

        virtual ~lv_soot_MONO(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif


