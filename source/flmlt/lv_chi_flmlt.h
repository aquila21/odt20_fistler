/**
 * @file lv_chi_flmlt.h
 * Header file for class lv_chi_flmlt
 */

#ifndef LV_CHI_FLMLT_H
#define LV_CHI_FLMLT_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_chi_flmlt of parent lv object.
 *
 *  @author David O. Lignell
 */

class lv_chi_flmlt : public lv {


    //////////////////// DATA MEMBERS //////////////////////

    public:

        double chiStoic;

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:
        virtual void setVar(const int ipt=-1);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_chi_flmlt(){}
        lv_chi_flmlt(odtline     *line,
               const string       s,
               const bool         Lt,
               const bool         Lo=true);

        virtual ~lv_chi_flmlt(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
