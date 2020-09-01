/**
 * @file lv_enth.h
 * Header file for class lv_enth
 */

#ifndef LV_ENTH_H
#define LV_ENTH_H

#include "lv.h"
#include "radiation.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_enth of parent lv object.
 *
 *  @author David O. Lignell
 */

class lv_enth : public lv {


    //////////////////// DATA MEMBERS //////////////////////

    private:
        int nspc;
        int iMe;
        radiation  *rad;

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsSrc(const int ipt=-1);
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:

        void setFlux(const vector<double> &gf,
                     const vector<double> &dxc);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_enth(){}
        lv_enth(odtline     *line,
                const string s,
                const bool   Lt,
                const bool   Lo=true);

        virtual ~lv_enth(){ if(rad) delete rad; }

};


////////////////////////////////////////////////////////////////////////////////

#endif
