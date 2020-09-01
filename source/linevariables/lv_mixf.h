/**
 * @file lv_mixf.h
 * Header file for class lv_mixf
 */

#ifndef LV_MIXF_H
#define LV_MIXF_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_mixf of parent lv object.
 *
 *  @author David O. Lignell
 */

class lv_mixf : public lv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

        double Dmf;     ///<  default mixture fraction diffusivity (read from input file)

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void setVar(const int ipt=-1);
        virtual void getRhsSrc(const int ipt=-1);
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:

        void setFlux(const vector<double> &gf,
                     const vector<double> &dxc);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_mixf(){}
        lv_mixf(odtline     *line,
                const string s,
                const bool   Lt,
                const bool   Lo=true);

        virtual ~lv_mixf(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
