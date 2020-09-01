/**
 * @file lv_sca.h
 * Header file for class lv_sca
 */

#ifndef LV_SCA_H
#define LV_SCA_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_sca of parent lv object.
 *
 *  @author David O. Lignell
 *  @author Marten Klein
 */

class lv_sca : public lv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        //virtual void setVar(const int ipt=-1);
        virtual void getRhsSrc(const int ipt=-1);
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_sca(){}
        lv_sca(odtline     *line,
                const string s,
                const bool   Lt,
                const bool   Lo=true);

        virtual ~lv_sca(){}

};

////////////////////////////////////////////////////////////////////////////////

#endif
