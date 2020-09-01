/**
 * @file lv_uvw.h
 * Header file for class lv_uvw
 */

#ifndef LV_UVW_H
#define LV_UVW_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_uvw of parent lv object.
 *
 *  @author David O. Lignell
 */

class lv_uvw : public lv {

    //////////////////// DATA MEMBERS //////////////////////

    public:

    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void getRhsSrc(const int ipt=-1);
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_uvw(){}
        lv_uvw(odtline     *line,
               const string s,
               const bool   Lt,
               const bool   Lo=true);

        virtual ~lv_uvw(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
