/**
 * @file lv_mixf_hips.h
 * Header file for class lv_mixf_hips
 */

#ifndef LV_MIXF_HIPS_H
#define LV_MIXF_HIPS_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_mixf_hips of parent lv object.
 *
 *  @author David O. Lignell
 */

class lv_mixf_hips : public lv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsSrc(const int ipt=-1);
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_mixf_hips(){}
        lv_mixf_hips(odtline     *line,
                     const string s,
                     const bool   Lt,
                     const bool   Lo=true);

        virtual ~lv_mixf_hips(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
