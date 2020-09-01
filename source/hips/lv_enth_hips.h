/**
 * @file lv_enth_hips.h
 * Header file for class lv_enth_hips
 */

#ifndef LV_ENTH_HIPS_H
#define LV_ENTH_HIPS_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_enth_hips of parent lv object.
 *
 *  @author David O. Lignell
 */

class lv_enth_hips : public lv {


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

        lv_enth_hips(){}
        lv_enth_hips(odtline     *line,
                const string       s,
                const bool         Lt,
                const bool         Lo=true);

        virtual ~lv_enth_hips(){ }

};


////////////////////////////////////////////////////////////////////////////////

#endif
