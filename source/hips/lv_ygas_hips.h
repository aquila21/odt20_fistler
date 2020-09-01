/**
 * @file lv_ygas_hips.h
 * Header file for class lv_ygas_hips
 */

#ifndef LV_YGAS_HIPS_H
#define LV_YGAS_HIPS_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_ygas_hips of parent lv object.
 *
 *  @author David O. Lignell
 */

class lv_ygas_hips : public lv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

        static int   nspc;        ///< number of gas species
        int          kMe;         ///< index of this spc in list: 0 to nspc-1; set from var_name

    public:

    //////////////////// MEMBER FUNCTIONS /////////////////


        virtual void getRhsSrc(const int ipt=-1);
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_ygas_hips(){}
        lv_ygas_hips(odtline     *line,
                const string       s,
                const bool         Lt,
                const bool         Lo=true);

        virtual ~lv_ygas_hips(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
