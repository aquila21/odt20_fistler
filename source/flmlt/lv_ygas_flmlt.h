/**
 * @file lv_ygas_flmlt.h
 * Header file for class lv_ygas_flmlt
 */

#ifndef LV_YGAS_FLMLT_H
#define LV_YGAS_FLMLT_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_ygas_flmlt of parent lv object.
 *
 *  @author David O. Lignell
 */

class lv_ygas_flmlt : public lv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

        static int                     nspc;        ///< number of gas species

        int                            kMe;         ///< index of this spc in list: 0 to nspc-1; set from var_name

    public:

    //////////////////// MEMBER FUNCTIONS /////////////////


        virtual void getRhsSrc(const int ipt=-1);
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_ygas_flmlt(){}
        lv_ygas_flmlt(odtline     *line,
                const string       s,
                const bool         Lt,
                const bool         Lo=true);

        virtual ~lv_ygas_flmlt(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
