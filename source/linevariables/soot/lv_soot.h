/**
 * @file lv_soot.h
 * Header file for class lv_soot
 */

#ifndef LV_SOOT_H
#define LV_SOOT_H

#include "lv.h"
#include "soot_chem.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_soot of parent lv object.
 *
 *  @author Victoria B. Lansinger
 */

class lv_soot : public lv {

    //////////////////// DATA MEMBERS //////////////////////

    public:

        static soot_chem*   chem;           ///< pointer to soot_chem object
        static int          N;              ///< iterator for assigning kMe values

        int                 kMe;            ///< kMe
        int                 nsvar;          ///< number of soot variables

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsSrc(const int ipt=-1) = 0;
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:

        void setFlux(const vector<double> &gf,
                     const vector<double> &dxc);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_soot(odtline     *line,
                const string s,
                const bool   Lt,
                const bool   Lo=true);

        virtual ~lv_soot(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif


