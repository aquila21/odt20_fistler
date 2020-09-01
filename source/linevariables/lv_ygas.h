/**
 * @file lv_ygas.h
 * Header file for class lv_ygas
 */

#ifndef LV_YGAS_H
#define LV_YGAS_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_ygas of parent lv object.
 *
 *  @author David O. Lignell
 */

class lv_ygas : public lv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

        static int                     nspc;        ///< number of gas species

        int                            kMe;         ///< index of this spc in list: 0 to nspc-1; set from var_name

        double                         aP;
        double                         aP_x;

    public:

    //////////////////// MEMBER FUNCTIONS /////////////////


        virtual void getRhsSrc(const int ipt=-1);
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:

        void setFlux(const vector<double> &gf,
                     const vector<double> &dxc);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_ygas(){}
        lv_ygas(odtline     *line,
                const string s,
                const bool   Lt,
                const bool   Lo=true);

        virtual ~lv_ygas(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
