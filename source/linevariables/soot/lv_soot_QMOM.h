/**
 * @file lv_soot_QMOM.h
 * Header file for class lv_soot_QMOM
 */

#ifndef lv_soot_QMOM_H
#define lv_soot_QMOM_H

#include "lv_soot.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_soot_QMOM of parent lv object.
 *
 *  @author Victoria B. Lansinger
 */

class lv_soot_QMOM : virtual public lv_soot {

    //////////////////// DATA MEMBERS //////////////////////

    private:

        vector<double>        wts;        ///< weights from inversion algorithm
        vector<double>        absc;       ///< abscissas from inversion algoritm

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsSrc(const int ipt=-1);

    private:

        int     factorial(int n);
        double  binomCoeff(int r, int k);
        double  Mk(double exp);
        void    getWtsAbs(vector<double> M, vector<double> &wts, vector<double> &abs, const int ipt);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_soot_QMOM(odtline      *line,
                     const string  s,
                     const bool    Lt,
                     const bool    Lo=true) : lv_soot(line, s, Lt, Lo) {}

        virtual ~lv_soot_QMOM(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif


