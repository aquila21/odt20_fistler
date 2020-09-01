/**
 * @file lv_soot_MOMIC.h
 * Header file for class lv_soot_MOMIC
 */

#ifndef lv_soot_MOMIC_H
#define lv_soot_MOMIC_H

#include "lv_soot.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_soot_MOMIC of parent lv object.
 *
 *  @author Victoria B. Lansinger
 */

class lv_soot_MOMIC : virtual public lv_soot {

    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsSrc(const int ipt=-1);

    private:

        double  lagrangeInterp(double x_i, vector<double> x, vector<double> y);
        double  MOMIC(double p, vector<double> M);
        double  f_grid(int x, int y, vector<double> M);
        int     factorial(int n);
        double  binomCoeff(int r, int k);
        double  beta(int p, int q, int ipt);
        double  getCoag(double T, double P, double mu, vector<double> M, int r);
        void    downselectIfNeeded(int ipt, vector<double> &M, int &N);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_soot_MOMIC(odtline      *line,
                     const string  s,
                     const bool    Lt,
                     const bool    Lo=true) : lv_soot(line, s, Lt, Lo) {}

        virtual ~lv_soot_MOMIC(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif


