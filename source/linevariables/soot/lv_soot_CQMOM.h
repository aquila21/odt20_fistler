/**
 * @file lv_soot_CQMOM.h
 * Header file for class lv_soot_CQMOM
 */

#ifndef lv_soot_CQMOM_H
#define lv_soot_CQMOM_H

#include "lv_soot.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_soot_CQMOM of parent lv object.
 *
 *  @author Victoria B. Lansinger
 */

class lv_soot_CQMOM : virtual public lv_soot {

    //////////////////// DATA MEMBERS //////////////////////

    private:

        vector<vector<double> >     wts;        ///< weights
        vector<vector<double> >     absc_v;     ///< abscissas, v direction
        vector<vector<double> >     absc_s;     ///< abscissas, s direction

        int                         nused_v;    ///<
        vector<int>                 nused_s;    ///<

        int                         nsvar_v;    ///<
        int                         nsvar_s;    ///<

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsSrc(const int ipt=-1);

    private:

        double Mkj(double v, double s);
//        void adaptive_CQMOM(int N_v, int N_s, const vector<vector<double> > &m,
//                    const double &eabs_v, const double &eabs_s,
//                    const vector<double> &rmin_v, const vector<double> &rmin_s,
//                    int &Nuse_v, vector<int> &Nuse_s, vector<vector<double> > &W,
//                    vector<vector<double> > &X_v, vector<vector<double> > &X_s );
//        void vandermonde_solver(const vector<double> &x, vector<double> &w, const vector<double> &q, const int &n);
//        void adaptive_wheeler(const vector<double> &m, int N, const vector<double> &rmin, const double &eabs,
//                      int &Nout, vector<double> &w, vector<double> &x );

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_soot_CQMOM(odtline      *line,
                     const string  s,
                     const bool    Lt,
                     const bool    Lo=true) : lv_soot(line, s, Lt, Lo) {}

        virtual ~lv_soot_CQMOM(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif


