/**
 * @file solver_hips.h
 * Header file for class solver_hips
 */

#ifndef SOLVER_HIPS_H
#define SOLVER_HIPS_H

#include <vector>
#include <utility>    // pair
#include "solver.h"

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing solver_hips object
 *
 *  @author David O. Lignell
 *
 */

class solver_hips : public solver {


    //////////////////// DATA MEMBERS //////////////////////

    private:

        vector<double> levelRates;     ///< list of eddy event rates at each level
        vector<pair<double,int> > eTL; ///< list of eddy times and levels
        int Nm1;                       ///< Nlevels - 1
        int Nm3;                       ///< Nlevels - 3



    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void calculateSolution();

    private:

        void setEddyEventTimes();
        void selectAndSwapTwoSubtrees(const int iLevel, int &Qstart, int &Rstart, int &nPswap);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        solver_hips(){}
        virtual void init(odtline *p_odtl);
        virtual ~solver_hips();

};


////////////////////////////////////////////////////////////////////////////////

#endif

