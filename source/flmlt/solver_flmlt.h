/**
 * @file solver_flmlt.h
 * Header file for class solver_flmlt
 */

#ifndef SOLVER_FLMLT_H
#define SOLVER_FLMLT_H

#include "solver.h"

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing solver_flmlt object
 *
 *  @author David O. Lignell
 */

class solver_flmlt : public solver {


    //////////////////// DATA MEMBERS //////////////////////

    public:


    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void calculateSolution();

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        solver_flmlt(){}
        virtual void init(odtline *p_odtl);
        virtual ~solver_flmlt(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif

