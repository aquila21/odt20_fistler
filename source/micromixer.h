/**
 * @file micromixer.h
 * Header file for class micromixer
 */

#ifndef MICROMIXER_H
#define MICROMIXER_H

#include <vector>
#include <string>
#include "cvodeDriver.h"

class odtline;
class particle;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing micromixer object
 *
 *  @author David O. Lignell
 */

class micromixer {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        odtline        *odtl;          ///< pointer to line object
        cvodeDriver    *cvode;         ///< pointer to cvode driver object for implicit ODE integration (stiff)
		particle	   *part;          ///< pointer to particle object

        double         tstart;
        double         time;           ///< current time
        double         tend;
        double         dtStepNominal;  ///< nominal step size
        double         dt;             ///< actual step size (shortened based on output or tend)

        vector<double> dxc;            ///< abs(\Delta(x^c))
        vector<double> dx;             ///< abs(\Delta(x))
        vector<double> gf;             ///< grid factor for derivatives: (df/dx) = gf * (f - f)

        bool           LdoDump;        ///<

        vector<double> uDL_1;          ///< for DL instability: old velocity
        vector<double> uDL_2;          ///< for DL instability: new velocity
        vector<double> xDL_1;          ///< for DL instability: = "old" cell center positions
        vector<double> xDL_2;          ///< for DL instability: = "new" cell center positions
        vector<double> posDL_old;      ///< for DL instability: = "new" cell center positions

        vector<double> oldrho_or_rhov; ///< store the old density for continuity

    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void advanceOdt(const double p_tstart, const double p_tend);

        void check_balance(int io);

    protected:

        virtual void setGf();                ///< sets the gf array
        virtual void setGridDxcDx();         ///< sets the dxc array
        virtual void set_oldrho_or_rhov();   ///< record old rho (or rho*u) for continuity
        virtual bool adaptGridsIfNeeded();   ///< expansion or contraction --> adapt
        virtual void setNominalStepSize();   ///< sets a nominal dt for the whole period

        void setStepSize();                  ///< set a local dt for interruptions (dump or tend)
        void updateGrid();                   ///< enforce the continuity condition: (e.g., rho*dx = const).
        void do_DL(string doWhat);

        void advanceOdtSingleStep_Explicit();
        void advanceOdtSingleStep_SemiImplicit();
        void advanceOdtSingleStep_StrangSplit();





    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        micromixer();
        void init(odtline *p_odtl);
        virtual ~micromixer(){ if(cvode) delete cvode; }

};


////////////////////////////////////////////////////////////////////////////////

#endif

