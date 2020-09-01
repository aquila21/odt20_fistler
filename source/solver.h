/**
 * @file solver.h
 * Header file for class solver
 */

#ifndef SOLVER_H
#define SOLVER_H

#include "eddy.h"

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing solver object
 *
 *  @author David O. Lignell
 */

class solver {


    //////////////////// DATA MEMBERS //////////////////////

    public:

        odtline        *odtl;          ///< pointer to line object

        eddy           *ed3;           ///< pointer to eddy object for thirds
        odtline        *eddl3;         ///< pointer to eddy line object

        double         time;           ///< odt time (during sampling)
        double         t0;             ///< time of last eddy event; diffusion left off here.
        double         dtSmean;        ///< initial mean eddy sample time
        double         dtCUmax;        ///< max time before catch up diff/eddy

        bool           LeddyAccepted;  ///< flag for accepted eddy
        int            iEtrials;       ///< number of eddy trials

        double         PaSum;          ///< sum of Pa of eddies
        int            nPaSum;         ///< number going into PaSum
        int            neddies;        ///< number of eddies accepted
        double         PaSumC;         ///< sum of Pa of eddies
        int            nPaSumC;        ///< number going into PaSum

        //---------- for hips interface (inherited)

        double         tMix;           ///< parcel mixing timescale
        vector<int>    pLoc;           ///< parcel index array for fast implementation of swaps

		//---------- stats --------
		int            iStat;
		int            ngrdS;     ///< Number of grid points. Constant grid assumed here
        int            ngrdfS;    ///< Number of grid faces. Constant grid assumed here
        double         LdomainS;       ///< stat domain size (same as odt)
		double 		   dtStat;
		
		vector<double> posS;           ///< vector of cell center positions
        vector<double> posfS;          ///< vector of cell face  positions


        vector<double> uMeanS;         ///< mean u velocities
        vector<double> vMeanS;         ///< mean v velocities
        vector<double> wMeanS;         ///< mean w velocities

        vector<double> vTrans;        ///< a vector to tranfer odt grid to stat grid

	vector< vector< vector< vector<double> > > >  edstat;
        // data accumulation for calculation of budget terms     

        vector< vector< vector< vector<double> > > >  cstat;
        // second variable for data accumulation                 

        vector<double>                                ctime;
        // accumulation of time; used for calculating the average

        vector< vector<double> >                      oldVars;
        // variable for saving the old state to calculate the 
        // differenc caused by an eddy or a diffusion step

		vector< vector< vector<double> > >             cstats;
        // temporary variable for statisic written to cstat

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void calculateSolution();
		
		void   initStats();
        void   odtGrd2statGrd(vector<double> odtposf, vector<double> odtvec);
        void   initCstats();
        void   statsTime(const double tStep);
		void   cstats2statGrd();
		void   statsSetOld();
        void   statsChange(const int &jj);
		void   statsOutput();

    private:

        bool   sampleEddyAndImplementIfAccepted();
        void   computeDtSmean();
        void   computeDtCUmax();
        double sampleDt();
        void   diffusionCatchUpIfNeeded(bool Ldoit=false);
        void   raiseDtSmean();
        void   lowerDtSmean();
		bool   testLES_elapsedTime(const double time, const double tauEddy);
        bool   testLES_fracDomain( const double eSize);
        bool   testLES_integralLength(const double time, const double eSize);
        bool   testLES_thirds();

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        solver(){ ed3=0; eddl3=0; }
        virtual void init(odtline *p_odtl);
        virtual ~solver();

};


////////////////////////////////////////////////////////////////////////////////

#endif

