/**
 * @file eddy.h
 * Header file for class eddy
 */

#ifndef EDDY_H
#define EDDY_H

#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing eddy object
 *
 *  @author David O. Lignell
 */

class eddy {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        odtline             *odtl;              ///< pointer to line object
        odtline             *eddl;              ///< pointer to eddy line object

        double              eddySize;           ///< size of eddy
        double              leftEdge;           ///< left edge location of eddy
        double              rightEdge;          ///< right edge location of eddy
        double              invTauEddy;         ///< inverse eddy timescale
        double              Pa;                 ///< eddy acceptance probability
        bool                LperiodicEddy;      ///< a wrap-around eddy
        vector<double>      cCoef;              ///< coefficient of K kernel
        vector<double>      bCoef;              ///< coefficient of J kernel
        vector<double>      K;                  ///< eddy kernel K
        vector<double>      dxc;                ///< \delta(x^cCoord) is prop. to cell "volume"
        vector<double>      pos0;               ///< initial eddy cell locations, for kernel

        double              esdp1;              ///< eddy size distribution parameters.
        double              esdp2;
        double              esdp3;
        double              esdp4;

        vector<double> cca,ccb,ccc,ccd;                 ///< polynomial coefficient arrays for cylindricalAnomalyHack

		double              Ufavre;             ///< save for particle TWC  
        double              Etot;
        double              eViscPenalty;
        double              ePeEddy;

        int                 eddyType;           ///< Eddy type for shear flow (HST)
        ///////////////// homogeneous shear flow //////////////////
        double              p3;                 ///< probability of eddy type 3 between x and z (p1 = p2 = (1-p3=/2))
    //////////////////// MEMBER FUNCTIONS /////////////////

        void   sampleEddySize();
        void   sampleEddyPosition();
        void   tripMap(odtline *line,    const int iS, int iE, const double C, const bool LsplitAtEddy=false);
        bool   eddyTau(const double Z_value, const double C);
        void   computeEddyAcceptanceProb(const double dtSample);
        bool   applyVelocityKernels(odtline *line, const int iS, const int iE);

    private:

        void   fillKernel();
        double eddyFavreAvgVelocity(const vector<double> &dxc);
        void   set_kernel_coefficients();

        double cylindricalAnomalyHack(double rL);



    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        eddy(){}
        void init(odtline *p_odtl, odtline *p_eddl);
        ~eddy(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif

