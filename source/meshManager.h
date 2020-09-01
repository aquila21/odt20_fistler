/**
 * @file meshManager.h
 * Header file for class meshManager
 */

#ifndef MESHMANAGER_H
#define MESHMANAGER_H

#include "lv.h"
#include <vector>

class odtline;

///////////////////////////////////////////////////////////////////////////////

/**
 * @class meshManager
 *
 * Mesh adaptor class.  This class will adapt odtlines or scalines etc.
 *  The specific line to adapt is set through a pointer.
 *  A single variable profile is used for adaption of a given line (e.g. velocity
 *  or mixture fraction).  The mesh adapter has several parts (see adaptODTgrid).
 *  This is meant to be called after eddy events or significant diffusion.  The
 *  adaption inserts or removes cells depending on cell size, gradients, curvature,
 *  and neighboring cell sizes.
 *
 *  @author David O. Lignell
 */

using namespace std;

class meshManager {

    ////////////////////// DATA MEMBERS /////////////////////

    public:

        odtline                  *odtl;     ///< pointer to odt line to adapt

        vector<lv*>              phi;       ///< vector of pointers to linevariable objects

        vector<double>           dx;        ///< vector of cell sizes
        vector<int>              mark;      ///< dummy small cell index array for sorting

        vector<double>           xf;        ///< vector of cell face positions
        vector<vector<double> >  yf;        ///< vector of cell values
        vector<double>           xnf;       ///< vector of new cell face positions
        vector<double>           X;         ///< vector of cell center positions

        int                      ngrd;      ///< local number of grid points
        int                      ngrdf;     ///< local number of grid faces

        int                      iLower;    ///< region of grid to adapt (the cell w/ left eddy edge)
        int                      iUpper;    ///< region of grid to adapt (the cell w/ right eddy edge)
        double                   posLower;  ///< physical region of eddy (lower bound)
        double                   posUpper;  ///< physical region of eddy (upper bound)

        vector<double>           lastDA;    ///< constant (unif) mesh to list time of last adapt

        ////////////////////// MEMBER FUNCTIONS  /////////////////////

        void adaptGrid(int iLower, int iUpper);       // Public interface

        inline bool operator()(const int &a, const int &b) const {
            return dx[mark[a]] < dx[mark[b]];
        }                                                // Functor used for sorting

        void adaptAfterSufficientDiffTime(const double &time,
                                          double       &tLastDA,
                                          int          &cLastDA,
                                          double       &dtCUmax);

        void adaptEddyRegionOfMesh(const double &time,
                                   double &tLastDA,
                                   int &cLastDA);

        void splitCell(const int isplt,
                       const int nsplt,
                       const vector<double> &cellFaces);

        void enforceDomainSize();
        //void removeFaceNearZero();

        //--------------- line operations

        void setGridDxc(const odtline *line, vector<double> &dxc, double C);
        void setGridDx(const odtline *line, vector<double>  &dx);

        void setGridFromDxc(const vector<double> &dxc2);



    private:

        void adaptGrid_details(const int iLower,
                               const int iUpper);        // called by adaptGrid
        void mergeSmallCells();
        void mergeSmallCellsMP();
        void impose2point5rule();
        void splitLargeCells();

        void fix2point5offender(const int mPos,
                                const int &iglobal);
        void setDxArray();

        int  findPos(const vector<double> &x,
                     const double         val,
                     const int            &istart);

        void interp1pt(const vector<double> &x,
                       const vector<double> &y,
                       const double         &xval,
                             double         &yval,
                             int            &istart);

        void interpVec(const vector<double> &x,
                       const vector<double> &y,
                       const vector<double> &xn,
                             vector<double> &yn);

        double calcDistance(const vector<double>          &x,
                            const vector<vector<double> > &y,
                                  vector<double>          &sDist);
        void set_iLower_iUpper();

        void updateDA(const double &time, double &tLastDA, int &cLastDA,
                      int iStart, int iEnd);

		//void splitCellWithZero();
        //--------------- line operations

        void merge2cells(const int imrg, const bool LconstVolume=false);

        ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public:

        void init(odtline *p_odtl, const vector<lv*> p_phi);
        meshManager(){};
        ~meshManager(){}

};

////////////////////////////////////////////////////////////////////////////////////

#endif

