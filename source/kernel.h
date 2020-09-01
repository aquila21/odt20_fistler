/**
 * @file kernel.h
 * Header file for class kernel
 */

#ifndef KERNEL_H
#define KERNEL_H

#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing eddy object
 *
 *  @author Marco Fistler
 */

class kernel {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        odtline             *odtl;              ///< pointer to line object

        int                 nkernels;           ///< no. of implemented kernel events
        double              injectionTime;      ///< time when next kernel event will be sampled
        double              meanTime;           ///< kernel event mean occurence time
        double              kernelTime;         ///< time intervall till next kernel event
        double              kernelSize;         ///< size of eddy
        double              leftEdge;           ///< left edge location of kernel event
        double              rightEdge;          ///< right edge location of kernel event
        bool                LperiodicKernel;    ///< a wrap-around kernel
        int                 iStart;
     	int                 iEnd;
    	int                 ngrd;
        double              E;
        double              DeltaE;

        vector<double>      pos0;               ///< initial cell position
        vector<double>      dCoef;              ///< coefficient of K kernel
        vector<double>      K;                  ///< kernel K
        vector<double>      dxc;                ///< \delta(x^cCoord) is prop. to cell volume
        vector<double>      rhoU2;
        vector<double>      DeltaU2TM;

    //////////////////// MEMBER FUNCTIONS /////////////////

        void   kernelEnergyInjection(double time);
    private:

	void   sampleKernelTime();
	void   sampleKernelPosition();
	void   tripMap();
	void   fillKernel();
	void   computeKernelCoeff();
	void   applyKernels();
    void   integrateTKE(int ii);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        kernel(){}
        void init(odtline *p_odtl);
        ~kernel(){}

};

////////////////////////////////////////////////////////////////////////////////

#endif

