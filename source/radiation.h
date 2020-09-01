/**
 * @file radiation.h
 * Header file for class radiation
 */

#ifndef RADIATION_H
#define RADIATION_H

#include <vector>
#include <string>

class odtline;

using namespace std;

///////////////////////////////////////////////////////////////////////////////

/** Class implementing radiation models: optically thin or two flux.
 *  Assumes gray gases
 *
 *  @author David O. Lignell
 */

class radiation {

  ////////////////////// DATA MEMBERS /////////////////////

    public:

        odtline                      *odtl;       ///< pointer to odtline

        double                       nRadSp;      ///< number of radiating species
        vector<vector<double> >      radCoefs;    ///< Radiation Coefficients [spc][coef]
        vector<int>                  iRadIndx;    ///< radiation species indicies: ch4 co2 h2o co: negative if not present

        double                       sigmaSB;     ///< Stefan Boltzman const

        double                       sootFactor;  ///< Ksoot = 1863 * fvsoot * T


    ////////////////////// MEMBER FUNCTIONS  /////////////////////

        void getRadHeatSource(const vector<vector<double> > &xMoleSp,
                              const vector<double>          &temp,
                              vector<double>                &radSource,
                              const vector<double>          &fvSoot = vector<double>(0,0.0));
    private:

        void opthinRadHeatSource(const vector<vector<double> > &xMoleSp,
                                 const vector<double>          &temp,
                                 vector<double>                &radSource,
                                 const vector<double>          &fvSoot=vector<double>(0,0.0));

        void twoFluxRadHeatSource(const vector<vector<double> > &xMoleSp,
                                  const vector<double>          &temp,
                                  vector<double>                &radSource_G,
                                  const vector<double>          &fvSoot=vector<double>(0,0.0));

        double getGasAbsorptionCoefficient(const vector<double> &xMole,
                                           const double         &T,
                                           const double         &pressure,
                                           const double         &fvSoot=0);


    ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public:

        radiation(){};   // constructor
        void init(odtline *p_odtl);
        ~radiation(){}

};

#endif
