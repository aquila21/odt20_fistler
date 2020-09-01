/**
 * @file lv_rho_mf.h
 * Header file for class lv_rho_mf
 */

#ifndef LV_RHO_MF_H
#define LV_RHO_MF_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_rho_mf of parent lv object.
 *  Density is based on simple mixing between two streams with densities rho0, rho1
 *  (Assumes ideal gas mixing with const MW, cp.)
 *
 *  @author David O. Lignell
 */

class lv_rho_mf : public lv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

        double rho0;         ///< read from input file (streams section)
        double rho1;         ///< read from input file (streams section)

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void setVar(const int ipt=-1);

        virtual void merge2cells(const int    imrg,
                                 const double m2,
                                 const double m1,
                                 const bool   LconstVolume=false);

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_rho_mf(){}
        lv_rho_mf(odtline      *line,
                  const string s,
                  const bool   Lt,
                  const bool   Lo=true);

        virtual ~lv_rho_mf(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
