/**
 * @file lv_rho_mf.cc
 * Header file for class lv_rho_mf
 */


#include "lv_rho_mf.h"
#include "odtline.h"

////////////////////////////////////////////////////////////////////////////////
/*! lv_rho_mf  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_rho_mf::lv_rho_mf(odtline    *line,
                           const      string s,
                           const bool Lt,
                           const bool Lo) {

    odtl          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(odtl->ngrd, odtl->odtp->rho0);

    rho0   = odtl->io->streamProps["rho0"].as<double>();
    rho1   = odtl->io->streamProps["rho1"].as<double>();

}

////////////////////////////////////////////////////////////////////////////////
/*! lv_rho_mf merger2cells function
 *
 * @param imrg    \input merge cells imrg and imrg+1
 * @param imrg \input merge cells imrg and imrg+1
 * @param m1   \input mass in cell imrg
 * @param m2   \input mass in cell imrg
 * @param LconstVolume \input (for posf, default is false)
 */

void lv_rho_mf::merge2cells(const int    imrg,
                               const double m1,
                               const double m2,
                               const bool   LconstVolume) {

    setVar(imrg);
    d.erase(d.begin() + imrg+1);
}

////////////////////////////////////////////////////////////////////////////////
/*! lv_rho_mf setVar function
 *  @param ipt \input optional point to compute at
 */

void lv_rho_mf::setVar(const int ipt) {

    d.resize(odtl->ngrd,odtl->odtp->rho0);
    if(ipt == -1)
        for(int i=0; i<odtl->ngrd; i++)
            d.at(i) = 1.0/( (1-odtl->mixf->d.at(i))/rho0 + odtl->mixf->d.at(i)/rho1 );
    else
        d.at(ipt) = 1.0/( (1-odtl->mixf->d.at(ipt))/rho0 + odtl->mixf->d.at(ipt)/rho1 );
}

