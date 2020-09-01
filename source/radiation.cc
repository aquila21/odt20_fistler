/**
 * @file radiation.cc
 * Source file for class radiation
 */

#include "radiation.h"
#include "odtline.h"
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
/** Initialization function
 *  @param p_odtl \input pointer to odtline object
 */

void radiation::init(odtline *p_odtl) {

    odtl = p_odtl;

    if(odtl->odtp->radType=="NONE")
        return;

    nRadSp = 4;                // CH4, CO2, H2O, CO

    //radCoefs.resize(nRadSp,vector<double>(6,0));  // kp (=) 1/atm*m
    radCoefs.resize(nRadSp+1,vector<double>(6,0));  // kp (=) 1/atm*m

    //radCoefs[0][0] =  1.017015E+1;    // ch4; kp=2.798 at 1150K
    //radCoefs[0][1] = -7.947312E-03;
    //radCoefs[0][2] =  4.342446E-7;
    //radCoefs[0][3] =  1.048611E-9;
    //radCoefs[0][4] = -2.287861E-13;
    //radCoefs[0][5] =  0.000000E+0;
    //radCoefs[1][0] =  3.24442E+1;     // co2; kp=29.197 at 925K
    //radCoefs[1][1] =  7.537513E-02;
    //radCoefs[1][2] = -1.535140E-04;
    //radCoefs[1][3] =  9.48794E-8;
    //radCoefs[1][4] = -2.509259E-11;
    //radCoefs[1][5] =  2.447995E-15;
    //radCoefs[2][0] =  6.86948E+1;     // h2o; kp=4.474  at 1119K
    //radCoefs[2][1] = -1.52349E-01;
    //radCoefs[2][2] =  1.417848E-04;
    //radCoefs[2][3] = -6.620996E-8;
    //radCoefs[2][4] =  1.52415E-11;
    //radCoefs[2][5] = -1.373456E-15;
    //radCoefs[3][0] =  1.56536E+0;    // co; kp=2.501 at 1007 K
    //radCoefs[3][1] =  1.483914E-02;
    //radCoefs[3][2] = -2.656035E-05;
    //radCoefs[3][3] =  1.68798E-8;
    //radCoefs[3][4] = -4.674473E-12;
    //radCoefs[3][5] =  4.767887E-16;

    radCoefs[0][0] =  6.6334;     // ch4
    radCoefs[0][1] = -0.0035686;
    radCoefs[0][2] =  1.6682E-08;
    radCoefs[0][3] =  2.5611E-10;
    radCoefs[0][4] = -2.6558E-14;
    radCoefs[0][5] =  0;
    radCoefs[1][0] =  18.741;     // co2
    radCoefs[1][1] = -121.310;
    radCoefs[1][2] =  273.500;
    radCoefs[1][3] = -194.050;
    radCoefs[1][4] =  56.310;
    radCoefs[1][5] = -5.8169;
    radCoefs[2][0] = -0.23093;     // h2o
    radCoefs[2][1] = -1.12390;
    radCoefs[2][2] =  9.41530;
    radCoefs[2][3] = -2.99880;
    radCoefs[2][4] =  0.51382;
    radCoefs[2][5] = -1.86840E-5;
    radCoefs[3][0] =  4.7869;       // co < 750 K
    radCoefs[3][1] = -0.06953;
    radCoefs[3][2] =  2.95775E-4;
    radCoefs[3][3] = -4.25732E-7;
    radCoefs[3][4] =  2.02894E-10;
    radCoefs[3][5] =  0.0;
    radCoefs[4][0] =  10.09;         // co > 750 K
    radCoefs[4][1] = -0.01183;
    radCoefs[4][2] = 4.7753E-6;
    radCoefs[4][3] = -5.87209E-10;
    radCoefs[4][4] = -2.5334E-14;
    radCoefs[4][5] =  0.0;

    sigmaSB = 5.670E-8;       // W/m2*K4

    sootFactor = 1863.0;     // 1220

   if(odtl->odtp->radType=="OPTHIN" && odtl->odtp->TBClo != odtl->odtp->TBChi)
       cout << "\n#********** WARNING in opthinRad TBClo != TBChi, using TBChi" << endl;

   //-------------- set iRadIndx: rad species not in the mechanism are -1

       bool fmissing = false;

       iRadIndx.resize(nRadSp);

       int isp;

       isp = odtl->gas->speciesIndex("CH4");
       isp = (isp > 0) ? isp : odtl->gas->speciesIndex("ch4");
       iRadIndx[0] = isp;
       if(isp < 0) fmissing = true;

       isp = odtl->gas->speciesIndex("CO2");
       isp = (isp > 0) ? isp : odtl->gas->speciesIndex("co2");
       iRadIndx[1] = isp;
       if(isp < 0) fmissing = true;

       isp = odtl->gas->speciesIndex("H2O");
       isp = (isp > 0) ? isp : odtl->gas->speciesIndex("h2o");
       iRadIndx[2] = isp;
       if(isp < 0) fmissing = true;

       isp = odtl->gas->speciesIndex("CO");
       isp = (isp > 0) ? isp : odtl->gas->speciesIndex("co");
       iRadIndx[3] = isp;
       if(isp < 0) fmissing = true;

       if(fmissing)
           cout << endl << "Warning one or more radiating species missing from mechanism" << endl;
}

///////////////////////////////////////////////////////////////////////////////
/** Function computes radiatative heat source (W/kg = W/m3 / rho).
 *  Use optically thin for Imode == 1, twoflux for Imod == 2.
 *  @param xMoleSp   \input species mole fractions
 *  @param &temp     \input temperature profile
 *  @param radSource \output gas radiation volumetric heat source (W/kg)
 *  @param fvSoot    \input optional soot volume fraction
 */
void radiation::getRadHeatSource(const vector<vector<double> > &xMoleSp,
                                 const vector<double>          &temp,
                                 vector<double>                &radSource,
                                 const vector<double>          &fvSoot) {

    if(odtl->odtp->radType == "NONE")
        return;
    else if(odtl->odtp->radType == "OPTHIN")
            opthinRadHeatSource(xMoleSp, temp, radSource, fvSoot);
    else if(odtl->odtp->radType == "TWOFLUX")
            twoFluxRadHeatSource(xMoleSp, temp, radSource, fvSoot);
    else {
        *odtl->io->ostrm << endl << "ERROR: radType not recognized" << endl;
        exit(0);
    }

}

///////////////////////////////////////////////////////////////////////////////
/** Function computes optically thin volumetric radiative heat source.
 *  @param xMoleSp   \input species mole fractions
 *  @param &temp     \input temperature profile
 *  @param radSource \output gas radiation volumetric heat source (W/kg)
 *  @param fvSoot    \input optional soot volume fraction
 */

void radiation::opthinRadHeatSource(const vector<vector<double> > &xMoleSp,
                                    const vector<double>          &temp,
                                    vector<double>                &radSource,
                                    const vector<double>          &fvSoot) {

   double Kabsorption;
   double TBChi4 = pow(odtl->odtp->TBChi, 4.0);

   for(int i=0; i<odtl->ngrd; i++) {

       if(fvSoot.size() == 0)
           Kabsorption = getGasAbsorptionCoefficient( xMoleSp[i], temp[i], odtl->odtp->pres);
       else
           Kabsorption = getGasAbsorptionCoefficient( xMoleSp[i], temp[i], odtl->odtp->pres, fvSoot[i] );

       radSource[i] = -4*sigmaSB*Kabsorption*(pow(temp[i],4.0) - TBChi4) / odtl->rho->d[i];
   }
}

///////////////////////////////////////////////////////////////////////////////
/** Function computes mean absorption coefficient of gas
 *  \cond
 *  Kp = sum x_k * P * K_k
 *  K_k = a5*T^5 + a4*T^4 + a3*T^3 + a2*T^2 + a1*T + a0
 *  \endcond
 *
 *  \f[
 *  K_p = \sum_k{x_k * P * K_k}
 *  \f]
 *  \f[
 *  K_k = a_5 T^5 + a_4 T^4 + a_3 T^3 + a_2 T^2 + a_1 T + a_0
 *  \f]
 *
 *  @param xMole \input xMole[species] - Species are all species in mechanism
 *  @param T \input Input temperature (K)
 *  @param pressure \input (Pa)
 *  @param fvSoot \input optional soot volume fraction
 *  @return Plank mean absorption coeficient (1/m)
 */

double radiation::getGasAbsorptionCoefficient(const vector<double> &xMole,
                                              const double         &T,
                                              const double         &pressure,
                                              const double         &fvSoot) {

    double P = pressure/101325.;         // Pa --> atm

    double Kabs;
    double Ktot = 0.0;

    //for(int k=0; k<nRadSp; k++) {
    //    if(iRadIndx[k] < 0)          // this radiation species k is not in the mechanism
    //        continue;
    //    Kabs = radCoefs[k][5];
    //    for(int j=4; j>=0; j--)
    //        Kabs = Kabs * T + radCoefs[k][j];

    //    Ktot += xMole[iRadIndx[k]]*P*Kabs;
    //}

    for(int k=0; k<nRadSp; k++) {
        if(iRadIndx[k] < 0)          // this radiation species k is not in the mechanism
            continue;
        if(k==0) {          // ch4
            Kabs = radCoefs[k][4];
            for(int j=3; j>=0; j--)
                Kabs = Kabs * T + radCoefs[k][j];
        }
        else if(k==1) {     // co2
            Kabs = radCoefs[k][5];
            for(int j=4; j>=0; j--)
                Kabs = Kabs * 1000/T + radCoefs[k][j];
        }
        else if(k==2) {     // h2o
            Kabs = radCoefs[k][5];
            for(int j=4; j>=0; j--)
                Kabs = Kabs * 1000/T + radCoefs[k][j];
        }
        else if(k==3){      // co
            int kk = (T<=750 ? k : k+1);
            Kabs = radCoefs[kk][4];
            for(int j=3; j>=0; j--)
                Kabs = Kabs * T + radCoefs[kk][j];
        }

        Ktot += xMole[iRadIndx[k]]*P*Kabs;
    }




    if(fvSoot > 0.0)
        Ktot += sootFactor * fvSoot * T;

    return Ktot;
}

///////////////////////////////////////////////////////////////////////////////
/** Function computes radiative source (W/m3) at grid points using the two-flux model
*
 *
 *
 *  \cond
 *  Solving (13.30a) and (13.30b) in Modest 1993 Radiative Heat Transfer P. 492:
 *  dIp/dx =  2*k*Ib - 2*k*Ip,    bc x=0: Ip = sig*Tbclo^4 / pi
 *  dIm/dx = -2*k*Kb + 2*k*Im,    bc x=0: Im = sig*Tbchi^4 / pi
 *
 *  with del dot q = rad source = k*(4*pi*Ib - G) (eqn 8.54) = 2*pik*(2*Ib-Ip-Im) using (eqn 13.32)
 *
 *  Here, written in terms of q+ and q- not I+ and I- (Ip, Im).
 *  With particles:
 *     dq+/ds =  2kg*sigma*Tg^4 + 2*sigma*sum_j(kp_j * Tp_j^4) - 2kg*q+ - 2*q+*sum_j(kp_j)
 *     dq-/ds = -2kg*sigma*Tg^4 - 2*sigma*sum_j(kp_j * Tp_j^4) + 2kg*q- + 2*q-*sum_j(kp_j)
 *  or: qout = qin + gas emission + particle emission - gas absorption - particle absorption
 *
 *  The solution grid is dx between cell centers:
 *  odt grid:      | * |        *        |                       *                     |
 *  solution grid: * *          *                                *                     *
 *  with dx between stars in the solution grid
 *
 *  Verified against Example 13.4.  Also computed del dot q directly using (13.33) (verif).
 *  This was done for constant dx, and constant k, and uniform T.
 *  The two ODEs are integrated at odtline cell center points using a finite difference grid
 *    using implicit euler.
 *  That is: q+_i = [q+_{i-1} + ds*(2Kg*sig*Tg^4 + 2*sig*sum_j(Kp_j*Tp_j^4)) ]_i / (1+ds*(2Kg 2sum_j(Kp_j)))_i
 *  The terms in this equation are computed and stored in arrays.
 *  First compute q+, q-, then Gas source (W/m3) is dq+/ds - dq-/ds.
 *  The particle source W/m3 / #/m3 --> W/particle.
 *
 *  Particles in a given cell are assumed to live at the grid center for purposes of computing radiative sources
 *  \endcond
 *
 *  @param xMoleSp   \input species mole fractions
 *  @param &temp     \input temperature profile
 *  @param radSource \output gas radiation volumetric heat source (W/kg)
 *  @param fvSoot    \input optional soot volume fraction
 */

void radiation::twoFluxRadHeatSource(const vector<vector<double> > &xMoleSp,
                                     const vector<double>          &temp,
                                     vector<double>                &radSource,
                                     const vector<double>          &fvSoot) {


   int npts = odtl->ngrd + 2;             // add in the two boundaries
   vector<double> qp(npts);
   vector<double> qm(npts);
   vector<double> Kabs(npts);
   vector<double> dx(npts-1);
   vector<double> gasEmmTerm(npts);       // 2*kg*sigma*Tg^4

   //------------- Get the dx array

   dx[0] = 0.5*(odtl->posf->d[1]-odtl->posf->d[0]);
   for(int i=1; i<npts-2; i++)
       dx[i] = 0.5*(odtl->posf->d[i+1]-odtl->posf->d[i-1]);
   dx[npts-2] = 0.5*(odtl->posf->d[odtl->ngrdf-1]-odtl->posf->d[odtl->ngrdf-2]);

   //------------- Get the gas absorption coefficient

   bool Lhave_soot = (fvSoot.size() > 0) ? true : false;

   for(int i=0; i<npts-2; i++) {
       if(Lhave_soot)
           Kabs[i+1] = getGasAbsorptionCoefficient( xMoleSp[i], temp[i], odtl->odtp->pres, fvSoot[i] );
       else
           Kabs[i+1] = getGasAbsorptionCoefficient( xMoleSp[i], temp[i], odtl->odtp->pres );
   }

   Kabs[0] = Kabs[1];
   Kabs[npts-1] = Kabs[npts-2];

   //------------- Get the gas emmision term

   for(int i=0; i<npts-2; i++)
       gasEmmTerm[i+1] = 2.0*Kabs[i+1] * sigmaSB*pow(temp[i],4.0);
   gasEmmTerm[0]      = gasEmmTerm[1];
   gasEmmTerm[npts-1] = gasEmmTerm[npts-2];

   //------------- Get radiative fluxes: qp, qm
   // Marching qp low to high, qm high to low

   double qpBClo = sigmaSB * pow(odtl->odtp->TBClo,4.0);
   double qmBChi = sigmaSB * pow(odtl->odtp->TBChi,4.0);

   qp[0] = qpBClo;
       for(int i=1; i<npts; i++)
       qp[i] = ( qp[i-1] + dx[i-1]*gasEmmTerm[i] ) / ( 1. + 2.*Kabs[i]*dx[i-1] );

   qm[npts-1] = qmBChi;
       for(int i=npts-2; i>=0; i--)
           qm[i] = ( qm[i+1] + dx[i]*gasEmmTerm[i] ) / ( 1. + 2.*Kabs[i]*dx[i] );

   //------------- Get gas radiative source: (div q) (=) W/kg (additive source to E-bal)

    for(int i=0, ip=1; i<npts-2; i++, ip++)
       radSource[i] = -2.*Kabs[ip] * (2.*sigmaSB*pow(temp[i],4.0) - qp[ip] - qm[ip]) /
                      odtl->rho->d[i];
}

////////////////////////////////DOXYGEN DOCUMENTATION//////////////////////////////////

/** \fn void radiation::twoFluxRadHeatSource(const vector<vector<double> > &xMoleSp,
                                     const vector<double> &temp,
                                     vector<double> &radSource,
                                     const vector<double> &fvSoot)
 *
 *  Solving (13.30a) and (13.30b) in Modest 1993 Radiative Heat Transfer P. 492:
 *  \f[
 *      \frac{d I_p}{dx} = 2 k I_b - 2 k I_p,    bc: x=0: Ip = \sigma \frac{{Tbc_{lo}}^4 }{ \pi}
 *  \f]
 *  \f[
 *      \frac{d I_m}{dx} = -2 k K_b + 2 k I_m,    bc: x=0: Im = \sigma \frac{{Tbc_{hi}}^4}{ \pi}
 *  \f]
 *
 *  with
 *  \f[\nabla \cdot q = \text{rad source} = k (4 \pi I_b - G) \f]
 *  \f[\text{(eqn 8.54)} = 2 \pi k (2 I_b-I_p-I_m) \text{ using (eqn 13.32)} \f]
 *
 *  Here, written in terms of \fun{q+} and \fun{q-} not \fun{I+} and \fun{I-} (\fun{Ip}, \fun{Im}).
 *  With particles:
 *     \f[
 *         \frac{dq+}{ds} =  2kg * \sigma * {Tg}^4 + 2*\sigma*\sum_j{kp_j * {{Tp}_j}^4} - 2kg*{q+} - 2*{q+}*\sum_j{kp_j}
 *     \f]
 *     \f[
 *         \frac{d{q-}}{ds} = -2kg*\sigma*{Tg}^4 - 2*\sigma*\sum_j{kp_j * {{Tp}_j}^4} + 2kg*{q-} + 2*{q-}*\sum_j{kp_j}
 *		 \f]
 *  or:
 *     \f[
 *          q_{out} = q_{in} + \text{ gas emission } + \text{ particle emission } - \text{ gas absorption } - \text{ particle absorption }
 *     \f]
 *
 *  The solution grid is dx between cell centers:\vc{
 *  odt grid:      | * |        *        |                       *                     |
 *  solution grid: * *          *                                *                     *
 *  }
 *  with dx between stars in the solution grid
 *
 *  Verified against Example 13.4.  Also computed \fun{\nabla \cdot q} directly using (13.33) (verif).
 *  This was done for constant \fun{dx}, constant \fun{k}, and uniform \fun{T}.
 *  The two ODEs are integrated at odtline cell center points using a finite difference grid
 *    using implicit euler.
 *  That is:
 *  \f[
 *      {q+}_i = \frac{[{q+}_{i-1} + ds*(2Kg*\sigma*{Tg}^4 + 2*\sigma*\sum_j{Kp_j*{{Tp}_j}^4}) ]_i }{ [1+ds*(2Kg 2\sum_j{Kp_j})]_i}
 *  \f]
 *  The terms in this equation are computed and stored in arrays.
 *  First compute \fun{q+}, \fun{q-}, then Gas source (\fun{W/m^3}) is \fun{dq+/ds - dq-/ds}.
 *  The particle source \fun{(W/m^3) / (\#/m^3) \to W/\text{particle}}.
 *
 *  Particles in a given cell are assumed to live at the grid center for purposes of computing radiative sources.
 */
