/**
 * @file streams.cc
 * Header file for class streams
 */

#include "streams.h"
#include "odtline.h"
#include "yaml-cpp/yaml.h"
#include <cassert>

using namespace std;

///////////////////////////////////////////////////////////////////////////////
/** Initialize
 *
 * @param p_odtl \input pointer to odtline object.
 */
void streams::init(odtline *p_odtl, const vector<double> &gammas) {

    odtl = p_odtl;

    nspc = odtl->gas->nSpecies();

    T0   = odtl->io->streamProps["T0"].as<double>();
    T1   = odtl->io->streamProps["T1"].as<double>();

    y0   = vector<double>(nspc, 0.0);
    y1   = vector<double>(nspc, 0.0);

    bool Lmole = odtl->io->streamProps["moleOrMass"].as<string>() == "MOLE" ? true : false;

    gCHON = gammas;

    //--------------- stream properties: T,h,y.

    double sum = 0.0;              // for normalizing
    vector<string> y0Labels;       // species names in list
    vector<double> y0short;        // value corresponding to y0Labels
    YAML::Node yy = odtl->io->streamProps["comp0"];
    for(YAML::const_iterator it = yy.begin(); it!=yy.end(); it++) {
        y0Labels.push_back(it->first.as<string>());
        y0short.push_back(it->second.as<double>());
        sum += it->second.as<double>();
    }
    for(int i=0; i<y0Labels.size(); i++)
        y0[odtl->gas->speciesIndex(y0Labels[i])] = y0short[i]/sum;

    sum = 0.0;
    vector<string> y1Labels;       // species names in list
    vector<double> y1short;        // value corresponding to y0Labels
    yy = odtl->io->streamProps["comp1"];
    for(YAML::const_iterator it = yy.begin(); it!=yy.end(); it++) {
        y1Labels.push_back(it->first.as<string>());
        y1short.push_back(it->second.as<double>());
        sum += it->second.as<double>();
    }
    for(int i=0; i<y1Labels.size(); i++)
        y1[odtl->gas->speciesIndex(y1Labels[i])] = y1short[i] / sum;

    if(Lmole) {
        odtl->gas->setMoleFractions( &y0[0] );
        odtl->gas->getMassFractions( &y0[0] );
        odtl->gas->setMoleFractions( &y1[0] );
        odtl->gas->getMassFractions( &y1[0] );
    }

    odtl->gas->setState_TPY( T0, odtl->odtp->pres, &y0[0] );
    h0 = odtl->gas->enthalpy_mass();
    odtl->gas->setState_TPY( T1, odtl->odtp->pres, &y1[0] );
    h1 = odtl->gas->enthalpy_mass();

    //---------------- stoich mixf

    setStoicMixf();

    getMixtureFraction(&y0[0], true);    // set beta0 and beta1

}
///////////////////////////////////////////////////////////////////////////////
/** Computes the temperature, enthalpy, and composition of mixing among streams.
 *
 *  @param mixf \input mixture fraction, defines elemental composition.
 *  @param ymix \output mass fractions of products of complete combustion.
 *  @param hmix \output enthalpy of products of complete combustion.
 *  @param Tmix \output temperature of products of complete combustion.
 */
void streams::getMixingState(const double mixf, vector<double> &ymix,
                             double &hmix, double &Tmix){

    hmix = h1*mixf + h0*(1.0-mixf);
    for(int k=0; k<nspc; k++)
        ymix[k] = y1[k]*mixf + y0[k]*(1.0-mixf);

    odtl->gas->setMassFractions( &ymix[0] );
    odtl->gas->setState_HP(hmix, odtl->odtp->pres, 1.E-10);    // get temperature as Tadiabatic

    Tmix = odtl->gas->temperature();
}

///////////////////////////////////////////////////////////////////////////////
/** Computes the temperature, enthalpy, and composition of complete combustion at the given mixf.
 *  For nonpremixed flames (don't do anything funny, like have oxygen in the fuel stream)
 *
 *  @param mixf \input mixture fraction, defines elemental composition.
 *  @param ypcc \output mass fractions of products of complete combustion.
 *  @param hpcc \output enthalpy of products of complete combustion.
 *  @param Tpcc \output temperature of products of complete combustion.
 */
void streams::getProdOfCompleteComb(const double mixf, vector<double> &ypcc,
                                    double &hpcc, double &Tpcc){

    //---------- Compute the mixing mass fractions and enthalpy

    ypcc.resize(nspc);
    double d1 = 1.0-mixf;
    for(int k=0; k<nspc; k++) {
        ypcc[k] = y1[k]*mixf + y0[k]*d1;
    }
    hpcc = h1*mixf + h0*d1;

    //--------- Set gas and element indicicies

    int iC = odtl->gas->elementIndex("C");
    int iH = odtl->gas->elementIndex("H");
    //int iO = odtl->gas->elementIndex("O"); // !!!!!  currently unusedvariable
    int iN = odtl->gas->elementIndex("N");
    int iCO2 = odtl->gas->speciesIndex("CO2");
    int iH2O = odtl->gas->speciesIndex("H2O");
    int iN2  = odtl->gas->speciesIndex("N2");
    int iO2  = odtl->gas->speciesIndex("O2");

    //---------- Set ypcc as the mixing mole fractions: Take a basis of one mole:
    // now we are working in moles
    // elemM are moles of each element
    // when stoic:
    // CxHyNz   + (x+y/4)O2  ==>  (z/2)N2 + (x)CO2 + (y/2)H2O
    // otherwise:
    // CxHyNz   + (beta)O2   ==>  (z/2)N2 + (x)CO2 + (y/2)H2O
    // Note this allows fuels with nitrogen and oxygen

    odtl->gas->setMassFractions( &ypcc[0] );
    odtl->gas->getMoleFractions( &ypcc[0] );

    double nOnotFromO2  = 0.0;
    double nHnotFromH2O = 0.0;
    double nCnotFromCO2 = 0.0;
    vector<double> elemM = getElementMoles( &ypcc[0], nOnotFromO2,
            nHnotFromH2O, nCnotFromCO2 );

    double x    = elemM[iC];
    double y    = elemM[iH];
    double z    = elemM[iN];
    double beta = elemM[odtl->gas->elementIndex("O")] * 0.5;        // moles of O as O2

    //--------------------------------------------------------------------------

    if(mixf < mixfStoic) {                        // lean: burn all fuel, leftover O2

        for(int k=0; k<nspc; k++)
            ypcc[k] = 0.0;

        if(iCO2 > 0)                              // CO2 is not in the H2 mechs
            ypcc[iCO2] = x;
        ypcc[iH2O] = y*0.5;
        ypcc[iN2]  = z*0.5;
        ypcc[iO2]  = beta - (x+y/4.0);

    }
    else{                                         // rich: burn all O2, leftover fuel

        //double eta = beta/(x+y/4.0); // extent of reaction
        double eta = (beta-nOnotFromO2/2)/(x+y/4-nOnotFromO2/2); // extent of reaction
        if(eta > 1.0)
            *odtl->io->ostrm << endl << "eta > 1.0" << endl;
        d1 = 1.0-eta;                            // fraction of fuel unburnt

        for(int k=0; k<nspc; k++)
            ypcc[k] *= d1;                       // initialize everything then correct
        if(iCO2 > 0)                             // CO2 is not in the H2 mechs
            ypcc[iCO2] = (x-nCnotFromCO2) + nCnotFromCO2*eta;       // init + formed
        ypcc[iH2O] = (y-nHnotFromH2O)*0.5 + nHnotFromH2O*0.5*eta;   // init + formed
        ypcc[iN2]  = z*0.5;
        ypcc[iO2]  = 0.0;

    }

    //--------------------------------------------------------------------------

    double sum = 0.0;                       // normalize moles
    for(int k=0; k<nspc; k++)
        sum += ypcc[k];
    assert(sum != 0.0);
    for(int k=0; k<nspc; k++)
        ypcc[k] /= sum;
    odtl->gas->setMoleFractions( &ypcc[0] );      // set mole fractions
    odtl->gas->getMassFractions( &ypcc[0] );      // set ypcc as mass fractions


    //--------------------------------------------------------------------------

    odtl->gas->setState_HP(hpcc, odtl->odtp->pres, 1.E-10);    // get temperature as Tadiabatic
    //odtl->gas->setState_HP(hpcc, odtl->odtp->pres);    // get temperature as Tadiabatic
    Tpcc = odtl->gas->temperature();

}

///////////////////////////////////////////////////////////////////////////////
/** Set the stoichiometric mixture fraction using Bilger's definition */

void streams::setStoicMixf() {

    vector<double> x_c(nspc,0.0);

    double mc0, mc1, mo0, mo1, mh0, mh1;

    vector<double> elemMassFrac0 = setElementMassFracs(&y0[0]);
    vector<double> elemMassFrac1 = setElementMassFracs(&y1[0]);

    mc0 = elemMassFrac0[odtl->gas->elementIndex("C")]/
        odtl->gas->atomicWeight(odtl->gas->elementIndex("C"));
    mc1 = elemMassFrac1[odtl->gas->elementIndex("C")]/
        odtl->gas->atomicWeight(odtl->gas->elementIndex("C"));
    mh0 = elemMassFrac0[odtl->gas->elementIndex("H")]/
        odtl->gas->atomicWeight(odtl->gas->elementIndex("H"));
    mh1 = elemMassFrac1[odtl->gas->elementIndex("H")]/
        odtl->gas->atomicWeight(odtl->gas->elementIndex("H"));
    mo0 = elemMassFrac0[odtl->gas->elementIndex("O")]/
        odtl->gas->atomicWeight(odtl->gas->elementIndex("O"));
    mo1 = elemMassFrac1[odtl->gas->elementIndex("O")]/
        odtl->gas->atomicWeight(odtl->gas->elementIndex("O"));

    mixfStoic = (2.0*mc0 + 0.5*mh0 - mo0) /
        (mo1-mo0 + 2.0*(mc0-mc1) + 0.5*(mh0-mh1));

    *odtl->io->ostrm << endl << "mixfStoic = m_fuel/(m_fuel+m_air) = " << mixfStoic << endl;

}

///////////////////////////////////////////////////////////////////////////////
/** Sets the elements to have the correct Mass Fractions based on the specified array.
 *  @param y \input mass fraction array to use to get corresponding element fractions.
 *  @return vector of element mass fractions.
 */

vector<double> streams::setElementMassFracs(const double *y) {


    vector<double> atomArr(odtl->gas->nElements());
    vector<double> elemMassFrac(odtl->gas->nElements(), 0.0);
    double sum = 0.0;

    odtl->gas->setMassFractions( &y[0] );

    for(int k=0; k<nspc; k++) {
        sum=0.0;
        odtl->gas->getAtoms(k, &atomArr[0]);           // [nelements] in sp k
        for(int m=0; m<(int)atomArr.size(); m++)
            sum += atomArr[m] * odtl->gas->atomicWeight(m);
        for(int m=0; m<(int)atomArr.size(); m++)
            elemMassFrac[m] += y[k] * atomArr[m]/sum * odtl->gas->atomicWeight(m);
                              // is * mass frac of elem in sp
    }

    return elemMassFrac;

}

///////////////////////////////////////////////////////////////////////////////
/** Sets the elements to have the correct Mole Fractions based on the specified array.
 *  @param y \input mass fraction array to use to get corresponding element fractions.
 *  @return vector of element mole fractions.
 */

vector<double> streams::setElementMoleFracs(const double *y) {


    vector<double> atomArr(odtl->gas->nElements());
    vector<double> elemMoleFrac(odtl->gas->nElements(), 0.0);

    odtl->gas->setMassFractions( &y[0] );

    for(int k=0; k<nspc; k++) {
        odtl->gas->getAtoms(k, &atomArr[0]);           // [nelements] in sp k
        for(int m=0; m<(int)atomArr.size(); m++)
            elemMoleFrac[m] += odtl->gas->moleFraction(k) * atomArr[m];
    }
    double sum = 0.0;
    for(int m=0; m<(int)atomArr.size(); m++)
        sum += elemMoleFrac[m];
    assert(sum != 0.0);
    for(int m=0; m<(int)atomArr.size(); m++)
        elemMoleFrac[m] /= sum;

    return elemMoleFrac;

}

///////////////////////////////////////////////////////////////////////////////
/** Get amount of moles for each element.
 *  @param x \input pointer to vector of species mole fractions.
 *  @param nOnotFromO2 \input number of moles of oxygen not from O2 (oxygen in the base fuel).
 *  @param nHnotFromH2O \input number of moles of hydrogen not from H2O.
 *  @param nCnotFromCO2 \input number of moles of carbon not from CO2.
 *  @return vector of element moles.
 */

vector<double> streams::getElementMoles(const double *x,
                                        double &nOnotFromO2,
                                        double &nHnotFromH2O,
                                        double &nCnotFromCO2) {


    vector<double> atomArr(odtl->gas->nElements());
    vector<double> elemM(odtl->gas->nElements(), 0.0);
    int iO2  = odtl->gas->speciesIndex("O2");
    int iO   = odtl->gas->elementIndex("O");
    int iCO2 = odtl->gas->speciesIndex("CO2");
    int iC   = odtl->gas->elementIndex("C");
    int iH2O = odtl->gas->speciesIndex("H2O");
    int iH   = odtl->gas->elementIndex("H");

    for(int k=0; k<nspc; k++) {
        odtl->gas->getAtoms(k, &atomArr[0]);           // [nelements] in sp k
        for(int m=0; m<(int)atomArr.size(); m++)
            elemM[m] += x[k] * atomArr[m];
        if(k != iO2)  nOnotFromO2  += atomArr[iO] * x[k];
        if(k != iCO2) nCnotFromCO2 += atomArr[iC] * x[k];
        if(k != iH2O) nHnotFromH2O += atomArr[iH] * x[k];
    }
    return elemM;

}

///////////////////////////////////////////////////////////////////////////////
/**Compute the mixture fraction from the mass fractions using Bilger's mixf.
 * Set doBeta01=true on first call to initialize members beta0, beta1.
 * Later calls of this function only use the first parameter.
 *
 * @param y \input vector of species mass fractions.
 * @param doBeta01 \input flag=true on first call to set members beta0, beta1.
 * @return mixture fraction
 */

double streams::getMixtureFraction(const double *y, const bool doBeta01) {

    vector<double> elemMF;
    double         beta;

    elemMF = setElementMassFracs(y);
    beta   = gCHON[0] * elemMF[odtl->gas->elementIndex("C")] +
             gCHON[1] * elemMF[odtl->gas->elementIndex("H")] +
             gCHON[2] * elemMF[odtl->gas->elementIndex("O")] +
             gCHON[3] * elemMF[odtl->gas->elementIndex("N")];

    if(doBeta01) {
        elemMF = setElementMassFracs(&y0[0]);
        beta0  = gCHON[0] * elemMF[odtl->gas->elementIndex("C")] +
                 gCHON[1] * elemMF[odtl->gas->elementIndex("H")] +
                 gCHON[2] * elemMF[odtl->gas->elementIndex("O")] +
                 gCHON[3] * elemMF[odtl->gas->elementIndex("N")];

        elemMF = setElementMassFracs(&y1[0]);
        beta1  = gCHON[0] * elemMF[odtl->gas->elementIndex("C")] +
                 gCHON[1] * elemMF[odtl->gas->elementIndex("H")] +
                 gCHON[2] * elemMF[odtl->gas->elementIndex("O")] +
                 gCHON[3] * elemMF[odtl->gas->elementIndex("N")];
    }

    if( beta1-beta0 == 0 )
        return 0.0;            // to avoid division by zero
    else
        return (beta - beta0)/(beta1 - beta0);

}

