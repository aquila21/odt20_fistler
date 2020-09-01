
#include "micromixer.h"
#include "odtline.h"
#include "particle.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <numeric> //accumulate

#include "interp_linear.h"

///////////////////////////////////////////////////////////////////////////////
/** micromixer constructor function
 */

micromixer::micromixer() {
    cvode = new cvodeDriver();
}

///////////////////////////////////////////////////////////////////////////////
/** micromixer initialization function
 *
 * @param p_odtl  \input set odtline pointer with.
 */

void micromixer::init(odtline *p_odtl) {
    odtl = p_odtl;

    bool LincludeRhsMix = (odtl->odtp->Lsolver=="SEMI-IMPLICIT") ? true : false;
    cvode->init(odtl, LincludeRhsMix);
}

///////////////////////////////////////////////////////////////////////////////
/** Advance ODT solution: diffusion and reaction
  */

void micromixer::advanceOdt(const double p_tstart, const double p_tend) {

    tstart = p_tstart;
    tend   = p_tend;

    setNominalStepSize();
    if(odtl->odtp->LdoDL) do_DL("init");

	odtl->solv->initCstats();

    for(time=tstart; time<tend; time+=dt) {

        if(adaptGridsIfNeeded())
           setNominalStepSize();
        setStepSize();
        if(odtl->odtp->Lsolver=="EXPLICIT")
            advanceOdtSingleStep_Explicit();
        else if(odtl->odtp->Lsolver=="SEMI-IMPLICIT")
            advanceOdtSingleStep_SemiImplicit();
        else if(odtl->odtp->Lsolver=="STRANG")
            advanceOdtSingleStep_StrangSplit();

        odtl->io->dumpLineIfNeeded();
		odtl->solv->statsTime(dt);
    }
	odtl->solv->cstats2statGrd();
    if(odtl->odtp->LdoDL) do_DL("calc a");

}
///////////////////////////////////////////////////////////////////////////////
/**Set the cell sizes vectors: dxc and dx */

void micromixer::setGridDxcDx() {

    odtl->mesher->setGridDxc(odtl, dxc, odtl->odtp->cCoord);
    odtl->mesher->setGridDx(odtl, dx);

}

///////////////////////////////////////////////////////////////////////////////
/** Set time step size.
 *  This is based on a diffusive (or other) timescale.
 *  This is a uniform step size for the given integration period.
 *  The actual step size will be rest based on a data dump or tend.
 */

void micromixer::setNominalStepSize() {

    if (odtl->odtp->Lspatial) {
        double velMin = *min_element(odtl->uvel->d.begin(), odtl->uvel->d.end());
        if (velMin <= 0.0) {
            cout << "\nError micromixer::setNominalStepSize: velMin = " << velMin << ": neg or 0" << endl;
            exit(0);
        }
    }

    odtl->mesher->setGridDx(odtl, dx);

    double coef = 0.0;
    double dmb;
    for (int i=0; i < odtl->ngrd; i++) {
        dmb = odtl->dvisc->d.at(i) / odtl->rho->d.at(i) / dx.at(i) / dx.at(i);
        if (odtl->odtp->Lspatial)
            dmb /= odtl->uvel->d.at(i);
        if (dmb > coef)
            coef = dmb;
    }
    dtStepNominal = odtl->odtp->diffCFL * 0.5 / coef;
}

///////////////////////////////////////////////////////////////////////////////
/** Set time step size.
 *  This is based on a diffusive (or other) timescale.
 *  This is a uniform step size for the given integration period.
 *  The actual step size will be rest based on a data dump or tend.
 */

void micromixer::setStepSize() {

    dt = dtStepNominal;

    if (time+dt > odtl->io->dumpTimes.at(odtl->io->iNextDumpTime)) {
        dt = odtl->io->dumpTimes.at(odtl->io->iNextDumpTime) - time;
        odtl->io->LdoDump = true;
    }

    if (time + dt > tend) {
        dt = (tend - time)*(1.0+1.0E-8);
        odtl->io->LdoDump = false;
    }
}

///////////////////////////////////////////////////////////////////////////////
/** Advance ODT solution: diffusion and reaction
 */

void micromixer::advanceOdtSingleStep_Explicit(){

    setGridDxcDx();
    setGf();
    if(odtl->odtp->LdoDL) do_DL("set DL_1");

    odtl->odtc->setCaseSpecificVars();

    set_oldrho_or_rhov();
    if(odtl->odtp->Lspatial) transform(oldrho_or_rhov.begin(), oldrho_or_rhov.end(), odtl->uvel->d.begin(), oldrho_or_rhov.begin(), multiplies<double>());

    for(int k=0; k<odtl->v.size(); k++)
        if(odtl->v.at(k)->L_transported) {
            odtl->v.at(k)->getRhsMix(gf, dxc);
            odtl->v.at(k)->getRhsSrc();
        }

	//------------ particle advancement MF -------------------                                          
	if(odtl->part->particleON) 
		odtl->part->advanceParticleAndCalcWeights(dxc, dt);                                       
	//--------------------------------------------------------   

    for(int k=0; k<odtl->v.size(); k++){
        if(odtl->v.at(k)->L_transported){
            for(int i=0; i < odtl->ngrd; i++){
                if(odtl->odtp->partCoupl == 2 && odtl->part->particleON){
                    if(odtl->v.at(k)->var_name=="uvel"|| odtl->v.at(k)->var_name=="wvel"){
                        if(odtl->part->sumWeights.at(i) != 0.0){
                            vector<double> vel0 = odtl->v.at(k)->d; 
                            int fac = ceil(odtl->part->sumWeights.at(i)*dt*2); //check timestep and split if necessary
                            for(int j=0; j<fac; j++){
                                if(odtl->v.at(k)->var_name=="uvel"){ 
                                    odtl->v.at(k)->d.at(i) = odtl->v.at(k)->d.at(i) + dt/fac*( odtl->v.at(k)->rhsMix.at(i) + odtl->v.at(k)->rhsSrc.at(i) + odtl->part->sumWeights.at(i)*(odtl->part->av_uP.at(i)-odtl->v.at(k)->d.at(i)));
                                }
                                else{
                                    odtl->v.at(k)->d.at(i) = odtl->v.at(k)->d.at(i) + dt/fac*( odtl->v.at(k)->rhsMix.at(i) + odtl->v.at(k)->rhsSrc.at(i) + odtl->part->sumWeights.at(i)*((odtl->part->av_wP.at(i)-odtl->v.at(k)->d.at(i))));
                                }
                            }
                        }
                        else{
                            odtl->v.at(k)->d.at(i) = odtl->v.at(k)->d.at(i) + dt*( odtl->v.at(k)->rhsMix.at(i) + odtl->v.at(k)->rhsSrc.at(i));
                        }
                    }
                    else{
                        odtl->v.at(k)->d.at(i) = odtl->v.at(k)->d.at(i) + dt*( odtl->v.at(k)->rhsMix.at(i) + odtl->v.at(k)->rhsSrc.at(i));
                    }
                }
                else{
                    odtl->v.at(k)->d.at(i) = odtl->v.at(k)->d.at(i) + dt*( odtl->v.at(k)->rhsMix.at(i) + odtl->v.at(k)->rhsSrc.at(i));
                }
            }
        }
    }
    if(odtl->odtp->probType == "SHEARFLOW"){
        if(abs(odtl->uvel->d.at(odtl->ngrd-1) - odtl->uvel->d.at(0) - odtl->odtp->uBChi+odtl->odtp->uBClo) > 0.1){
            odtl->uvel->d.at(odtl->ngrd-1) = odtl->odtp->uBChi;
            odtl->uvel->d.at(0) = odtl->odtp->uBClo;
        }
    }

    updateGrid();            // update cell sizes due to rho or rho*v variations (continuity)

    if(odtl->odtp->LdoDL) do_DL("set DL_2");

    odtl->mesher->enforceDomainSize();     // chop the domain

}

///////////////////////////////////////////////////////////////////////////////
/** Advance ODT solution: diffusion and reaction; Some terms are implicit, others explicit.
 *  Nominally the mixing terms are explicit. Calling the cvode driver.
 *  First order.
 *  dphi/dt =  D(phi_0) + S(phi) : solving from t0 to t1.
 *  Here, D() is the diffusive term, and S() is the (stiff) source term.
 *  We solve the whole RHS implicitly, but the D(phi_0) is fixed at time 0.
 */

void micromixer::advanceOdtSingleStep_SemiImplicit() {

    if(odtl->odtp->Lsolver!="SEMI-IMPLICIT")
        return;

    setGridDxcDx();
    setGf();
    odtl->odtc->setCaseSpecificVars();
    set_oldrho_or_rhov();
    if(odtl->odtp->Lspatial) transform(oldrho_or_rhov.begin(), oldrho_or_rhov.end(), odtl->uvel->d.begin(), oldrho_or_rhov.begin(), multiplies<double>());

    //--------------- Set the explicit (mixing) terms

    for(int k=0; k<odtl->v.size(); k++)
        if(odtl->v.at(k)->L_transported)
            odtl->v.at(k)->getRhsMix(gf, dxc);

    //--------------- Perform the implicit integration on each cell

    for(int i=0; i<odtl->ngrd; i++)
        cvode->integrateCell(i, dt);

    //---------------

    updateGrid();            // update cell sizes due to density variations (continuity)

    odtl->mesher->enforceDomainSize();     // chop the domain

}

///////////////////////////////////////////////////////////////////////////////
/** Advance ODT solution: diffusion and reaction; Some terms are implicit, others explicit.
 *  Nominally the mixing terms are explicit. Calling the cvode driver.
 *  First order.
 *  dphi/dt =  D(phi_0) + S(phi) : solving from t0 to t1.
 *  Here, D() is the diffusive term, and S() is the (stiff) source term.
 *  We solve the whole RHS implicitly, but the D(phi_0) is fixed at time 0.
 */

void micromixer::advanceOdtSingleStep_StrangSplit() {

    if(odtl->odtp->Lsolver!="STRANG")
        return;

    //--------------- First step: phi_1 = phi_0 + 0.5*dt*D(phi_0)

    setGridDxcDx();
    setGf();
    odtl->odtc->setCaseSpecificVars();
    set_oldrho_or_rhov();
    if(odtl->odtp->Lspatial) transform(oldrho_or_rhov.begin(), oldrho_or_rhov.end(), odtl->uvel->d.begin(), oldrho_or_rhov.begin(), multiplies<double>());

    for(int k=0; k<odtl->v.size(); k++)
        if(odtl->v.at(k)->L_transported)
            odtl->v.at(k)->getRhsMix(gf, dxc);

    for(int k=0; k<odtl->v.size(); k++)
        if(odtl->v.at(k)->L_transported)
            for(int i=0; i < odtl->ngrd; i++)
                odtl->v.at(k)->d.at(i) = odtl->v.at(k)->d.at(i) + 0.5*dt*odtl->v.at(k)->rhsMix.at(i);

    updateGrid();            // update cell sizes due to density variations (continuity)

    //--------------- Second step: phi_2 = phi_1 + dt*S(phi_2)
    // (actually: dphi/dt = S(phi), with initial condition phi=phi_1. Implicit.)

    setGridDxcDx();
    //not needed: setGf();
    odtl->odtc->setCaseSpecificVars();
    set_oldrho_or_rhov();
    if(odtl->odtp->Lspatial) transform(oldrho_or_rhov.begin(), oldrho_or_rhov.end(), odtl->uvel->d.begin(), oldrho_or_rhov.begin(), multiplies<double>());

    for(int k=0; k<odtl->v.size(); k++)
        if(odtl->v.at(k)->L_transported)
            odtl->v.at(k)->getRhsSrc();

    for(int i=0; i<odtl->ngrd; i++)
        cvode->integrateCell(i, dt);

    updateGrid();            // update cell sizes due to density variations (continuity)

    //--------------- Third step: phi_3 = phi_2 + 0.5*dt*D(phi_2)

    setGridDxcDx();
    setGf();
    odtl->odtc->setCaseSpecificVars();
    set_oldrho_or_rhov();
    if(odtl->odtp->Lspatial) transform(oldrho_or_rhov.begin(), oldrho_or_rhov.end(), odtl->uvel->d.begin(), oldrho_or_rhov.begin(), multiplies<double>());

    for(int k=0; k<odtl->v.size(); k++)
        if(odtl->v.at(k)->L_transported)
            odtl->v.at(k)->getRhsMix(gf, dxc);

    for(int k=0; k<odtl->v.size(); k++)
        if(odtl->v.at(k)->L_transported)
            for(int i=0; i < odtl->ngrd; i++)
                odtl->v.at(k)->d.at(i) = odtl->v.at(k)->d.at(i) + 0.5*dt*odtl->v.at(k)->rhsMix.at(i);

    updateGrid();            // update cell sizes due to density variations (continuity)

    //-------------------------

    odtl->mesher->enforceDomainSize();     // chop the domain (for strang, just at the final step)

}

///////////////////////////////////////////////////////////////////////////////
/** Set the grid factor array
 */

void micromixer::setGf(){

    gf.resize(odtl->ngrdf, 0.0);

    for (int i=1, im=0; i<odtl->ngrd; i++, im++) // interior
        gf.at(i) = 2.0 / (dx.at(im) + dx.at(i));

    gf.at(0)          = 2.0 / dx.at(0);                // lo boundary
    gf.at(odtl->ngrd) = 2.0 / dx.at(odtl->ngrd - 1);   // hi boundary
    if (odtl->odtp->bcType == "PERIODIC") {            // periodic
        gf.at(0)          = 2.0 / (dx.at(0) + dx.at(odtl->ngrd - 1));
        gf.at(odtl->ngrd) = gf.at(0);
    }

}

///////////////////////////////////////////////////////////////////////////////
/** Update grid for gas expansion
 */

void micromixer::updateGrid() {

    if(odtl->odtp->bcType=="WALL" || odtl->odtp->LisFlmlt || odtl->odtp->LisHips)
        return;

    //--------------

    odtl->rho->setVar();

    vector<double> dxc2(odtl->ngrd);
    for(int i=0; i<odtl->ngrd; i++)
        dxc2[i] = dxc.at(i)*(oldrho_or_rhov.at(i)/odtl->rho->d.at(i));
    if(odtl->odtp->Lspatial)
        for(int i=0; i<odtl->ngrd; i++)
            dxc2[i] /= odtl->uvel->d.at(i);

    double C    = odtl->odtp->cCoord;
    double invC = 1.0/odtl->odtp->cCoord;

    //-------------

    odtl->mesher->setGridFromDxc(dxc2);


}

///////////////////////////////////////////////////////////////////////////////

/** Adapt during diffusion for spatial cases for which grid contraction
 *  results in small grid cells
 */

bool micromixer::adaptGridsIfNeeded() {

    if (odtl->odtp->Lspatial && *min_element(dx.begin(), dx.end()) < 0.9*odtl->odtp->dxmin) {
        //*odtl->io->ostrm << endl << "#------- ADAPTING DURING DIFFUSION";
        odtl->mesher->adaptGrid(0, odtl->ngrd-1);
        return true;
    }
    return false;
}

///////////////////////////////////////////////////////////////////////////////
/** Processes the DL instability
 *  @param doWhat \input string that indicates what to do.
 * */

void micromixer::do_DL(string doWhat) {

    if(!odtl->odtp->LdoDL)
        return;

    //--------------------------------------------------
    if(doWhat == "init") {

        uDL_1     = vector<double>(odtl->ngrd, 0.0);
        uDL_2     = vector<double>(odtl->ngrd, 0.0);
        xDL_1     = odtl->pos->d;
        xDL_2     = odtl->pos->d;
        posDL_old = odtl->pos->d;
        odtl->aDL->d = vector<double>(odtl->ngrd, 0.0);
    }
    //--------------------------------------------------
    else if(doWhat == "set DL_1") {

        xDL_1 = xDL_2;
        uDL_1 = uDL_2;
        posDL_old = odtl->pos->d;

    }
    //--------------------------------------------------
    else if(doWhat == "set DL_2") {

        xDL_2.resize(odtl->ngrd);
        uDL_2.resize(odtl->ngrd);

        for(int i=0; i<odtl->ngrd; i++) {
            uDL_2.at(i) = (odtl->pos->d[i] - posDL_old.at(i)) / (odtl->odtp->Lspatial ? dt/odtl->uvel->d[i] : dt);
            xDL_2.at(i) = 0.5*(posDL_old.at(i) + odtl->pos->d[i]);
        }
    }
    //--------------------------------------------------
    else if(doWhat == "calc a") {

        vector<double> dmb;

        dmb = uDL_1;
        Linear_interp Linterp(xDL_1, dmb);
        uDL_1.resize(odtl->ngrd);
        for(int i=0; i<odtl->ngrd; i++)
            uDL_1.at(i)= Linterp.interp(odtl->pos->d[i]);

        dmb = uDL_2;
        Linterp = Linear_interp(xDL_2, dmb);
        uDL_2.resize(odtl->ngrd);
        for(int i=0; i<odtl->ngrd; i++)
            uDL_2.at(i)= Linterp.interp(odtl->pos->d[i]);

        for(int i=0; i<odtl->ngrd; i++)
            odtl->aDL->d.at(i) = (uDL_2.at(i) - uDL_1.at(i)) / (odtl->odtp->Lspatial ? dt/odtl->uvel->d[i] : dt);
    }
}

///////////////////////////////////////////////////////////////////////////////
/** Store the old density for continuity.
  * This is a function for generality on inheritance.
  */

void micromixer::set_oldrho_or_rhov() {
    oldrho_or_rhov = odtl->rho->d;
}


#include <iomanip>
void micromixer::check_balance(int io) {

    setGridDxcDx();
    double mom = 0.0;

    for(int i=0; i<odtl->ngrd; i++)
        mom += odtl->rho->d[i] * odtl->uvel->d[i] * (odtl->uvel->d[i]-0.0) * dxc[i];

    cout << scientific;
    cout << setprecision(13);
    cout << endl << "check: " << io << " " << mom;
}

