
#include "solver.h"
#include "odtline.h"
#include "particle.h"
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include <numeric> //accumulate

///////////////////////////////////////////////////////////////////////////////
/** solver initialization function
 *
 * @param p_odtl  \input set odtline pointer with.
 */

void solver::init(odtline *p_odtl) {

    odtl       = p_odtl;

    if(odtl->odtp->LES_type=="THIRDS"){
        ed3   = new eddy;
        eddl3 = new odtline(odtl, odtl->odtp);
        eddl3->init(NULL,NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, true);
        ed3->init(odtl, eddl3);
    }

}

///////////////////////////////////////////////////////////////////////////////
/** Destructor
 */

solver::~solver(){
    if(ed3)   delete ed3;
    if(eddl3) delete eddl3;
}

///////////////////////////////////////////////////////////////////////////////
/** The actual solver
 *
 * Advance in time only sampling eddies till get one (EE).  Then diffuse the system
 * until you catch up to the eddy time:
 * <code><pre>
 *     last EE    this EE       diffuse to catch up
 * .....|..........|           ---------------------->  ................|
 *      t0         time                                                t0,time
 *                          (last EE)   this EE         next EE
 * Then advance some more ......;.........|................|         etc.
 *                                        t0              time
 * </pre></code>
 * Eddy timesteps are smaller than required diffusive steps, so if we were to
 * lock-step the two processes we would do too much diffusive work (i.e., take
 * more diffusive steps than needed)
 *
 */

void solver::calculateSolution() {


    //-------------------------------------------------------------------------

    double tLastDA         = 0.0;            ///< time of last diffusive mesh adaption
    int    cLastDA         = 0;              ///< for adaption
    bool   LeddyAccepted   = false;

    stringstream ss1;
    string       s1;

    time     = odtl->odtp->trst;    // 0.0 unless you are restarting
    t0       = odtl->odtp->trst;    // 0.0 unless you are restarting

    PaSum    = 0.0;
    nPaSum   = 0;
    neddies  = 0;
    PaSumC   = 0.0;
    nPaSumC  = 0;

    iEtrials = 0;

    //-------------------------------------------------------------------------

    computeDtSmean();
    computeDtCUmax();

    odtl->io->writeDataFile("odt_init.dat", time);

    initStats(); 

    odtl->mesher->adaptGrid(0, odtl->ngrd-1);

    odtl->io->writeDataFile("odt_init_adpt.dat", time);

    odtl->io->outputHeader();

    //-------------------------------------------------------------------------

    while(time <= odtl->odtp->tEnd) {
		
	    statsSetOld();
        diffusionCatchUpIfNeeded();
	    statsChange(1);

	    statsSetOld();
        odtl->mesher->adaptAfterSufficientDiffTime(time, tLastDA, cLastDA, dtCUmax);
	    statsChange(3);
	
	    computeDtCUmax();

	    statsSetOld();        
        LeddyAccepted = sampleEddyAndImplementIfAccepted();  // may reduce dtSmean; sets Pa

        iEtrials++;

        //----------------------

	    if(LeddyAccepted) {
	        statsChange(0);

		    if(odtl->odtp->partCoupl == 1 && odtl->part->particleON)                                                          
			    odtl->part->oneWayCoupl(time); // handling one way coupling

            if(++neddies % (odtl->odtp->modDisp*50) == 0)
          	    odtl->io->outputHeader();
            if(neddies   % odtl->odtp->modDisp ==0)
           	    odtl->io->outputProgress();

		    statsSetOld();
            odtl->mesher->adaptEddyRegionOfMesh(time, tLastDA, cLastDA);
		    statsChange(3);

		    statsSetOld();
            diffusionCatchUpIfNeeded(true);
		    statsChange(1);
		    //statsChange(2);

		    statsSetOld();
            odtl->mesher->adaptEddyRegionOfMesh(time, tLastDA, cLastDA);
		    statsChange(3);

	 	    if (neddies % odtl->odtp->modDump == 0) {
         		ss1.clear();  ss1 << setfill('0') << setw(4) << neddies; ss1 >> s1;
               	odtl->io->writeDataFile("odt_"+s1+"_adptDif.dat", time);
                odtl->part->writeData("part_eddy_"+s1+".dat", time);
		    }
        }
        if(odtl->odtp->kernEv){
	        if(time > odtl->kern->injectionTime){ // external TKE injection by kernel event (for forcing HIT) 	
	            statsSetOld();
	            odtl->kern->kernelEnergyInjection(time); 	
	            statsChange(4);

                odtl->io->outputKernHeader();
                odtl->io->outputKernProgress();
            
	            statsSetOld();
	            odtl->mesher->adaptEddyRegionOfMesh(time, tLastDA, cLastDA);
	            statsChange(3);

	            statsSetOld();
                diffusionCatchUpIfNeeded(true);
                statsChange(1);

                statsSetOld();
                odtl->mesher->adaptEddyRegionOfMesh(time, tLastDA, cLastDA);
	            statsChange(3);
            }
   	    }

        //----------------------

        time += sampleDt();             // advance the time
        raiseDtSmean();                 // may reset PaSum, nPaSum, dtSmean
    }

    time = odtl->odtp->tEnd;
    if(t0 < time){
	    statsSetOld();
        diffusionCatchUpIfNeeded(true);
	    statsChange(1);
    }
    //-------------------------------------------------------------------------

    odtl->io->writeDataFile("odt_end.dat", time);
    statsOutput();

}

///////////////////////////////////////////////////////////////////////////////
/**dtSmean is computed, which is the mean eddy sample time.  The Poisson
 * process draws dt's with this mean.
 */

void solver::computeDtSmean() {

    if(odtl->odtp->Llem) ;
    //    dtSmean = 1.0/(odtl->odtp->lemRateParam * odtl->Ldomain());

    else if(!odtl->odtp->Lspatial)
        dtSmean = 0.1*odtl->odtp->Pav * odtl->Ldomain() * odtl->Ldomain() /
                  odtl->odtp->kvisc0 / odtl->ngrd / odtl->ngrd / odtl->ngrd;
    else
        dtSmean = 10.*odtl->odtp->Pav * odtl->Ldomain() / odtl->ngrd / odtl->ngrd;

}

///////////////////////////////////////////////////////////////////////////////
/**Sample the eddy trial time step with mean dtSmean.
 * @return Poisson sampled timestep.
 */

double solver::sampleDt() {

    return -dtSmean*log( max(1.0e-14, odtl->rand->getRand()) );

}

///////////////////////////////////////////////////////////////////////////////
/**Computes dtCUmax (as the name suggests).  This variable is the time
 * increment of eddy trial time advancement before we diffuse to catch
 * up to that time in the event of no eddy before that time.
 */

void solver::computeDtCUmax() {

    double dxmin = 1.0E10;
    double d1;
    for(int i=1; i<odtl->ngrdf; i++) {
        d1 = odtl->posf->d.at(i) - odtl->posf->d.at(i-1);
        if(d1 < dxmin)
            dxmin = d1;
    }
    if(!odtl->odtp->Lspatial)
        dtCUmax = odtl->odtp->tdfac * dxmin * dxmin / odtl->odtp->kvisc0;
    else
        dtCUmax = 50.0 * odtl->odtp->tdfac * dxmin;

}

///////////////////////////////////////////////////////////////////////////////
/** Diffuse the line to catch up t0 to time if we have not had eddies for a while.
 *  @param Ldoit  \input Flag with default false
 */

void solver::diffusionCatchUpIfNeeded(bool Ldoit) {

    if(!Ldoit && time-t0 < dtCUmax)
        return;

#ifndef SILENT
//    *odtl->io->ostrm << endl << "# Catching up diffusion to eddies: time = " << time;
#endif
	
    odtl->mimx->advanceOdt(t0, time);

    t0 = time;

}

///////////////////////////////////////////////////////////////////////////////
/**Every once in a while (nDtSmeanWait) test the mean acceptance probability.
 * Increase dtSmean if its too small.
 */

void solver::raiseDtSmean() {

    if(odtl->odtp->Llem)
        return;

    if(iEtrials % odtl->odtp->nDtSmeanWait != 0)
        return;                              // only do this once in a while

    if(nPaSum > 0)
        PaSum /= nPaSum;                     // PaSum is now average Pa of eddies

    if(PaSum < odtl->odtp->Pav) {
       if(PaSum < odtl->odtp->Pav/odtl->odtp->dtfac)
           dtSmean *= odtl->odtp->dtfac;     // increase by at most dtfac
       else
           dtSmean *= odtl->odtp->Pav/PaSum; // increase dtSmean to target value

       *odtl->io->ostrm << endl << "#  increasing dtSmean to " << dtSmean
                        << " (neddies = " << neddies << "  eddy trails = " << iEtrials << ")";
    }

    PaSum  = 0.0;                          // reset values
    nPaSum = 0;
    if (iEtrials > 10000*odtl->odtp->nDtSmeanWait) {
        *odtl->io->ostrm << endl << "#  reset iEtrials, PaSumC, nPaSumC after "
                         << 10000*odtl->odtp->nDtSmeanWait << " eddy trials.";
        iEtrials = 0;
        PaSumC   = 0.0;
        nPaSumC  = 0;
    }

}

///////////////////////////////////////////////////////////////////////////////
/**Sample an eddy size and position. Fill the eddl from odtl.
 * Then triplet map the eddy, compute the eddy timescale, then the
 * acceptance probability.  Roll the dice and if you win (rand # < prob)
 * then accept the eddy.  This means, apply velocity kernels, then
 * insert the eddl into the odtl.
 * Note, this function may be better as a member of eddy
 *
 * @return true if the sampled eddy was implemented.
 */

bool solver::sampleEddyAndImplementIfAccepted() {

	stringstream ss1;
    string       s1;

    odtl->ed->sampleEddySize();

    if(!testLES_fracDomain(odtl->ed->eddySize) && !odtl->odtp->Llem)
        return false;

    odtl->ed->sampleEddyPosition();


    //---------- extract the eddy segment from odtl into eddl

    int iStart = odtl->linePositionToIndex(odtl->ed->leftEdge,  true, 4);
	int iEnd = odtl->linePositionToIndex(odtl->ed->rightEdge, false, 5);

    if(!odtl->ed->LperiodicEddy) {
        if(iEnd - iStart < odtl->odtp->eddyMinCells && !odtl->odtp->Llem) // only for odt because numerical shear can be >> physical shear
            return false;
    }
    else {
        if( (odtl->ngrd-iStart)+(iEnd)+1 < odtl->odtp->eddyMinCells && !odtl->odtp->Llem)
            return false;
		if( iEnd < 1)
            return false;
    }

    odtl->eddl->setLineFromRegion(iStart, iEnd);
    //---------- invoke the eddy

    odtl->ed->tripMap(odtl->eddl, 0, odtl->eddl->ngrd-1, odtl->odtp->cCoord);

    if(!odtl->odtp->Llem) {

    	if(!odtl->ed->eddyTau(odtl->odtp->Z_param, odtl->odtp->cCoord))
    		return false;

        double eddyTauOrSizeForLES = (odtl->odtp->Lspatial ? odtl->ed->eddySize : 1.0/odtl->ed->invTauEddy);
        if(!testLES_elapsedTime(time, eddyTauOrSizeForLES))
            return false;

        if(!testLES_integralLength(time, odtl->ed->eddySize))
            return false;

	//-------------------Two-way coupling particles (MF)                                            
  		if(odtl->odtp->partCoupl > 1 && odtl->part->particleON){                                                              
  			if(!odtl->part->twoWayCoupl(time, iStart, iEnd))                                       
  				return false;                                                                       
  	                                                                                       
  			eddyTauOrSizeForLES = (odtl->odtp->Lspatial ? odtl->ed->eddySize : 1.0/odtl->ed->invTauEddy);
  			if(!testLES_elapsedTime(time, eddyTauOrSizeForLES))                                     
  				return false;                                                                       
  			if(!testLES_integralLength(time, odtl->ed->eddySize))                                   
  				return false;                                                                       
  		}                                                                                           
  //--------------------------------------------------

        odtl->ed->computeEddyAcceptanceProb(dtSmean);

        //---------- lower dtSample if Pa too high; changes Pa, dtSmean

        lowerDtSmean();
    }
    //---------- apply the eddy if accepted

    double rnd = odtl->odtp->Llem ? -1 : odtl->rand->getRand();

    if(rnd <= odtl->ed->Pa) {
        if(!odtl->odtp->Llem && !testLES_thirds())  // large eddy supression test
            return false;
        ss1.clear();  ss1 << setfill('0') << setw(4) << neddies; ss1 >> s1;
        if(odtl->ed->LperiodicEddy) {
            *odtl->io->ostrm << endl << "#   periodic eddy ";
            double cycleDistance = odtl->cyclePeriodicLine(iEnd);
            double bkp_rightEdge = odtl->ed->rightEdge;
            odtl->ed->rightEdge = odtl->ed->leftEdge + odtl->ed->eddySize;
            iStart = odtl->linePositionToIndex(odtl->ed->leftEdge,  true, 6);
            iEnd   = odtl->linePositionToIndex(odtl->ed->rightEdge, false, 7);
	    	odtl->ed->tripMap(odtl, iStart, iEnd, odtl->odtp->cCoord, true);
            iEnd   = odtl->ngrd-1;;
			if(!odtl->ed->applyVelocityKernels(odtl, iStart+1, iEnd-1))
                return false;
    		odtl->backCyclePeriodicLine(cycleDistance);
            odtl->ed->rightEdge = bkp_rightEdge;
        }
        else{
            odtl->ed->tripMap(odtl, iStart, iEnd, odtl->odtp->cCoord, true);
            iStart = odtl->linePositionToIndex(odtl->ed->leftEdge,  true, 8);
            iEnd   = odtl->linePositionToIndex(odtl->ed->rightEdge, false, 9);
            if(odtl->odtp->Lspatial) {     // vary domain to conserve mass flux for u changes from kernels
                vector<double> dxc_or_rhoUdxc(odtl->ngrd);
                odtl->mesher->setGridDxc(odtl, dxc_or_rhoUdxc, odtl->odtp->cCoord);
                for(int i=0; i<odtl->ngrd; i++)
                    dxc_or_rhoUdxc[i] = dxc_or_rhoUdxc[i] * odtl->rho->d[i] * odtl->uvel->d[i];

                if(!odtl->ed->applyVelocityKernels(odtl, iStart, iEnd))
				    return false;

                for(int i=0; i<odtl->ngrd; i++)
                    dxc_or_rhoUdxc[i] = dxc_or_rhoUdxc[i] / (odtl->rho->d[i] * odtl->uvel->d[i]);
                odtl->mesher->setGridFromDxc(dxc_or_rhoUdxc);
                odtl->mesher->enforceDomainSize();     // chop the domain
            }
            else{
                if(!odtl->ed->applyVelocityKernels(odtl, iStart, iEnd))
					    return false;
			}
		}
        return true;
    }

    return false;
}

///////////////////////////////////////////////////////////////////////////////
/** Apply a large eddy suppression test based on elapsed time (Echekki 2001)
 *  @param time    \input current time.
 *  @param tauEddy \input eddy timescale, or in spatial cases, the eddy size
 *  Note, this function may be better as a member of eddy
 */

bool solver::testLES_elapsedTime(const double time, const double tauEddy) {

    if(odtl->odtp->LES_type != "ELAPSEDTIME")
        return true;

    if( (time-odtl->odtp->x0virtual) < odtl->odtp->Z_LES * tauEddy )
        return false;

    return true;
}

///////////////////////////////////////////////////////////////////////////////
/** Large eddy suppression test on fraction of domain size
 *  @param eSize \input eddy size
 *  Note, this function may be better as a member of eddy
 */
bool solver::testLES_fracDomain(const double eSize) {
    if(odtl->odtp->LES_type != "FRACDOMAIN")
        return true;
    if(eSize/odtl->Ldomain() > odtl->odtp->Z_LES)
        return false;
    return true;
}

///////////////////////////////////////////////////////////////////////////////
/** Apply a large eddy suppression test based on integral length scale
 *  integral length scale L = L0 * (t/t0)^(1-n/2)
 *  t is elapsed time;
 *  n is between 1.15 and 1.45, usually 1.3
 *  @param time  \input current time.
 *  @param eSize \input eddy size
 *  Guangyuan Sun 06/2013
 *  Note, this function may be better as a member of eddy
 */

bool solver::testLES_integralLength(const double time, const double eSize) {

    if(odtl->odtp->LES_type != "INTEGRALSCALE")
        return true;

    double n = 1.1;
    double t0 = 0.15899;
    double L0 = 0.028323;
    double integralLength = odtl->odtp->Z_LES * L0 * pow(time/t0, 1-0.5*n);

    if(eSize > integralLength)
        return false;
    return true;
}

///////////////////////////////////////////////////////////////////////////////
/** Apply the large-eddy suppression test.
 *  Note, this function may be better as a member of eddy
 */

bool solver::testLES_thirds() {

    if(odtl->odtp->LES_type != "THIRDS")
        return true;

    ed3->eddySize = odtl->ed->eddySize/3.0;
    ed3->LperiodicEddy = false;

    double leftThird = odtl->ed->leftEdge;
    double rightThird = leftThird + odtl->ed->eddySize/3.0;

    for(int third=0; third<3; third++) {

        ed3->leftEdge  = leftThird;
        ed3->rightEdge = rightThird;

        int iStart = odtl->linePositionToIndex(leftThird,  true,10);
        int iEnd   = odtl->linePositionToIndex(rightThird, false,11);

        eddl3->setLineFromRegion(iStart, iEnd);

        //---------- invoke the eddy

        ed3->tripMap(eddl3, 0, eddl3->ngrd-1, odtl->odtp->cCoord);                       // apply trip map to eddyLine

        if(!ed3->eddyTau(odtl->odtp->Z_LES, 1))
            return false;

        leftThird  = rightThird;
        rightThird = leftThird + odtl->ed->eddySize/3.0;

    }

    return true; // Eddy is allowed (not suppressed).
}

///////////////////////////////////////////////////////////////////////////////
/** Reduce dtSmean if it is resulting in too high an acceptance probability.
 */

void solver::lowerDtSmean() {

    if(odtl->ed->Pa > odtl->odtp->Pmax) {
        dtSmean *= odtl->odtp->Pmax/odtl->ed->Pa;
        odtl->ed->Pa = odtl->odtp->Pmax;
        *odtl->io->ostrm << endl << "#   reducing dtSmean to " << dtSmean;
        *odtl->io->ostrm << endl << "#     ed.Pa, odtp.Pmax = " << odtl->ed->Pa << " " << odtl->odtp->Pmax;
        *odtl->io->ostrm << endl << "#     time " << time;
    }

    //---------- update properties used to raise dtSample

    PaSum += odtl->ed->Pa;
    nPaSum++;

    PaSumC += odtl->ed->Pa;
    nPaSumC++;

}

///////////////////////////////////////////////////////////////////////////////
/** Initialize stats grid and stats vectors.
 */

void solver::initStats() {
	
	iStat = 1.0;
	LdomainS = odtl->odtp->domainLength;
    ngrdS    = odtl->odtp->ngrdStat;
    ngrdfS   = ngrdS+1;
	dtStat  = odtl->odtp->tEnd/odtl->odtp->nStat;
    posS     = vector<double>(ngrdS,  0.0);
    posfS    = vector<double>(ngrdfS, 0.0);

    uMeanS   = vector<double>(ngrdS,  0.0);
    vMeanS   = vector<double>(ngrdS,  0.0);
    wMeanS   = vector<double>(ngrdS,  0.0);

    vTrans   = vector<double>(ngrdS, 0.0);

    double dp = LdomainS/ngrdS;
 	posfS.at(0) = -1*LdomainS/2;
    for(int i=1; i<ngrdfS; i++)
        posfS.at(i) = posfS.at(i-1) + dp;
    posS.at(0) = posfS.at(0) + dp/2;
    for(int i=1; i<ngrdS; i++)
        posS.at(i) = posS.at(i-1) + dp;

    // edstat(j,k,l,i): j = 0.5*del_U_eddy^2, 0.5*del_U_flow^2, 0.5*del_U_all^2, 0.5*del_U_adpt^2, 0.5*del_U_kern^2
    //                  k = uv, u2, v2, w2 
    //                  l = average period
    //                  i = coordinate
    edstat      = vector<vector<vector<vector<double> > > >
                    (5, vector<vector<vector<double> > >
                    (4, vector<vector<double> >
                    (odtl->odtp->nStat, vector<double>(ngrdS, 0.0) ) ) );
    // cstat(m,k,l,i):  i = coordinate
    //                  m = u, u^2, du   
    //                  k = u, v, w, 
    //                  l = average period
    // ctime(l):        l = average period
    cstat       = vector<vector<vector<vector<double> > > >
                    (3, vector<vector<vector<double> > >
                    (3, vector<vector<double> >
                    (odtl->odtp->nStat, vector<double>(ngrdS, 0.0) ) ) );
    ctime       = vector<double>(odtl->odtp->nStat, 0.0);
    // oldVars(k,i)     i = coordinate
    //                  k = u, u2, v2, w2
    oldVars     = vector<vector<double> > (4, vector<double>(ngrdS, 0.0) );

	for(int i=0; i<ngrdS; i++){
        uMeanS.at(i) = 0.0; //odtl->io->initParams["Srate"].as<double>()*posfS.at(i);
		vMeanS.at(i) = 0.0;
		wMeanS.at(i) = 0.0;
    }
	if(odtl->odtp->partCoupl > 0) odtl->part->initStats();

}

///////////////////////////////////////////////////////////////////////////////
/** Initialize stats grid and stats vectors.
 */

void solver::odtGrd2statGrd(vector<double> odtposf, vector<double> odtvec){

	int odtngrd = odtl->ngrd;
	
	//---------- Now transfer grids

	int i, ip;
	int j, jp;
    double pLeft;
	
	//---------- transfer grids

    j     = 0;
    jp    = 1;
    pLeft = posfS.at(0);

    vTrans = vector<double>(ngrdS, 0.0);
    
    for(i=0, ip=1; i<ngrdS; i++, ip++) {                     // loop over stat grd
        while( odtposf.at(jp) < posfS.at(ip) && jp!=odtngrd) {     // loop over odt grd
            vTrans.at(i) += odtvec.at(j) * (odtposf.at(jp)-pLeft);
            j++;
            jp++;
            pLeft = odtposf.at(j);
        }
        vTrans.at(i) +=  odtvec.at(j) * (posfS.at(ip) - pLeft);        // get the leftovers
        pLeft = posfS.at(ip);

    }

    //---------- normalize

    for(i=0, ip=1; i<ngrdS; i++, ip++)
        vTrans.at(i) /= (posfS.at(ip)-posfS.at(i));
}
///////////////////////////////////////////////////////////////////////////////
/**
 *  This function initializes a temporary vector of the same size as the odtline
 *  used for the data gathering during the time steps . Afterwards, the gathered
 *  date is converted to the stats-grid and added to the original data gathering.
 *  This is only done for temporal flows to leave out the conversion to a stats 
 *  grid. For spacial flows this wont work and each time step has to be converted
 *
 */
void solver::initCstats(){
	cstats = vector<vector<vector<double> > > (cstat.size(), 
                    vector<vector<double> > (cstat[0].size(), 
                    vector<double> (odtl->ngrd, 0.0) ) );
}
//////////////////////////////////////////////////////////////////////////////
/**
 *  This function switches between a temporal and a spazial simulation.
 *  Depending on the simulation case the odtline is saved in the cstats and 
 *  phstats arrays (temporal case) or the information of the odtline is 
 *  converted to the stats grid and added directly to the cstat and phstat 
 *  arrays (spazial case).
 * 
 * @param tStep \input the time step for the gathering
 *
 */
void solver::statsTime(const double tStep){

    int N = odtl->ngrd;
    double delx, delxp;

    // time t
    ctime.at(iStat-1) += tStep;
    
    if(odtl->odtp->probType != "SHEARFLOW"){
    // u velocity
        for(int i=0; i<N; i++){
            cstats[0][0][i] += tStep *     odtl->uvel->d.at(i);
            cstats[1][0][i] += tStep * pow(odtl->uvel->d.at(i),2);
            cstats[0][1][i] += tStep *     odtl->vvel->d.at(i);
            cstats[1][1][i] += tStep * pow(odtl->vvel->d.at(i),2);
            cstats[0][2][i] += tStep *     odtl->wvel->d.at(i);
            cstats[1][2][i] += tStep * pow(odtl->wvel->d.at(i),2);
        }
    }
    else{
        double shearRate = odtl->io->initParams["Srate"].as<double>();
        for(int i=0; i<N; i++){
            cstats[0][0][i] += tStep *     (odtl->uvel->d.at(i));
            cstats[1][0][i] += tStep * pow(odtl->uvel->d.at(i)-(odtl->pos->d.at(i)*shearRate),2);
            cstats[0][1][i] += tStep *     odtl->vvel->d.at(i);
            cstats[1][1][i] += tStep * pow(odtl->vvel->d.at(i),2);
            cstats[0][2][i] += tStep *     odtl->wvel->d.at(i);
            cstats[1][2][i] += tStep * pow(odtl->wvel->d.at(i),2);
        }
    }

    // derivative along y of u,v,w
    delx  = (odtl->posf->d.at(3)-odtl->posf->d.at(1))/2.0; // second distance between cell centers
    delxp = (odtl->posf->d.at(2)-odtl->posf->d.at(0))/2.0; // first distance between cell centers
    // first cell, forward approximation second order
    cstats[2][0][0] += pow( (odtl->uvel->d.at(0)-odtl->uvel->d.at(2)) * delxp / (delx*(delxp+delx))
                                    +(odtl->uvel->d.at(1)-odtl->uvel->d.at(0)) * (delxp+delx) / (delxp*delx)
                                    ,2) * tStep;
    cstats[2][1][0] += pow( (odtl->vvel->d.at(0)-odtl->vvel->d.at(2)) * delxp / (delx*(delxp+delx))
                                    +(odtl->vvel->d.at(1)-odtl->vvel->d.at(0)) * (delxp+delx) / (delxp*delx)
                                    ,2) * tStep;
    cstats[2][2][0] += pow( (odtl->wvel->d.at(0)-odtl->wvel->d.at(2)) * delxp / (delx*(delxp+delx))
                                    +(odtl->wvel->d.at(1)-odtl->wvel->d.at(0)) * (delxp+delx) / (delxp*delx)
                                    ,2) * tStep;

    for(int i=1; i<N-1; i++){
        // internal field, central approximation second order
        delx  = delxp; // copy following cell center distance to current center
        delxp = (odtl->posf->d.at(i+2)-odtl->posf->d.at(i))/2.0; // calculate next cell center distance
        cstats[2][0][i] += pow( (odtl->uvel->d.at(i+1)-odtl->uvel->d.at(i)) * delx / (delxp*(delx+delxp))
                                        +(odtl->uvel->d.at(i)-odtl->uvel->d.at(i-1)) * delxp / (delx*(delx+delxp))
                                        ,2) * tStep;
        cstats[2][1][i] += pow( (odtl->vvel->d.at(i+1)-odtl->vvel->d.at(i)) * delx / (delxp*(delx+delxp))
                                        +(odtl->vvel->d.at(i)-odtl->vvel->d.at(i-1)) * delxp / (delx*(delx+delxp))
                                        ,2) * tStep;
        cstats[2][2][i] += pow( (odtl->wvel->d.at(i+1)-odtl->wvel->d.at(i)) * delx / (delxp*(delx+delxp))
                                        +(odtl->wvel->d.at(i)-odtl->wvel->d.at(i-1)) * delxp / (delx*(delx+delxp))
                                        ,2) * tStep;
    }
    // last cell, backward approximation second order
    delx  = delxp; // copy last distance between cell centers
    delxp = (odtl->posf->d.at(N-1)-odtl->posf->d.at(N-3))/2.0; // calculate second last distance
    cstats[2][0][N-1] += pow( (odtl->uvel->d.at(N-3)-odtl->uvel->d.at(N-1)) * delx / (delxp*(delx+delxp))
                                    +(odtl->uvel->d.at(N-1)-odtl->uvel->d.at(N-2)) * (delx+delxp) / (delx*delxp)
                                    ,2) * tStep;
    cstats[2][1][N-1] += pow( (odtl->vvel->d.at(N-3)-odtl->vvel->d.at(N-1)) * delx / (delxp*(delx+delxp))
                                    +(odtl->vvel->d.at(N-1)-odtl->vvel->d.at(N-2)) * (delx+delxp) / (delx*delxp)
                                    ,2) * tStep;
    cstats[2][2][N-1] += pow( (odtl->wvel->d.at(N-3)-odtl->wvel->d.at(N-1)) * delx / (delxp*(delx+delxp))
                                    +(odtl->wvel->d.at(N-1)-odtl->wvel->d.at(N-2)) * (delx+delxp) / (delx*delxp)
                                    ,2) * tStep;
	
	if(odtl->odtp->partCoupl > 0 && odtl->part->particleON) odtl->part->statsTime(tStep);

    if(odtl->odtp->partCoupl > 0 && !odtl->part->particleON && odtl->odtp->tparticleON <= odtl->mimx->time){
        odtl->part->particleON = true;
        odtl->part->initTime = odtl->mimx->time;
        odtl->part->setVelocity();
    }

	if(ctime.at(iStat-1) > dtStat && iStat < odtl->odtp->nStat){
		cstats2statGrd();
        iStat++;
    }
    return;
}

///////////////////////////////////////////////////////////////////////////////
/**
 *  This function is only used if the simulation is temporal. In this case the 
 *  gathered data is temporally stored in cstats and phstats. These arrays are
 *  now converted to the stats grid and added to the arrays cstat and phstat.
 *
 */
void solver::cstats2statGrd(){

    int N = odtl->ngrd;
    vector<double> temp; // temporal vector for use in odtGrd2statGrd()
    temp = vector<double>(N,0.0);

    // adding cstats to cstat
    // k = 0 => uvel   k = 1 => vvel   k = 2 => wvel   
    for(int k=0; k<(int)cstats[0].size(); k++){
        for(int i=0; i<N; i++){
            temp[i] = cstats[0][k][i];
        }
        odtGrd2statGrd(odtl->posf->d, temp);
        for(int i=0; i<(int)vTrans.size(); i++){
            cstat[0][k][iStat-1][i] += vTrans[i];}
    }

    // k = 0 => uvel^2   k = 1 => vvel^2   k = 2 => wvel^2   
    for(int k=0; k<(int)cstats[0].size(); k++){
        for(int i=0; i<N; i++){
            temp[i] = cstats[1][k][i];}
        odtGrd2statGrd(odtl->posf->d, temp);
        for(int i=0; i<(int)vTrans.size(); i++){
            cstat[1][k][iStat-1][i] += vTrans[i];
		}
    }

    // k = 0 => del uvel   k = 1 => del vvel   k = 2 => del wvel   
    double bcl, bcu;
    double dx10, dx20, dx30, dx21, dx31, dx32;
    double dx01, dx02, dx03, dx12, dx13, dx23;
    dx10 = odtl->pos->d.at(0)-odtl->posf->d.at(0); dx01 = odtl->posf->d.at(N)-odtl->pos->d.at(N-1);
    dx20 = odtl->pos->d.at(1)-odtl->posf->d.at(0); dx02 = odtl->posf->d.at(N)-odtl->pos->d.at(N-2);
    dx30 = odtl->pos->d.at(2)-odtl->posf->d.at(0); dx03 = odtl->posf->d.at(N)-odtl->pos->d.at(N-3);
    dx21 = odtl->pos->d.at(1)-odtl->pos->d.at(0);  dx12 = odtl->pos->d.at(N-1)-odtl->pos->d.at(N-2);
    dx31 = odtl->pos->d.at(2)-odtl->pos->d.at(0);  dx13 = odtl->pos->d.at(N-1)-odtl->pos->d.at(N-3);
    dx32 = odtl->pos->d.at(2)-odtl->pos->d.at(1);  dx23 = odtl->pos->d.at(N-2)-odtl->pos->d.at(N-3);
    for(int k=0; k<(int)cstats[0].size(); k++){
        for(int i=0; i<N; i++){
            temp[i] = cstats[2][k][i];
        }
        bcl = ( (dx20*temp[0]-dx10*temp[1]) *dx30 *dx31
                +(dx10*temp[2]-dx30*temp[0]) *dx20 *dx21 )
                / (dx21*dx31*dx32);
        bcu = ( (dx02*temp[N-1]-dx01*temp[N-2]) *dx03 *dx13
                -(dx03*temp[N-1]-dx01*temp[N-3]) *dx02 *dx12 )
                / (dx12 *dx13 *dx23);
        odtGrd2statGrd(odtl->posf->d, temp);
		for(int i=0; i<(int)vTrans.size(); i++){
            cstat[2][k][iStat-1][i] += vTrans[i];
        }
    }

	return;
}
///////////////////////////////////////////////////////////////////////////////
/**
 *  BSetOld saves the current odtline converted to the stats grid to be used in 
 *  BChange to calculate the difference generated through an eddy event or a 
 *  diffusion step. The conversion to the stats grid is needed due to the change 
 *  of the odt grid during the process.
 *
 */
void solver::statsSetOld(){

	odtGrd2statGrd(odtl->posf->d, odtl->uvel->d);
    for (int i=0; i<(int)vTrans.size(); i++){
        oldVars[0][i] = vTrans[i]-uMeanS[i];
        oldVars[1][i] = 0.5*pow(vTrans[i]-uMeanS[i], 2);
    }
    odtGrd2statGrd(odtl->posf->d, odtl->vvel->d);
    for (int i=0; i<(int)vTrans.size(); i++)
        oldVars[2][i] = 0.5*pow(vTrans[i]-vMeanS[i], 2);

    odtGrd2statGrd(odtl->posf->d, odtl->wvel->d);
    for (int i=0; i<(int)vTrans.size(); i++)
        oldVars[3][i] = 0.5*pow(vTrans[i]-wMeanS[i], 2);

	if(odtl->odtp->partCoupl > 0 && odtl->part->particleON) odtl->part->statsSetOld();

    return;
}

///////////////////////////////////////////////////////////////////////////////
/**
 *  statsChange saves and gathers the changes generated through a process (eddy 
 *  event, diffusion, ...) into edstat.
 *
 *  @param jj    \input switch: 0: eddy event  2: diffusion  4: all processes
 */
void solver::statsChange(const int &jj){
    // jj = 0  ==>  differences caused by eddies
    // jj = 1  ==>  differences caused by diffusion
    // jj = 2  ==>  differences caused by eddies and diffusion (all changes) 
    // jj = 3  ==>  differences caused by adaption after diffusion and eddies
    // jj = 4  ==>  differences caused by kernel events

    double dx = odtl->odtp->domainLength/odtl->odtp->ngrdStat; 

    odtGrd2statGrd(odtl->posf->d, odtl->uvel->d);
    for (int i=0; i<(int)vTrans.size(); i++){
        edstat[0+jj][0][iStat-1][i] += posS[i]*(vTrans[i]-uMeanS[i]-oldVars[0][i]); 
        edstat[0+jj][1][iStat-1][i] += 0.5*pow(vTrans[i]-uMeanS[i], 2) - oldVars[1][i];
    }

    odtGrd2statGrd(odtl->posf->d, odtl->vvel->d);
    for (int i=0; i<(int)vTrans.size(); i++){
        edstat[0+jj][2][iStat-1][i] += 0.5*pow(vTrans[i]-vMeanS[i], 2) - oldVars[2][i];
    }

    odtGrd2statGrd(odtl->posf->d, odtl->wvel->d);
    for (int i=0; i<(int)vTrans.size(); i++){
		edstat[0+jj][3][iStat-1][i] += 0.5*pow(vTrans[i]-wMeanS[i], 2) - oldVars[3][i];
    }

    double E = 0.0;
    for(int i=0; i<(int)vTrans.size(); i++){
        E += dx*(edstat[0+jj][0][iStat-1][i]+edstat[0+jj][1][iStat-1][i]+edstat[0+jj][2][iStat-1][i]);
    }

    if(odtl->odtp->partCoupl > 0 && odtl->part->particleON && jj < 4) odtl->part->statsChange(jj);

    return;
}

///////////////////////////////////////////////////////////////////////////////
/**
 *  The function BSnap uses the gatherd data to calculat the mean values, 
 *  the variances, the tke, the budget terms of tke and further data using the 
 *  stats arrays ctime, cstat, edstat, phstat. After calculation the information 
 *  is stored in a file located in the data folder.
 *
 *  @input myid \input the id of the current processor (only used for mpi runs)
 */
void solver::statsOutput(){

    cout << "Start subroutine statsOutput" << endl;

    int    N  = 0;
    double dz = 0.0;
    double fluxfac2 = 0.0;

    // deklaration and initialisation of needed arrays
    vector< vector< vector< vector<double> > > > eavg;
    vector< vector< vector< vector<double> > > > cavg;

    // eavg(j,k,l,i):   i = coordinate
    //                  j = del_0.5*U^2_eddy, del_0.5*U^2_flow, del_0.5*U^2_all,
    //                      del_0.5*U^2_adpt, 0.5*del_U_kern^2
    //                  k = u, u2, v2, w2
    //                  l = average periode
    eavg = vector<vector<vector<vector<double> > > > (5,
           vector<vector<vector<double> > > (4,
           vector<vector<double> > (odtl->odtp->nStat,
           vector<double>(ngrdS ,0.0) ) ) );
    // cavg(m,k,l,i):     i = coordinate
    //                    m = u, u^2, du
    //                    k = u, v, w
    //                    l = average periode
    cavg = vector<vector<vector<vector<double> > > > (3,
           vector<vector<vector<double> > > (3,
           vector<vector<double> > (odtl->odtp->nStat,
           vector<double>(ngrdS ,0.0) ) ) );

    vector< vector<double> > tke;
    vector< vector<double> > d_diff;
    vector< vector<double> > d_eddy;
    vector< vector<double> > d_all;
    vector< vector<double> > d_adpt;
    vector< vector<double> > d_kern;
    vector< vector<double> > bal;
	vector<double> temp;

    // tke(i,l):        i = coordinate
    //                  l = average periode
    // tv, d, ta, p have the same as tke
    tke    = vector<vector<double> > (odtl->odtp->nStat, vector<double>(ngrdS ,0.0) );
	d_diff = vector<vector<double> > (odtl->odtp->nStat, vector<double>(ngrdS ,0.0) );
	d_eddy = vector<vector<double> > (odtl->odtp->nStat, vector<double>(ngrdS ,0.0) ); 
	d_all  = vector<vector<double> > (odtl->odtp->nStat, vector<double>(ngrdS ,0.0) );
	d_adpt = vector<vector<double> > (odtl->odtp->nStat, vector<double>(ngrdS ,0.0) );
	d_kern = vector<vector<double> > (odtl->odtp->nStat, vector<double>(ngrdS ,0.0) );
	bal    = vector<vector<double> > (odtl->odtp->nStat, vector<double>(ngrdS ,0.0) );
	temp = vector<double> (ngrdS, 0.0);
	
    cout << "statsOutput :: calculate cavg" << endl;
    for(int i = 0; i < odtl->odtp->nStat; i++){
        for(int k = 0; k < 3; k++){
            for(int l = 0; l < 3; l++){
                for(int j = 0; j < ngrdS; j++){
                    cavg[k][l][i][j] = cstat[k][l][i][j] / ctime[i];
                }
            }
        }
    }

    cout << "statsOutput :: calculate eavg" << endl;
    for(int k = 0; k < odtl->odtp->nStat; k++){
        for(int i = 0; i < 5; i++){
            for(int l = 0; l < 4; l++){
                for(int j = 0; j < ngrdS; j++){
                    eavg[i][l][k][j] = edstat[i][l][k][j] / ctime[k];
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    // Calculation of the butget terms of the turbulent kinetic energy
    //
    // --> assumptions: statistically stationary problem like channel flow
    //
    // ToDo: The derivations at the boundaries are currently calculated as if 
    //       the grid is a finite differences grid -> 1st point is on the wall.
    //--------------------------------------------------------------------------

    // begin calculation of turbulent kinetic energy
    for (int k = 0; k < odtl->odtp->nStat; k++){
        N        = ngrdS;
        dz       = LdomainS / (N*1.0);

        for (int j = 0; j < N; j++){

            tke[k][j]     = 0.5*cavg[1][0][k][j] ; // <u'^2>
            tke[k][j]    += 0.5*cavg[1][1][k][j] ; // <v'^2>
            tke[k][j]    += 0.5*cavg[1][2][k][j] ; // <w'^2>
            d_eddy[k][j]  = eavg[0][1][k][j] + eavg[0][2][k][j] + eavg[0][3][k][j]; // eddy event
	        d_diff[k][j]  = eavg[1][1][k][j] + eavg[1][2][k][j] + eavg[1][3][k][j]; // diffusion
	        d_all[k][j]   = eavg[2][1][k][j] + eavg[2][2][k][j] + eavg[2][3][k][j]; // eddy + diffusion
	        d_adpt[k][j]  = eavg[3][1][k][j] + eavg[3][2][k][j] + eavg[3][3][k][j]; // mesh adaption
	        d_kern[k][j]  = eavg[4][1][k][j] + eavg[4][2][k][j] + eavg[4][3][k][j]; // kernel event
        }
    }

    for (int k = 0; k < odtl->odtp->nStat; k++){
        for (int j = 0; j < ngrdS; j++){
            bal[k][j] = d_eddy[k][j] + d_diff[k][j] + d_adpt[k][j] + d_kern[k][j];
        }
    }

// start data output

    cout << "statsOutput :: start writing output" << endl;
    // -- writing ASCII output
    double temp2 = 0.0;
    double temp3 = 0.0;
	N = ngrdS;

    for (int k = 0; k < odtl->odtp->nStat; k++){
        stringstream ss1; string s1;
        ss1.clear(); s1.clear();
		ss1 << setfill('0') << setw(5) << k; ss1 >> s1;
		string name = "../data/test/data/data-average-output_" + s1 + ".dat";
		ofstream file_avg(name.c_str(),ios::app);
		assert(file_avg.is_open());
        file_avg << "# grid points = " << N;
        file_avg << "\n# Domain Size = " << LdomainS;
        // could be used for single values to be in the output file
        // the values to be written in the order of the previous line
        file_avg << "\n#          l_min          l_ave"
                 <<    "          l_max            rho"
                 <<    "          molec           visc"
                 <<    "              t";
        file_avg << setprecision(8);
        file_avg << "\n# " << setw(14) << odtl->odtp->Lmin * odtl->odtp->ngrdStat
                 <<    " " << setw(14) << odtl->odtp->Lp   * odtl->odtp->ngrdStat
                 <<    " " << setw(14) << odtl->odtp->Lmax * odtl->odtp->ngrdStat
                 <<    " " << setw(14) << odtl->odtp->rho0
                 <<    " " << setw(14) << odtl->odtp->kvisc0
                 <<    " " << setw(14) << odtl->odtp->kvisc0/odtl->odtp->rho0
                 <<    " " << setw(14) << ctime.at(k);
//                 <<    " " << setw(14) << odtl->odtp->tEnd * (k+1) / odtl->odtp->nStat;
        file_avg <<  "\n#  y                      "
                 <<     "  u                      "
                 <<     "  v                      "
                 <<     "  w                      "
                 <<     "  uP2                     "
                 <<     "  vP2                    "
                 <<     "  wP2                    "
                 <<     "  tke                    " 
                 <<     "  d_eddy                 " 
                 <<     "  d_diff                 "
                 <<     "  d_all                  "
                 <<     "  d_adpt                 "
                 <<     "  d_kern                 "
                 <<     "  bal                    ";
        file_avg << scientific;
        file_avg << setprecision(16);
        temp2 = 0.0;
        temp3 = 0.0;
        for (int j = 0; j < N; j++){
            temp2 += dz * eavg[0][0][k][j];
            file_avg << endl;
            file_avg << setw(25) << posS[j]
                     << setw(25) << cavg[0][0][k][j]
                     << setw(25) << cavg[0][1][k][j]
                     << setw(25) << cavg[0][2][k][j]
                     << setw(25) << cavg[1][0][k][j]
                     << setw(25) << cavg[1][1][k][j]
                     << setw(25) << cavg[1][2][k][j]
                     << setw(25) << tke[k][j]
                     << setw(25) << d_eddy[k][j]
                     << setw(25) << d_diff[k][j]
                     << setw(25) << d_all[k][j]
                     << setw(25) << d_adpt[k][j]
                     << setw(25) << d_kern[k][j]
                     << setw(25) << bal[k][j];
        }
        file_avg.close();
    }
    cout << "statsOutput :: output in ASCII files finished" << endl;
	if(odtl->odtp->partCoupl > 1 && odtl->part->particleON)
		odtl->part->statsOutput();

return;
}



