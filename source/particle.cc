#include "particle.h"
#include "odtline.h"
#include "odtparam.h"
//#include <math.h>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <numeric> //accumulate

/////////////////////////////////////////////////////////////////////
 

void particle::init(odtline *_odtl){ 
        	
    odtl = _odtl;
	initTime = 0;
    particleON = false;

	if(odtl->odtp->partCoupl == 0)
		return;
	//initializing of particle position and velocity from input file
	yPos.resize(odtl->odtp->Nparticle, 0.0);
	linePos.resize(odtl->odtp->Nparticle, 0.0);
	uvelP.resize(odtl->odtp->Nparticle, 0.0);
    vvelP.resize(odtl->odtp->Nparticle, 0.0);
    wvelP.resize(odtl->odtp->Nparticle, 0.0);
    av_uP.resize(odtl->odtp->Nparticle, 0.0);
    av_vP.resize(odtl->odtp->Nparticle, 0.0);
    av_wP.resize(odtl->odtp->Nparticle, 0.0);
	uvelG.resize(odtl->odtp->Nparticle,0);
    vvelG.resize(odtl->odtp->Nparticle,0);
    wvelG.resize(odtl->odtp->Nparticle,0);
	NPartPerParcel.resize(odtl->odtp->Nparticle, odtl->io->initParticleParams[0][0].as<int>());
	pDistr.resize(240,0.0);

	if(odtl->odtp->probType == "COLDJET"){
		double djeti = odtl->io->initParams["djeti"].as<double>();
		for(int i=0; i<odtl->odtp->Nparticle; i++){
			yPos.at(i) = djeti*odtl->rand->getRand()-djeti*0.5;
		}
	}
	else if(odtl->odtp->probType == "SHEARFLOW"){
		double domainLength  = odtl->io->params["domainLength"].as<double>();
    	double shearRate     = odtl->io->initParams["Srate"].as<double>();
		for(int i =0; i< odtl->odtp->Nparticle; i++){
			yPos.at(i) = (odtl->rand->getRand()-0.5)*domainLength;
		}
	}
	else if(odtl->odtp->probType == "CHANNEL"){
    	double domainLength  = odtl->io->params["domainLength"].as<double>();
        for(int i =0; i< odtl->odtp->Nparticle; i++){
        	yPos.at(i) = (odtl->rand->getRand()-0.5)*odtl->odtp->domainLength;
        }
    }
	else if(odtl->odtp->probType == "HIT"){
		double domainLength  = odtl->io->params["domainLength"].as<double>();
		for(int i =0; i< odtl->odtp->Nparticle; i++){
			yPos.at(i) = (odtl->rand->getRand()-0.5)*domainLength;
		}
	}
	else{
		cout<<"\n#Please add initial particle distribution for elected case!";
	}

	for(int i = 0; i < odtl->odtp->Nparticle; i++){
	    linePos.at(i) = yPos.at(i);
    }

    //finding near by cell index for each particle on odtline (upper one)
    indP.resize(odtl->odtp->Nparticle, 0);
	for(int i = 0; i < odtl->odtp->Nparticle; i++){
		indP.at(i) = odtl->linePositionToIndex(linePos.at(i), false, 1);
	}
    // displacement of particles
    displ.resize(odtl->odtp->Nparticle, 0.0);
	//type of particle
    if(odtl->odtp->typeP == 0 && odtl->odtp->partCoupl == 2)
        odtl->odtp->partCoupl = 1;
    typeP.resize(odtl->odtp->Nparticle, odtl->odtp->typeP);
	//activate particle
	activeP.resize(odtl->odtp->Nparticle, true);
	//initializing of particle density (same) 
    dens0P.resize(odtl->odtp->Nparticle, odtl->odtp->DensiParticle);
	//initializing of particle diameter (same)
	diamP.resize(odtl->odtp->Nparticle, odtl->odtp->Dparticle);
	//initializing of particle mass (spherical)
	massP.resize(odtl->odtp->Nparticle, 3.1416*odtl->odtp->Dparticle*odtl->odtp->Dparticle*odtl->odtp->Dparticle/6*odtl->odtp->DensiParticle); 
	//activate particle
	activeP.resize(odtl->odtp->Nparticle, true);
	//initialize gravity vector
	gravity.resize(3,0);
	gravity.at(0) = odtl->odtp->Gx;
	gravity.at(1) = odtl->odtp->Gy;
	gravity.at(2) = odtl->odtp->Gz;
	//initialize Stokes and Reynolds numbers
	f.resize(odtl->odtp->Nparticle, 1);
	tauP.resize(odtl->odtp->Nparticle, 1);

	for(int i =0; i< odtl->odtp->Nparticle; i++){
        uvelP.at(i) = 0.0;
        vvelP.at(i) = 0.0; 
    	wvelP.at(i) = 0.0;
	}

	fDivTauP.resize(odtl->odtp->Nparticle, 1.0);
	uvelPnew.resize(odtl->odtp->Nparticle, 0.0);
    vvelPnew.resize(odtl->odtp->Nparticle, 0.0);
    wvelPnew.resize(odtl->odtp->Nparticle, 0.0);
	uEddy.resize(odtl->odtp->Nparticle, 0.0);
    vEddy.resize(odtl->odtp->Nparticle, 0.0);
    wEddy.resize(odtl->odtp->Nparticle, 0.0);
    yPosnew.resize(odtl->odtp->Nparticle, 0.0);
	zPosnew.resize(odtl->odtp->Nparticle, 0.0);
   	linePosnew.resize(odtl->odtp->Nparticle, 0.0);
 	deltaLinePos.resize(odtl->odtp->Nparticle, 0.0);
	partMomSrc.resize(3, 0.0);
	partEnergSrc.resize(3, 0.0);

	crmax = odtl->odtp->crmax; // Maximum relative velocity between particles
	
	pJumpU = 0.0;
	pJumpL = 0.0;

    if(odtl->odtp->partCoupl > 0 && odtl->odtp->tparticleON == 0.0){
        particleON = true;
        initTime = odtl->mimx->time;
        setVelocity();
    }
    
	writeData("particle_init_00000.dat", 0.0);
}
/////////////////////////////////////////////////////////////////////
/** Set velocity for each particle to local gas velocity  
*
*/
void particle::setVelocity(){
	for(int i =0; i< odtl->odtp->Nparticle; i++){
        getGasVelocity(i, odtl->odtp->cCoord);
        uvelP.at(i) = uvelG.at(i);
        vvelP.at(i) = vvelG.at(i);
        wvelP.at(i) = wvelG.at(i);
    }
}

/////////////////////////////////////////////////////////////////////
/** Particle output 
*@param fname \input Name for output file
*@param time \input Time step of Simulation
*(Very depending on test case)
*/ 
void particle::writeData(string fname, const double time){
	
	if(odtl->odtp->partCoupl == 0)
                return;

	stringstream ss1;
	string s1;

	ss1.clear();  ss1 << setfill('0') << setw(5) << odtl->io->iNextDumpTime; ss1 >> s1;
    string name = "../data/test/particle_data/"+fname;

    ofstream write_output(name.c_str(),ios::app);
    assert(write_output.is_open());
	write_output << "# Particle properties";
	write_output << "\n# Time = "   << time ;
    write_output << "\n# Number of particles = "   << odtl->odtp->Nparticle ;
    write_output << "\n#";

	write_output <<  setw(15) << "0_particleID"<< setw(19) << "1_yPos" << setw(19) << "2_displ" << setw(19) << "3_massP" << setw(19) <<"4_NPartPerParcel" << setw(19) << "5_fDivTauP" << setw(19) << "6_uvelP" << setw(19) << "7_vvelP" << setw(19) << "8_wvelP" << setw(19) << "9_uvelG" << setw(19) << "10_vvelG" << setw(19) << "11_wvelG" << setw(19) << "12_uRhsSrc" << setw(19) << "13_vRhsSrc" << setw(19) << "14_wRhsSrc" << endl;
	for(int i=0; i<odtl->odtp->Nparticle; i++) {
			write_output << scientific;
            write_output << setprecision(10);
            write_output << setw(16) << i;      
            write_output << setw(19) << yPos.at(i);
	        write_output << setw(19) << displ.at(i);
 			write_output << setw(19) << massP.at(i);
			write_output << setw(19) << NPartPerParcel.at(i);
            write_output << setw(19) << fDivTauP.at(i);
            write_output << setw(19) << uvelP.at(i);
            write_output << setw(19) << vvelP.at(i);
            write_output << setw(19) << wvelP.at(i); 
			write_output << setw(19) << uvelG.at(i);
            write_output << setw(19) << vvelG.at(i);
            write_output << setw(19) << wvelG.at(i)<< "\n";
	}
	write_output.close();

}
/////////////////////////////////////////////////////////////////////
/** Advance of particle governed by drag law. 
* Computes also the momentum exchange with gas phase (Two-way coupling). 
* Computes weights for gas phase momentum equation due to explicit Euler scheme. 
* Otherwise the RHS can be higher as the physical limit of two-way coupling.
* @param dxc \input Vector contains cell volumes
* @param dt \input Time step for advancement
*/
    
void particle::advanceParticleAndCalcWeights(const vector<double> &dxc, double dt){
		
	pJumpU = 0.0;
	pJumpL = 0.0;
	
  	if(odtl->odtp->partCoupl == 0)
                return;
	// initialization of particle distribution/collision grid	
	int N = 240;
	for(int i=0; i<N;i++){
		pDistr.at(i)=0.0;
	}
	double dx = odtl->odtp->domainLength/N;
	
    /// local gas velocities at particle position (linear interpolation)
	for(int k=0; k < odtl->odtp->Nparticle; k++){
		getGasVelocity(k, odtl->odtp->cCoord);
	}
    if(odtl->odtp->partCoupl == 2){
        /// calculate weights and average particle velocity for each cell (check limits of mom tranfer)
        sumWeights.resize(dxc.size(), 0.0);  
        av_uP.resize(dxc.size(), 0.0);
        av_vP.resize(dxc.size(), 0.0);
        av_wP.resize(dxc.size(), 0.0); 
        for(int i=0; i<dxc.size();i++){
            sumWeights.at(i)=0.0;
            av_uP.at(i)=0.0;
            av_vP.at(i)=0.0;
            av_wP.at(i)=0.0;
        } 
    
        double mc;
        for(int i=0; i<odtl->odtp->Nparticle; i++){ 
            set_fDivTauP(i, uvelP.at(i), vvelP.at(i), wvelP.at(i), uvelG.at(i), vvelG.at(i), wvelG.at(i));
            mc = (odtl->rho->d.at(indP.at(i))*dxc.at(indP.at(i)))/odtl->odtp->cCoord;
            double w  = massP.at(i)*NPartPerParcel.at(i)*fDivTauP.at(i);
            sumWeights.at(indP.at(i)) += w;
            av_uP.at(indP.at(i)) += w*uvelP.at(i); 
            av_vP.at(indP.at(i)) += w*vvelP.at(i);
            av_wP.at(indP.at(i)) += w*wvelP.at(i);
        }
    
        for(int i=0; i<dxc.size(); i++){
            if(sumWeights.at(i) != 0.0){
                av_uP.at(i) /= sumWeights.at(i);
                av_vP.at(i) /= sumWeights.at(i);
                av_wP.at(i) /= sumWeights.at(i);
                mc = (odtl->rho->d.at(i)*dxc.at(i))/odtl->odtp->cCoord;
                sumWeights.at(i) /= mc;
            }
        }
    }
      
	/// Particle advancement governed by drag law
	/// explicit Euler & seperation tracer(0), inertial(1), ballistic(2)
        vector<double> y0 = yPos;
	for(int i = 0; i<odtl->odtp->Nparticle; i++){ 
	    if(activeP.at(i)){
			double timeStep = dt;
            if(typeP.at(i) == 0){   //tracer particle (tauP -> 0)
                uvelP.at(i) = uvelG.at(i);
                wvelP.at(i) = wvelG.at(i);
                vvelP.at(i) = vvelG.at(i);
                yPos.at(i) = yPos.at(i) + timeStep*vvelP.at(i);
            }
            else if(typeP.at(i) == 1){  // inertial particle
                set_fDivTauP(i, uvelP.at(i), vvelP.at(i), wvelP.at(i), uvelG.at(i), 0.0, wvelG.at(i));
                int fac = ceil(fDivTauP.at(i)*timeStep);
                for(int j=0; j<fac; j++){
                    double u0 = uvelP.at(i);
                    double w0 = wvelP.at(i);
                    uvelP.at(i) = uvelP.at(i) + timeStep/fac*(-(uvelP.at(i)-uvelG.at(i))*fDivTauP.at(i)+ gravity.at(0));
                    wvelP.at(i) = wvelP.at(i) + timeStep/fac*(-(wvelP.at(i)-wvelG.at(i))*fDivTauP.at(i)+ gravity.at(2));
                    vvelP.at(i) = vvelP.at(i) + timeStep/fac*(-(vvelP.at(i))*fDivTauP.at(i)+ gravity.at(1));
                    yPos.at(i) = yPos.at(i) + timeStep/fac*vvelP.at(i);
                    if(yPos.at(i) < -0.5*odtl->odtp->domainLength || yPos.at(i) > 0.5*odtl->odtp->domainLength)
                        cout <<"Stop!"; 
                }
            }
            else if(typeP.at(i) == 2){    // ballistic particles (tauP -> inf)
                uvelP.at(i) = uvelP.at(i) + timeStep*gravity.at(0);
                vvelP.at(i) = vvelP.at(i) + timeStep*gravity.at(1);
                wvelP.at(i) = wvelP.at(i) + timeStep*gravity.at(2);
                yPos.at(i)  = yPos.at(i)  + timeStep*vvelP.at(i);
            }
            else
                cout << "\n Error: No selectable particle type!" << endl;
       
			checkBCandSetLinePos(i);				
			partLinePosToIndex(i);

			// Particle distribution vector
			int iCell = 0;
			if(linePos.at(i) > 0)
				iCell = int(N/2)+int(linePos.at(i)/dx);
			else
				iCell = int(N/2)+int(linePos.at(i)/dx)-1;
			if(iCell > N-1)
				iCell = N-1;
			pDistr.at(iCell) += 1*NPartPerParcel.at(i);
		}
	}
	getDisplacement(y0, yPos);
//------particle collision
	if(odtl->odtp->partCoupl == 3){
		collision(dxc, dt);
	}
}
/////////////////////////////////////////////////////////////////////
/** check boundary conditions and set new line position after diffusion 
* @param i \input Partile ID
 */
void particle::checkBCandSetLinePos(int i){

    //---------------------Set odtline position 
    linePos.at(i) = yPos.at(i);

    //--------------------- Boundary conditions (0 = Outflow, 1 = Periodic, 2 = Wall (Elastic))
    if(odtl->odtp->BCparticle == 0){
        if(linePos.at(i) > odtl->posf->d.at(odtl->ngrdf-2) || linePos.at(i) < odtl->posf->d.at(2) || uvelP.at(i) < 0.3){
            if(typeP.at(i) == 1)
                activeP.at(i) = false;
        }
    }
    else if(odtl->odtp->BCparticle == 1){
        if(linePos.at(i) > odtl->posf->d.at(odtl->ngrdf-1)){
			pJumpL = odtl->posf->d.at(odtl->ngrdf-1)- odtl->posf->d.at(0);
            pJumpU = odtl->uvel->d.at(odtl->ngrd-1) - odtl->uvel->d.at(0);
            //pJumpU = odtl->odtp->probType == "SHEARFLOW" ? odtl->io->initParams["Srate"].as<double>()*odtl->io->params["domainLength"].as<double>() : 0.0;
			yPos.at(i) -= pJumpL;
			uvelP.at(i) -= uvelG.at(i) + pJumpU;
        }
        else if(linePos.at(i) < odtl->posf->d.at(0)){
			pJumpL = odtl->posf->d.at(odtl->ngrdf-1)- odtl->posf->d.at(0);
            pJumpU = odtl->uvel->d.at(odtl->ngrd-1) - odtl->uvel->d.at(0);
			yPos.at(i) += pJumpL;
			uvelP.at(i) += pJumpU;
        }
    }
    else if(odtl->odtp->BCparticle == 2){
        if(linePos.at(i) > odtl->posf->d.at(odtl->ngrdf-1)){
            yPos.at(i) = odtl->pos->d.at(odtl->ngrd-1);
            vvelP.at(i) = -1.0*vvelP.at(i);
        }
        if(linePos.at(i) < odtl->posf->d.at(1)){
            yPos.at(i) = odtl->pos->d.at(1);
            vvelP.at(i) = -1.0*vvelP.at(i);
        }
    }
    else{
        cout << endl << "Unknown boundary condition for particular phase.";
        exit(0);
    }


  //---------------------Update odtline position

    linePos.at(i) = yPos.at(i);
}
/////////////////////////////////////////////////////////////////////
/** Find cell index for particle line position 
* @param i \input Particle ID
 */
void particle::partLinePosToIndex(int i){

	if(activeP.at(i))
    	indP.at(i) = odtl->linePositionToIndex(linePos.at(i), false, 2);
}
/////////////////////////////////////////////////////////////////////
/** Get particle displacement 
* @param y0 \input initial particle position
* @param yNew \input new particle position
 */
void particle::getDisplacement(vector<double> y0, vector<double> yNew){
	
    double dy;	
    for(int i = 0; i<odtl->odtp->Nparticle; i++){
        if(activeP.at(i)){
	        dy = yNew.at(i) -y0.at(i);
            if(abs(dy) > 0.5*odtl->odtp->domainLength)
                dy = dy < 0 ? odtl->odtp->domainLength + dy : dy-odtl->odtp->domainLength;
	        displ.at(i) += dy;
	    }
    }
}
/////////////////////////////////////////////////////////////////////
/** Find corresponding gas phase velocity for particle line position 
* @param i \input Particle ID
* @param cCoord \input coordinate system (1 = planar, 2 = cylindrical, 3 = spherical)
 */
void particle::getGasVelocity(int k, int cCoord){
	    //find index of gas phase cell 
		partLinePosToIndex(k);
		if(activeP.at(k)){
			int l = indP.at(k);
			double posFraction;
			double uG;
            double vG; 
            double wG; 
            // if particle pos smaller than cell center 
            if(linePos.at(k) < odtl->pos->d.at(l)){
                if(l == 0){ 
                    //--------------------- Boundary conditions (0 = Outflow, 1 = Periodic, 2 = Wall (Elastic))
                    if(odtl->odtp->BCparticle == 0){
                        posFraction = (odtl->pos->d.at(l)-linePos.at(k))/(odtl->pos->d.at(l+1)-odtl->pos->d.at(l));
                        uG = odtl->uvel->d.at(l) + (odtl->uvel->d.at(l+1)-odtl->uvel->d.at(l))*posFraction;
                        vG = odtl->vvel->d.at(l) + (odtl->vvel->d.at(l+1)-odtl->vvel->d.at(l))*posFraction;
                        wG = odtl->wvel->d.at(l) + (odtl->wvel->d.at(l+1)-odtl->wvel->d.at(l))*posFraction;
                    }
                    else if(odtl->odtp->BCparticle == 1){
                        pJumpU =  odtl->odtp->probType == "SHEARFLOW" ? odtl->io->initParams["Srate"].as<double>()*odtl->io->params["domainLength"].as<double>() : 0.0;
                        posFraction = (odtl->odtp->domainLength + linePos.at(k)-odtl->pos->d.at(odtl->ngrd-1))/(odtl->odtp->domainLength + odtl->pos->d.at(l)-odtl->pos->d.at(odtl->ngrd-1)); 
                        uG = odtl->uvel->d.at(odtl->ngrd-1)-pJumpU + (odtl->uvel->d.at(l)-odtl->uvel->d.at(odtl->ngrd-1)+pJumpU)*posFraction;
                        vG = odtl->vvel->d.at(odtl->ngrd-1) + (odtl->vvel->d.at(l)-odtl->vvel->d.at(odtl->ngrd-1))*posFraction;
                        wG = odtl->wvel->d.at(odtl->ngrd-1) + (odtl->wvel->d.at(l)-odtl->wvel->d.at(odtl->ngrd-1))*posFraction;
                    }
                    else{
                        posFraction = (linePos.at(k)-odtl->posf->d.at(l))/(odtl->pos->d.at(l)-odtl->posf->d.at(l));
                        uG = 0.0 + (odtl->uvel->d.at(l)-0.0)*posFraction;
                        vG = 0.0 + (odtl->vvel->d.at(l)-0.0)*posFraction;
                        wG = 0.0 + (odtl->wvel->d.at(l)-0.0)*posFraction;
                    }
                }
                else{
                    posFraction = (linePos.at(k)-odtl->pos->d.at(l-1))/(odtl->pos->d.at(l)-odtl->pos->d.at(l-1));
                    uG = odtl->uvel->d.at(l-1) + (odtl->uvel->d.at(l)-odtl->uvel->d.at(l-1))*posFraction;
                    vG = odtl->vvel->d.at(l-1) + (odtl->vvel->d.at(l)-odtl->vvel->d.at(l-1))*posFraction;
                    wG = odtl->wvel->d.at(l-1) + (odtl->wvel->d.at(l)-odtl->wvel->d.at(l-1))*posFraction;
                }
            }
            else{
                if(l == odtl->ngrd-1){
                    if(odtl->odtp->BCparticle == 0){
                        posFraction = (linePos.at(k)-odtl->pos->d.at(l-1))/(odtl->pos->d.at(l)-odtl->pos->d.at(l-1)); 
                        uG = odtl->uvel->d.at(l-1) + (odtl->uvel->d.at(l)-odtl->uvel->d.at(l-1))*posFraction;
                        vG = odtl->vvel->d.at(l-1) + (odtl->vvel->d.at(l)-odtl->vvel->d.at(l-1))*posFraction;
                        wG = odtl->wvel->d.at(l-1) + (odtl->wvel->d.at(l)-odtl->wvel->d.at(l-1))*posFraction; 
                    }
                    else if(odtl->odtp->BCparticle == 1){
                        pJumpU =  odtl->odtp->probType == "SHEARFLOW" ? odtl->io->initParams["Srate"].as<double>()*odtl->io->params["domainLength"].as<double>() : 0.0;
                        posFraction = (linePos.at(k)-odtl->pos->d.at(l))/(odtl->odtp->domainLength + odtl->pos->d.at(0)-odtl->pos->d.at(l)); 
                        uG = odtl->uvel->d.at(l) + (odtl->uvel->d.at(0)+pJumpU-odtl->uvel->d.at(l))*posFraction;
                        vG = odtl->vvel->d.at(l) + (odtl->vvel->d.at(0)-odtl->vvel->d.at(l))*posFraction;
                        wG = odtl->wvel->d.at(l) + (odtl->wvel->d.at(0)-odtl->wvel->d.at(l))*posFraction;
                    }
                    else{
                        posFraction = (linePos.at(k)-odtl->pos->d.at(l))/(odtl->posf->d.at(l)-odtl->pos->d.at(l));                         
                        uG = odtl->uvel->d.at(l) + (0.0-odtl->uvel->d.at(l))*posFraction;
                        vG = odtl->vvel->d.at(l) + (0.0-odtl->vvel->d.at(l))*posFraction;
                        wG = odtl->wvel->d.at(l) + (0.0-odtl->wvel->d.at(l))*posFraction;
                    }
                }
                else{
                    posFraction = (linePos.at(k)-odtl->pos->d.at(l))/(odtl->pos->d.at(l+1)-odtl->pos->d.at(l));
                    uG = odtl->uvel->d.at(l) + (odtl->uvel->d.at(l+1)-odtl->uvel->d.at(l))*posFraction;
                    vG = odtl->vvel->d.at(l) + (odtl->vvel->d.at(l+1)-odtl->vvel->d.at(l))*posFraction;
                    wG = odtl->wvel->d.at(l) + (odtl->wvel->d.at(l+1)-odtl->wvel->d.at(l))*posFraction;
                }
            } 

			if(typeP.at(k) == 1 && odtl->odtp->typeI)	// inertial particle and type I
				vG = 0;
        	uvelG.at(k) =  uG;
        	vvelG.at(k) =  vG;
        	wvelG.at(k) =  wG;
		}
}
////////////////////////////////////////////////////////////////////
/** Particle-eddy interaction for one-way coupling 
*   Step 1: Computes eddy velocities for drag law 
*	Step 2: Computes particle-eddy interaction time step and new particle properties
*	Step 3: Applies new properties of particles
*	Step 4: Deactivates particles which are out of domain
*	@param time \input current time of ODT line
*/
void particle::oneWayCoupl(double time){
    
	if(odtl->odtp->partCoupl == 1 && odtl->odtp->typeP < 2){
    	if(computeEddyVelocities()){ 
        	timePEIandNewPartProp(odtl->ed->invTauEddy, time);
			applyPartProp();
       	}
	}
}
////////////////////////////////////////////////////////////////////
/** Particle-eddy interaction for two-way coupling
*	Step 1: Computes eddy velocities for drag law
*   Step 2: Iteration till convergence of inverse of eddy time 
*		Step 2.1: Computes particle-eddy interaction time step, new particle properties and momentum change
*		Step 2.2: Recomputes inverse of eddy time considering momentum transfer particle-gas phase
* @param time \input current time of ODT line
* @param iStart \input index of start cell of eddy line
* @param iEnd \input index of end cell of eddy line
* @param C \input coordinate system 1 - planar, 2 - cylindrical, 3 - spherical
* @return true - if eddy fullfills requirements 	  
*/
bool particle::twoWayCoupl(double time,int iStart, int iEnd) {

    if(odtl->odtp->typeI){            // only for type I
        if(computeEddyVelocities()){		// Step 1
            double invTau = odtl->ed->invTauEddy;
            double invTau0 = 1000;
            int it = 0;
            while(abs(invTau0 - invTau)/invTau0 > 1E-3 && it < 1000){ //Step 2
                it += 1;
                invTau0 = invTau;
                if(!timePEIandNewPartProp(invTau0, time)){ //Step 2.1
					for (int i = 0; i < odtl->odtp->Nparticle; i++){
                        uvelPnew.at(i) = 0.0;
                        vvelPnew.at(i) = 0.0;
                        wvelPnew.at(i) = 0.0;
                        yPosnew.at(i) = 0.0;
                        zPosnew.at(i) = 0.0;
                    }
					return false; 
				}
				if(typeP.at(0) != 1) break; // only for inertial particles
                if(!eddyTauPartSrc(iStart,iEnd,odtl->odtp->Z_param)){ // Step 2.2 
                    for (int i = 0; i < odtl->odtp->Nparticle; i++){
                        uvelPnew.at(i) = 0.0;
                        vvelPnew.at(i) = 0.0;
                        wvelPnew.at(i) = 0.0;
                        yPosnew.at(i) = 0.0;
						zPosnew.at(i) = 0.0;
                    }
                    return false;
                }
                invTau = odtl->ed->invTauEddy;
            }
            if(it == 1000){
            //    cout << "\n#---------- Two-way coupling not converged in 1000 iterations! Difference: " <<  invTau0 - invTau ;
				return false;
			}
        }
    }
return true;
}

/////////////////////////////////////////////////////////////////////
/** Computes eddy velocities for drag law
*	Step 1: Computes radial displacement of fluid particle (tracer) due to triplet map
*	Step 2: Set x,z-component of eddy velocity equal gas velocity
*	Step 3: Save radial displacement of fluid particle for future use
*	@return if particle is in eddy region
*/
bool particle::computeEddyVelocities(){
	
    int nPartInEddy = 0;
	periodicEddy = false;
	
	if(odtl->ed->rightEdge > odtl->posf->d.at(odtl->ngrdf-1))
		periodicEddy = true;

	for (int i = 0; i < odtl->odtp->Nparticle; i++){
			
		uEddy.at(i) = 0.0;
		vEddy.at(i) = 0.0;
		wEddy.at(i) = 0.0;
		linePosnew.at(i) = 0.0;
		deltaLinePos.at(i) = 0.0;
		pJumpU = 0.0;
		pJumpL = 0.0;

		if(periodicEddy && linePos.at(i) < odtl->ed->rightEdge - odtl->odtp->domainLength){
			pJumpL = odtl->odtp->domainLength;
            pJumpU = odtl->uvel->d.at(odtl->ngrd-1)-odtl->uvel->d.at(0);
			//pJumpU = odtl->odtp->probType == "SHEARFLOW" ? odtl->io->initParams["Srate"].as<double>()*odtl->io->params["domainLength"].as<double>() : 0.0;
		}
		if(activeP.at(i) && typeP.at(i) != 2){
		//Type I (already in the eddy box, not entering during te)
			if(linePos.at(i)+pJumpL >= odtl->ed->leftEdge && linePos.at(i)+pJumpL <= odtl->ed->rightEdge){
				
				nPartInEddy++;	

				// triplet map for every coordinates (as in eddy.cc/tripMap)
				
				int C = odtl->odtp->cCoord;
				double invC = 1.0/C;
				int    pm1;
				double dmb;
				double fracVleft;
    			double fracVmidl;
    			double fracVrght;
				//double oneThird = 1.0/3.0;

				//--------- Set fracVlft, fracVmid, fracVrgt
				double r1 = (2*odtl->ed->leftEdge+odtl->ed->rightEdge)/3;
    			double r2 = (odtl->ed->leftEdge+2*odtl->ed->rightEdge)/3;

				pm1 =odtl->ed->leftEdge*r1 < 0 ? -1.0 : 1.0;             // handles cells that split x=0
   	 			double Vleft = abs(pow(abs(odtl->ed->leftEdge), C) - pm1*pow(abs(r1), C));

    			pm1 =r1*r2 < 0 ? -1.0 : 1.0;                   // handles cells that split x=0
    			double Vmidl = abs(pow(abs(r1), C) - pm1*pow(abs(r2), C));

    			pm1 =r2*odtl->ed->rightEdge < 0 ? -1.0 : 1.0;            // handles cells that split x=0
    			double Vrght = abs(pow(abs(r2), C) - pm1*pow(abs(odtl->ed->rightEdge), C));

    			double volTot = Vleft + Vmidl + Vrght;
    			fracVleft = Vleft / volTot;
    			fracVmidl = Vmidl / volTot;
    			fracVrght = Vrght / volTot;

				// calculating volume between leftEdge and particle
				pm1 = (linePos.at(i)+pJumpL)*odtl->ed->leftEdge < 0 ? -1.0 : 1.0;
				double volToPart = abs(pow(abs(linePos.at(i)+pJumpL),C) - pm1*pow(abs(odtl->ed->leftEdge),C));

				double rand = odtl->rand->getRand();
				int x;

				if(rand < 0.3333)
                    x = 1;
                else if(rand > 0.66666)
                	x = 2;
                else                    
                	x = 3;

				// choosing place for particle
				if(x==1){
        			pm1 = odtl->ed->leftEdge < 0 ? -1 : 1;
					dmb = pow(abs(odtl->ed->leftEdge),C) + pm1*fracVleft*volToPart;
        			linePosnew.at(i) = dmb < 0.0 ? pow(-dmb, invC) : pm1*pow(dmb, invC);
    			}
				else if(x==2){
					pm1 = r1 < 0 ? -1 : 1;
                    dmb = pow(abs(r1),C) + pm1*fracVmidl*(volTot-volToPart);
                    linePosnew.at(i) = dmb < 0.0 ? pow(-dmb, invC) : pm1*pow(dmb, invC);
     			}   
				else{ 		//x==3
					pm1 = r2 < 0 ? -1 : 1;
					dmb = pow(abs(r2),C) + pm1*fracVrght*(volToPart);
                    linePosnew.at(i) = dmb < 0.0 ? pow(-dmb, invC) : pm1*pow(dmb, invC);
				}

                if(typeP.at(i) == 0){
                    yPos.at(i) = linePosnew.at(i);
                    return false;
                }

          		/// eddy velocities
				getGasVelocity(i,1);
				uEddy.at(i) = uvelG.at(i)+pJumpU;
				wEddy.at(i) = wvelG.at(i);
				deltaLinePos.at(i) = linePosnew.at(i)-(linePos.at(i)+pJumpL);
			}
		}
	}
	if(nPartInEddy == 0)
        return false;
	else	
		return true;
}
/////////////////////////////////////////////////////////////////////
/** Simulates particle-eddy interaction time and new particle properties after PEI
*	PEI time: An eddy is assumed to be represented by an eddy box [L x L x L] (L - eddy length).
*			  PEI time is the time a particle exits the eddy box calculated with an analytical solution of the drag law.
*	       	  (Critical point in the following means extrema of particle trajectory)\
*	PEI model is only governing lateral displacement and velocity change during eddy event!!! 
* 	Step 1: Calculate lateral eddy velocity (displacement divided by eddy life
*	Step 2: Compute eddy box exit time of particle in all three directions 
*	Step 3: Minimum of exit times is PEI time
*	Step 4: Compute new lateral position and velocity
*	Step 5: Compute momentum change during eddy event	
*	@param invTauEddy \input inverse of eddy time
* 	@param time	\input current time od ODT line
*/
bool particle::timePEIandNewPartProp(double invTauEddy, double time){

	for (int i = 0; i < odtl->odtp->Nparticle; i++){
		vvelPnew.at(i) = 0.0;
        uvelPnew.at(i) = 0.0;
        wvelPnew.at(i) = 0.0;
		yPosnew.at(i) = 0.0;
		zPosnew.at(i) = 0.0;
	}

	partMomSrc.at(0) = 0.0;
	partMomSrc.at(1) = 0.0;
    partMomSrc.at(2) = 0.0;
    partEnergSrc.at(0) = 0.0;
	partEnergSrc.at(1) = 0.0;
	partEnergSrc.at(2) = 0.0;
	
	double eddySize = odtl->ed->eddySize;
	double leftEdge = odtl->ed->leftEdge;
	double rightEdge = odtl->ed->rightEdge;

	double Txmin;         // interaction time in x, y, z, directions.
    double Tymin;
    double Tzmin;
	double eddyScale = 1.0/invTauEddy; // [s] (temporal) 

	for(int i = 0; i < odtl->odtp->Nparticle; i++){
		double eddyLife = eddyScale; 
        vEddy.at(i)= deltaLinePos.at(i)/(odtl->odtp->betaP*eddyLife); //Step 1
    }

    for (int i = 0; i < odtl->odtp->Nparticle; i++){
		double eddyLife = 0;
	
		pJumpL = 0.0;
		pJumpU = 0.0;

        if(periodicEddy && linePos.at(i) < odtl->ed->rightEdge - odtl->odtp->domainLength){
            pJumpL = odtl->odtp->domainLength;
            pJumpU = odtl->uvel->d.at(odtl->ngrd-1)-odtl->uvel->d.at(0);
            //pJumpU = odtl->odtp->probType == "SHEARFLOW" ? odtl->io->initParams["Srate"].as<double>()*odtl->io->params["domainLength"].as<double>() : 0.0;
        }

		eddyLife = odtl->odtp->betaP*eddyScale; 

		if(typeP.at(i) == 0){  //tracer particle
          	vvelPnew.at(i) = vEddy.at(i);
            uvelPnew.at(i) = uEddy.at(i);
            wvelPnew.at(i) = wEddy.at(i);
       	}
        else{		
			if(linePos.at(i)+pJumpL >= odtl->ed->leftEdge && linePos.at(i)+pJumpL <= odtl->ed->rightEdge && activeP.at(i)){
				set_fDivTauP(i, uvelP.at(i)+pJumpU, vvelP.at(i), wvelP.at(i), uEddy.at(i), vEddy.at(i), wEddy.at(i));
				Txmin = 0.0;         // interaction time in x, y, z, directions. (Step 2)
    			Tymin = 0.0;
    			Tzmin = 0.0;
			
				double Tcrit = -1;
            	double dmb;
            	double pPos_relEddyCenter;
            	//=====================================================================================================
                //-  Compute Txmin = minimum of (eddy life) and (time to cross x-eddy box boundary) (z is similar, not y)
                //=====================================================================================================
				Tcrit = -1;
				dmb   = (uEddy.at(i) + gravity.at(0)/fDivTauP.at(i) - (uvelP.at(i)+pJumpU))/(uEddy.at(i) + gravity.at(0)/fDivTauP.at(i)); 
				if(dmb > 0.0)
                	Tcrit = log(dmb)/fDivTauP.at(i);
                	//-------------------
				if(Tcrit < 0.0 || Tcrit >= eddyLife) {       // no critical point	
					pPos_relEddyCenter = pPosRelToEdge(uEddy.at(i), uvelG.at(i), uvelP.at(i)+pJumpU, fDivTauP.at(i), gravity.at(0), 0.0, 0.0, eddyLife);
					if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) )     // in eddy box
                   		Txmin = eddyLife;
                	else if (pPos_relEddyCenter < -eddySize*0.5)   {                                         // outside on left
                   		Txmin = eddyCrossingTime(uEddy.at(i), uvelG.at(i), uvelP.at(i)+pJumpU, fDivTauP.at(i), gravity.at(0), 0.0, -0.5*eddySize, 0.0, eddyLife);
                	}
                	else {                                                                                  // outside on right
                   		Txmin = eddyCrossingTime(uEddy.at(i), uvelG.at(i), uvelP.at(i)+pJumpU, fDivTauP.at(i), gravity.at(0), 0.0,  0.5*eddySize, 0.0, eddyLife);
                	}
          		}
				//------------------
				else {                  // has a critical point
                	pPos_relEddyCenter = pPosRelToEdge(uEddy.at(i), uvelG.at(i), uvelP.at(i)+pJumpU, fDivTauP.at(i), gravity.at(0), 0.0, 0.0, Tcrit);
                	if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) ) {   // in eddy box
                		pPos_relEddyCenter = pPosRelToEdge(uEddy.at(i), uvelG.at(i), uvelP.at(i)+pJumpU, fDivTauP.at(i), gravity.at(0), 0.0, 0.0, eddyLife);
                		if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) ) // in eddy box
                			Txmin = eddyLife;
                		else if (pPos_relEddyCenter < -eddySize*0.5)                                       // outside on left
                			Txmin = eddyCrossingTime(uEddy.at(i), uvelG.at(i), uvelP.at(i)+pJumpU, fDivTauP.at(i), gravity.at(0), 0.0, -0.5*eddySize, 0.0, eddyLife);
                		else                                                                               // outside on right
                   			Txmin = eddyCrossingTime(uEddy.at(i), uvelG.at(i), uvelP.at(i)+pJumpU, fDivTauP.at(i), gravity.at(0), 0.0,  0.5*eddySize, 0.0, eddyLife);
                	}
                	else if (pPos_relEddyCenter < -eddySize*0.5)  {                                         // outside on left
                  		Txmin = eddyCrossingTime(uEddy.at(i), uvelG.at(i), uvelP.at(i)+pJumpU, fDivTauP.at(i), gravity.at(0), 0.0, -0.5*eddySize, 0.0, Tcrit);
                	}
                	else {                                                                                 // outside on right
                   		Txmin = eddyCrossingTime(uEddy.at(i), uvelG.at(i), uvelP.at(i)+pJumpU,fDivTauP.at(i), gravity.at(0), 0.0,  0.5*eddySize, 0.0, Tcrit);
                	}
            	}
				
				//=====================================================================================================
				//-------------------  Compute Tymin = minimum of (eddy life) and (time to cross y-eddy box boundary)
				//=====================================================================================================
				Tcrit = -1;    // init to neg, then if a critical point exists, resets to positive; use +/- as test
				dmb   = (vEddy.at(i) + gravity.at(1)/fDivTauP.at(i) - vvelP.at(i))/(vEddy.at(i) + gravity.at(1)/fDivTauP.at(i));
		        if(dmb > 0.0)
                	Tcrit = log(dmb)/fDivTauP.at(i); // derivative of trajectory equals zero (turning point of trajectory)
				//-------------------
				if(Tcrit < 0.0 || Tcrit >= eddyLife) {
					pPos_relEddyCenter = pPosRelToEdge(vEddy.at(i), vvelG.at(i), vvelP.at(i), fDivTauP.at(i), gravity.at(1), linePos.at(i)+pJumpL, 0.5*(leftEdge + rightEdge), eddyLife);
                	if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) )     // in eddy box
                		Tymin = eddyLife;
                	else if (pPos_relEddyCenter < -eddySize*0.5)   {                                         // outside on left
                		Tymin = eddyCrossingTime(vEddy.at(i), vvelG.at(i), vvelP.at(i), fDivTauP.at(i), gravity.at(1), linePos.at(i)+pJumpL, leftEdge, 0.0, eddyLife);
                	}
                	else {                                                                                  // outside on right
                		Tymin = eddyCrossingTime(vEddy.at(i), vvelG.at(i), vvelP.at(i), fDivTauP.at(i), gravity.at(1), linePos.at(i)+pJumpL, rightEdge, 0.0, eddyLife);
                	}
            	}
            	//-------------------
				else {                  // has a critical point
					pPos_relEddyCenter = pPosRelToEdge(vEddy.at(i), vvelG.at(i), vvelP.at(i), fDivTauP.at(i), gravity.at(1), linePos.at(i)+pJumpL, 0.5*(leftEdge + rightEdge), Tcrit);
                	if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) ) {   // in eddy box
				    	pPos_relEddyCenter = pPosRelToEdge(vEddy.at(i), vvelG.at(i), vvelP.at(i), fDivTauP.at(i), gravity.at(1), linePos.at(i)+pJumpL, 0.5*(leftEdge + rightEdge), eddyLife);
                    	if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) ) // in eddy box
                     		Tymin = eddyLife;
                    	else if (pPos_relEddyCenter < -eddySize*0.5) // outside on left
                       		Tymin = eddyCrossingTime(vEddy.at(i), vvelG.at(i), vvelP.at(i), fDivTauP.at(i), gravity.at(1), linePos.at(i)+pJumpL, leftEdge, 0.0, eddyLife);
                    	else                                         // outside on right
                       		Tymin = eddyCrossingTime(vEddy.at(i), vvelG.at(i), vvelP.at(i), fDivTauP.at(i), gravity.at(1), linePos.at(i)+pJumpL,  rightEdge, 0.0, eddyLife);
                	}
                	else if (pPos_relEddyCenter < -eddySize*0.5)  {                                         // outside on left
						Tymin = eddyCrossingTime(vEddy.at(i), vvelG.at(i), vvelP.at(i), fDivTauP.at(i), gravity.at(1), linePos.at(i)+pJumpL, leftEdge, 0.0, Tcrit);
                	}
                	else {                                                                          // outside on right
						Tymin = eddyCrossingTime(vEddy.at(i), vvelG.at(i), vvelP.at(i),fDivTauP.at(i), gravity.at(1), linePos.at(i)+pJumpL, rightEdge, 0.0, Tcrit);
               		}
           		}
				//=====================================================================================================
				//-  Compute Tzmin = minimum of (eddy life) and (time to cross z-eddy box boundary) (z is similar, not y)
				//=====================================================================================================
				Tcrit = -1;
				dmb   = (wEddy.at(i) + gravity.at(2)/fDivTauP.at(i) - wvelP.at(i))/(wEddy.at(i) + gravity.at(2)/fDivTauP.at(i)); 
	            if(dmb > 0.0)
                	Tcrit = log(dmb)/fDivTauP.at(i);
                //-------------------			
				if(Tcrit < 0.0 || Tcrit >= eddyLife) {       // no critical point
					pPos_relEddyCenter = pPosRelToEdge(wEddy.at(i), wvelG.at(i), wvelP.at(i), fDivTauP.at(i), gravity.at(2), 0.0, 0.0, eddyLife);
                    if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) )     // in eddy box
                    	Tzmin = eddyLife;
                    else if (pPos_relEddyCenter < -eddySize*0.5)   {                                         // outside on left
                    	Tzmin = eddyCrossingTime(wEddy.at(i), wvelG.at(i), wvelP.at(i), fDivTauP.at(i), gravity.at(2), 0.0, -0.5*eddySize, 0.0, eddyLife);
                   	}
                    else {                                                                                  // outside on right
                    	Tzmin = eddyCrossingTime(wEddy.at(i), wvelG.at(i), wvelP.at(i), fDivTauP.at(i), gravity.at(2), 0.0,  0.5*eddySize, 0.0, eddyLife);
                    }
            	}
				//-------------------
				else {                  // has a critical point
                	pPos_relEddyCenter = pPosRelToEdge(wEddy.at(i), wvelG.at(i), wvelP.at(i), fDivTauP.at(i), gravity.at(2), 0.0, 0.0, Tcrit);
                    if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) ) {   // in eddy box
                    	pPos_relEddyCenter = pPosRelToEdge(wEddy.at(i), wvelG.at(i), wvelP.at(i), fDivTauP.at(i), gravity.at(2), 0.0, 0.0, eddyLife);
                    	if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) ) // in eddy box
                        	Tzmin = eddyLife;
                        else if (pPos_relEddyCenter < -eddySize*0.5)                                       // outside on left
                        	Tzmin = eddyCrossingTime(wEddy.at(i), wvelG.at(i), wvelP.at(i), fDivTauP.at(i), gravity.at(2), 0.0, -0.5*eddySize, 0.0, eddyLife);
                        else                                                                               // outside on right
                        	Tzmin = eddyCrossingTime(wEddy.at(i), wvelG.at(i), wvelP.at(i), fDivTauP.at(i), gravity.at(2), 0.0,  0.5*eddySize, 0.0, eddyLife);
              		}
                    else if (pPos_relEddyCenter < -eddySize*0.5)  {                                         // outside on left
                    	Tzmin = eddyCrossingTime(wEddy.at(i), wvelG.at(i), wvelP.at(i), fDivTauP.at(i), gravity.at(2), 0.0, -0.5*eddySize, 0.0, Tcrit);
                  	}
                    else {                                                                                 // outside on right
                    	Tzmin = eddyCrossingTime(wEddy.at(i), wvelG.at(i), wvelP.at(i),fDivTauP.at(i), gravity.at(2), 0.0,  0.5*eddySize, 0.0, Tcrit);
                    }
            	}
				
				double theta;    // minimal eddy interaction time (Step 3)
				//theta = min(eddyLife,min(min(Txmin,Tymin),min(Tymin,Tzmin))); 
				theta = min(eddyLife,Tymin); 
                        
				/// Suns code line 2126-2136 (Step 4)
				uvelPnew.at(i) = uvelP.at(i);
                vvelPnew.at(i) = vvelP.at(i)+vEddy.at(i)*(1-exp(-theta*fDivTauP.at(i)));
                wvelPnew.at(i) = wvelP.at(i);
                yPosnew.at(i) = yPos.at(i) + pJumpL + vEddy.at(i)*theta - vEddy.at(i)/fDivTauP.at(i)*(1-exp(-theta*fDivTauP.at(i)));
				/// MF       (Step 5)
				if(odtl->odtp->partCoupl > 1){
					partMomSrc.at(0) += massP.at(i)*NPartPerParcel.at(i)*(uvelPnew.at(i) - uvelP.at(i)); //should be always 0.0
					partMomSrc.at(1) += massP.at(i)*NPartPerParcel.at(i)*(vvelPnew.at(i) - vvelP.at(i));
					partMomSrc.at(2) += massP.at(i)*NPartPerParcel.at(i)*(wvelPnew.at(i) - wvelP.at(i)); //should be always 0.0
					partEnergSrc.at(0) += 0.5*massP.at(i)*NPartPerParcel.at(i)*(uvelPnew.at(i)*uvelPnew.at(i) - uvelP.at(i)*uvelP.at(i)); //should be always 0.0
				    partEnergSrc.at(1) += 0.5*massP.at(i)*NPartPerParcel.at(i)*(vvelPnew.at(i)*vvelPnew.at(i) - vvelP.at(i)*vvelP.at(i));
        			partEnergSrc.at(2) += 0.5*massP.at(i)*NPartPerParcel.at(i)*(wvelPnew.at(i)*wvelPnew.at(i) - wvelP.at(i)*wvelP.at(i));			        //should be always 0.0
				}
			}
		}
	}
return true;
}
/////////////////////////////////////////////////////////////////////
/*  Compute eddyTau with regard to the particle source terms
* Step 1: Compute kernel integrals
* Step 2: Consider energy and momentum exchange between particle-gas phase
* Step 3: Compute max available energy for each component
* Step 4: Compute new inverse of eddy tau
* @param iStart \input index of start cell of eddy line
* @param iEnd \input index of end cell of eddy line
* @param Z_value \input viscous coefficent
* @return true if eddy has enough energy to exist
*/
bool particle::eddyTauPartSrc(int iStart, int iEnd, const double Z_value) {

	vector<double> dxc(odtl->eddl->ngrd);

    double         KK=0;                                  // equivalent of the 4/27 fac
	double         rhoKK=0;
    double         rhoJK=0;
	double         rhoK=0;
	double         rhoJ=0;
    vector<double> uRhoJ(3,0.0), uRhoK(3,0.0);

    double          eKinEddy     = 0.0;
    double   		eViscPenalty = 0.0;
    double         	ePartEddy    = 0.0;

    vector<double> 	intRhoKi(odtl->eddl->ngrd);
    vector<double> 	intRhoJi(odtl->eddl->ngrd);
	double 			intRho = 0;
    int            	i;

	double rhoEddy=0;
    double viscEddy=0;
    double eLength2 = odtl->ed->eddySize*odtl->ed->eddySize;

   	//////////// cell sizes
	
	dxc = odtl->ed->dxc;
	double VolE = accumulate(dxc.begin(), dxc.end(), 0.0);
 
	///////////// Fill in integral quantities

	//---------- rhoK, rhoJ, UrhoK, UrhoJ
	
	for(i=0; i<odtl->eddl->ngrd; i++) {
		intRhoKi.at(i)  = odtl->ed->K.at(i)*odtl->eddl->rho->d.at(i)*dxc.at(i);                  
       	intRhoJi.at(i)  = abs(intRhoKi.at(i));
		KK       += odtl->ed->K.at(i)*odtl->ed->K.at(i)*dxc.at(i);   
		rhoEddy  += odtl->eddl->rho->d.at(i)  *dxc.at(i);
        viscEddy += odtl->eddl->dvisc->d.at(i)*dxc.at(i);
			rhoK     += intRhoKi.at(i);  
			rhoJ     += intRhoJi.at(i);
			rhoKK    += odtl->ed->K.at(i)*intRhoKi.at(i);
			rhoJK    += odtl->ed->K.at(i)*abs(intRhoKi.at(i));
			uRhoK.at(0) += intRhoKi.at(i)*odtl->eddl->uvel->d.at(i);
			uRhoK.at(1) += intRhoKi.at(i)*odtl->eddl->vvel->d.at(i);
			uRhoK.at(2) += intRhoKi.at(i)*odtl->eddl->wvel->d.at(i);
			uRhoJ.at(0) += intRhoJi.at(i)*odtl->eddl->uvel->d.at(i);
			uRhoJ.at(1) += intRhoJi.at(i)*odtl->eddl->vvel->d.at(i);
			uRhoJ.at(2) += intRhoJi.at(i)*odtl->eddl->wvel->d.at(i);
	}

	/////////// Particle coupling
	double A;
	double S;
    vector<double> P(3,0.0);
	vector<double> T(3,0.0);
	vector<double> M(3,0.0);
	vector<double> Q(3,0.0);
	vector<double> deltaE(3,0.0);
	
	A = rhoK/rhoJ;
    S = 0.5*(A*A+1.0)*rhoKK - A*rhoJK;
    for(i=0; i<3; i++) {
		M.at(i) = partMomSrc.at(i)/rhoJ;
		P.at(i) = uRhoK.at(i) - A * uRhoJ.at(i) - M.at(i)*rhoJK + M.at(i)*A*rhoKK;
		T.at(i) = 0.5*M.at(i)*M.at(i)*rhoKK - M.at(i)*uRhoJ.at(i);
        Q.at(i) = 0.25*P.at(i)*P.at(i)/S - T.at(i);
    }

    if(odtl->odtp->probType != "SHEARFLOW"){
        deltaE.at(0) = odtl->odtp->A_param*(0.5*Q.at(1)+0.5*Q.at(2)-Q.at(0))-partEnergSrc.at(0);
        deltaE.at(1) = odtl->odtp->A_param*(0.5*Q.at(0)+0.5*Q.at(2)-Q.at(1))-partEnergSrc.at(1);
        deltaE.at(2) = odtl->odtp->A_param*(0.5*Q.at(1)+0.5*Q.at(0)-Q.at(2))-partEnergSrc.at(2);
        eKinEddy = KK/(eLength2*VolE)*(Q.at(0)+Q.at(1)+Q.at(2)-partEnergSrc.at(1));
    }
    else{
        if(odtl->ed->eddyType == 1){
            deltaE.at(0) = odtl->odtp->A_param*(Q.at(1)-Q.at(0));
            deltaE.at(1) = odtl->odtp->A_param*(Q.at(0)-Q.at(1))-partEnergSrc.at(1);
            deltaE.at(2) = 0.0;
            eKinEddy = KK/(eLength2*VolE)*(Q.at(0)+Q.at(1)-partEnergSrc.at(1));
        }
        else if(odtl->ed->eddyType == 2){
            deltaE.at(0) = 0.0;
            deltaE.at(1) = odtl->odtp->A_param*(Q.at(2)-Q.at(1))-partEnergSrc.at(1);
            deltaE.at(2) = odtl->odtp->A_param*(Q.at(1)-Q.at(2));
            eKinEddy = KK/(eLength2*VolE)*(Q.at(1)+Q.at(2)-partEnergSrc.at(1));
        }
        else{
            if(Q.at(1) < -partEnergSrc.at(1)) 
                return false;
            deltaE.at(0) = odtl->odtp->A_param*(Q.at(2)-Q.at(0)); 
            deltaE.at(1) = -partEnergSrc.at(1);
            deltaE.at(2) = odtl->odtp->A_param*(Q.at(0)-Q.at(2));
            eKinEddy = KK/(eLength2*VolE)*(Q.at(0)+Q.at(2)-partEnergSrc.at(1));
        }
    }
	
    odtl->ed->invTauEddy = 0.0;
	double Etot = eKinEddy - odtl->ed->eViscPenalty + odtl->ed->ePeEddy;
	double energyDiff = odtl->ed->Etot - Etot;
	//cout << "\n Energy diff: " << energyDiff << " Etot before: " << odtl->ed->Etot;
    if(Etot < 0.0) return false;
    
	odtl->ed->invTauEddy = sqrt(2.0*KK/(rhoKK*VolE*eLength2) * Etot);

    //////////// Compute Kernel Coefficients
    if(odtl->odtp->probType != "SHEARFLOW"){
        for(i=0; i<3; i++) {
            if(4*S*(T.at(i) - deltaE.at(i)) > P.at(i)*P.at(i))
                return false;
            odtl->ed->cCoef.at(i) = 0.5/S * (-1*P.at(i) + (P.at(i)>0 ? 1.0 : -1.0)
                 * sqrt( P.at(i)*P.at(i) - 4*S*(T.at(i) - deltaE.at(i))));
            odtl->ed->bCoef.at(i) = - M.at(i) - odtl->ed->cCoef.at(i)*A;
        }
    }
    else{
        if(odtl->ed->eddyType == 1){
            if(4*S*(T.at(0) - deltaE.at(0)) > P.at(0)*P.at(0) || 4*S*(T.at(1) - deltaE.at(1)) > P.at(1)*P.at(1))
                return false;
            odtl->ed->cCoef.at(0) = 0.5/S * (-P.at(0) + (P.at(0)>0 ? 1.0 : -1.0)
                * sqrt( P.at(0)*P.at(0) - 4*S*(T.at(0) - deltaE.at(0))));
            odtl->ed->cCoef.at(1) = 0.5/S * (-P.at(1) + (P.at(1)>0 ? 1.0 : -1.0)
                * sqrt( P.at(1)*P.at(1) - 4*S*(T.at(1) - deltaE.at(1))));
            odtl->ed->cCoef.at(2) = 0.0;
        }
        else if(odtl->ed->eddyType == 2){
            if(4*S*(T.at(1) - deltaE.at(1)) > P.at(1)*P.at(1) || 4*S*(T.at(2) - deltaE.at(2)) > P.at(2)*P.at(2))
                return false;
            odtl->ed->cCoef.at(0) = 0.0;
            odtl->ed->cCoef.at(1) = 0.5/S * (-P.at(1) + (P.at(1)>0 ? 1.0 : -1.0)
                * sqrt( P.at(1)*P.at(1) - 4*S*(T.at(1) - deltaE.at(1))));
            odtl->ed->cCoef.at(2) = 0.5/S * (-P.at(2) + (P.at(2)>0 ? 1.0 : -1.0)
                * sqrt( P.at(2)*P.at(2) - 4*S*(T.at(2) - deltaE.at(2))));
        }
        else{
            if(4*S*(T.at(0) - deltaE.at(0)) > P.at(0)*P.at(0) || 4*S*(T.at(2) - deltaE.at(2)) > P.at(2)*P.at(2) || 4*S*(T.at(1) - deltaE.at(1)) > P.at(1)*P.at(1))
                return false;
            odtl->ed->cCoef.at(0) = 0.5/S * (-P.at(0) + (P.at(0)>0 ? 1.0 : -1.0)
                * sqrt( P.at(0)*P.at(0) - 4*S*(T.at(0) - deltaE.at(0))));
            odtl->ed->cCoef.at(1) = 0.5/S * (-P.at(1) + (P.at(1)>0 ? 1.0 : -1.0)
                * sqrt( P.at(1)*P.at(1) - 4*S*(T.at(1) - deltaE.at(1))));
            odtl->ed->cCoef.at(2) = 0.5/S * (-P.at(2) + (P.at(2)>0 ? 1.0 : -1.0)
                * sqrt( P.at(2)*P.at(2) - 4*S*(T.at(2) - deltaE.at(2))));
        }
    }

    for(i=0; i<3; i++)
        odtl->ed->bCoef.at(i) = -A*odtl->ed->cCoef.at(i);

return true;
}
//////////////////////////////////////////////////////////////////////
/*  Compute eddyTau with regard to the particle source terms
* Step 1: Compute kernel integrals
* Step 2: Consider energy and momentum exchange between particle-gas phase
* Step 3: Compute kernel coefficients
*/
bool particle::set_kernel_coefficients(){

	vector<double> dxc(odtl->eddl->ngrd);

    double         KK=0;                                  // equivalent of the 4/27 fac
	double         rhoKK=0;
    double         rhoJK=0;
	double         rhoK=0;
	double         rhoJ=0;
    vector<double> uRhoJ(3,0.0), uRhoK(3,0.0);

    double          eKinEddy     = 0.0;
    double   		eViscPenalty = 0.0;
    double         	ePartEddy    = 0.0;

    vector<double> 	intRhoKi(odtl->eddl->ngrd);
    vector<double> 	intRhoJi(odtl->eddl->ngrd);
	double 			intRho = 0;
    int            	i;

	double rhoEddy=0;
    double viscEddy=0;
    double eLength2 = odtl->ed->eddySize*odtl->ed->eddySize;

   	//////////// cell sizes
	
	dxc = odtl->ed->dxc;
	double VolE = accumulate(dxc.begin(), dxc.end(), 0.0);
 
	///////////// Fill in integral quantities

	//---------- rhoK, rhoJ, UrhoK, UrhoJ
	
	for(i=0; i<odtl->eddl->ngrd; i++) {
		intRhoKi.at(i)  = odtl->ed->K.at(i)*odtl->eddl->rho->d.at(i)*dxc.at(i);                  
       	intRhoJi.at(i)  = abs(intRhoKi.at(i));
		KK       += odtl->ed->K.at(i)*odtl->ed->K.at(i)*dxc.at(i);   
		rhoEddy  += odtl->eddl->rho->d.at(i)  *dxc.at(i);
        viscEddy += odtl->eddl->dvisc->d.at(i)*dxc.at(i);
		rhoK     += intRhoKi.at(i);  
		rhoJ     += intRhoJi.at(i);
		rhoKK    += odtl->ed->K.at(i)*intRhoKi.at(i);
		rhoJK    += odtl->ed->K.at(i)*abs(intRhoKi.at(i));
		uRhoK.at(0) += intRhoKi.at(i)*odtl->eddl->uvel->d.at(i);
		uRhoK.at(1) += intRhoKi.at(i)*odtl->eddl->vvel->d.at(i);
		uRhoK.at(2) += intRhoKi.at(i)*odtl->eddl->wvel->d.at(i);
		uRhoJ.at(0) += intRhoJi.at(i)*odtl->eddl->uvel->d.at(i);
		uRhoJ.at(1) += intRhoJi.at(i)*odtl->eddl->vvel->d.at(i);
		uRhoJ.at(2) += intRhoJi.at(i)*odtl->eddl->wvel->d.at(i);
	}

	/////////// Particle coupling
	double A;
	double S;
    vector<double> P(3,0.0);
	vector<double> T(3,0.0);
	vector<double> M(3,0.0);
	vector<double> Q(3,0.0);
	vector<double> deltaE(3,0.0);
	
	A = rhoK/rhoJ;
    S = 0.5*(A*A+1.0)*rhoKK - A*rhoJK;
    for(i=0; i<3; i++) {
		M.at(i) = partMomSrc.at(i)/rhoJ;
		P.at(i) = uRhoK.at(i) - A * uRhoJ.at(i) - M.at(i)*rhoJK + M.at(i)*A*rhoKK;
		T.at(i) = 0.5*M.at(i)*M.at(i)*rhoKK - M.at(i)*uRhoJ.at(i);
        Q.at(i) = 0.25*P.at(i)*P.at(i)/S - T.at(i);
    }
 
    if(odtl->odtp->probType != "SHEARFLOW"){
        deltaE.at(0) = odtl->odtp->A_param*(0.5*Q.at(1)+0.5*Q.at(2)-Q.at(0))-partEnergSrc.at(0);
        deltaE.at(1) = odtl->odtp->A_param*(0.5*Q.at(0)+0.5*Q.at(2)-Q.at(1))-partEnergSrc.at(1);
        deltaE.at(2) = odtl->odtp->A_param*(0.5*Q.at(1)+0.5*Q.at(0)-Q.at(2))-partEnergSrc.at(2);
    }
    else{
        if(odtl->ed->eddyType == 1){
            deltaE.at(0) = odtl->odtp->A_param*(Q.at(1)-Q.at(0));
            deltaE.at(1) = odtl->odtp->A_param*(Q.at(0)-Q.at(1))-partEnergSrc.at(1);
            deltaE.at(2) = 0.0;
        }
        else if(odtl->ed->eddyType == 2){
            deltaE.at(0) = 0.0;
            deltaE.at(1) = odtl->odtp->A_param*(Q.at(2)-Q.at(1))-partEnergSrc.at(1);
            deltaE.at(2) = odtl->odtp->A_param*(Q.at(1)-Q.at(2));
        }
        else{
            if(Q.at(1) < -partEnergSrc.at(1)) return false;
            deltaE.at(0) = odtl->odtp->A_param*(Q.at(2)-Q.at(0));
            deltaE.at(1) = -partEnergSrc.at(1);
            deltaE.at(2) = odtl->odtp->A_param*(Q.at(0)-Q.at(2));
        }
    }

    //////////// Compute Kernel Coefficients
    if(odtl->odtp->probType != "SHEARFLOW"){
        for(i=0; i<3; i++) {
            odtl->ed->cCoef.at(i) = 0.5/S * (-1*P.at(i) + (P.at(i)>0 ? 1.0 : -1.0)
                 * sqrt( P.at(i)*P.at(i) - 4*S*(T.at(i) - deltaE.at(i))));
            odtl->ed->bCoef.at(i) = - M.at(i) - odtl->ed->cCoef.at(i)*A;
        }
    }
    else{
        if(odtl->ed->eddyType == 1){
            if(4*S*(T.at(0) - deltaE.at(0)) > P.at(0)*P.at(0) || 4*S*(T.at(1) - deltaE.at(1)) > P.at(1)*P.at(1))
                return false;
            odtl->ed->cCoef.at(0) = 0.5/S * (-P.at(0) + (P.at(0)>0 ? 1.0 : -1.0)
                * sqrt( P.at(0)*P.at(0) - 4*S*(T.at(0) - deltaE.at(0))));
            odtl->ed->cCoef.at(1) = 0.5/S * (-P.at(1) + (P.at(1)>0 ? 1.0 : -1.0)
                * sqrt( P.at(1)*P.at(1) - 4*S*(T.at(1) - deltaE.at(1))));
            odtl->ed->cCoef.at(2) = 0.0;
        }
        else if(odtl->ed->eddyType == 2){
            if(4*S*(T.at(1) - deltaE.at(1)) > P.at(1)*P.at(1) || 4*S*(T.at(2) - deltaE.at(2)) > P.at(2)*P.at(2))
                return false;
            odtl->ed->cCoef.at(0) = 0.0;
            odtl->ed->cCoef.at(1) = 0.5/S * (-P.at(1) + (P.at(1)>0 ? 1.0 : -1.0)
                * sqrt( P.at(1)*P.at(1) - 4*S*(T.at(1) - deltaE.at(1))));
            odtl->ed->cCoef.at(2) = 0.5/S * (-P.at(2) + (P.at(2)>0 ? 1.0 : -1.0)
                * sqrt( P.at(2)*P.at(2) - 4*S*(T.at(2) - deltaE.at(2))));
        }
        else{
            if(4*S*(T.at(0) - deltaE.at(0)) > P.at(0)*P.at(0) || 4*S*(T.at(2) - deltaE.at(2)) > P.at(2)*P.at(2) || 4*S*(T.at(1) - deltaE.at(1)) > P.at(1)*P.at(1))
                return false;
            odtl->ed->cCoef.at(0) = 0.5/S * (-P.at(0) + (P.at(0)>0 ? 1.0 : -1.0)
                * sqrt( P.at(0)*P.at(0) - 4*S*(T.at(0) - deltaE.at(0))));
            odtl->ed->cCoef.at(1) = 0.5/S * (-P.at(1) + (P.at(1)>0 ? 1.0 : -1.0)
                * sqrt( P.at(1)*P.at(1) - 4*S*(T.at(1) - deltaE.at(1))));
            odtl->ed->cCoef.at(2) = 0.5/S * (-P.at(2) + (P.at(2)>0 ? 1.0 : -1.0)
                * sqrt( P.at(2)*P.at(2) - 4*S*(T.at(2) - deltaE.at(2))));
        }
    }

    for(i=0; i<3; i++)
        odtl->ed->bCoef.at(i) = -A*odtl->ed->cCoef.at(i);

return true;
}
////////////////////////////////////////////////////////////////////
/** Apply new particle properties after two way coupling and if eddy is accepted 
*/
void particle::applyPartProp(){
    vector<double> y0 = yPos;
	if(!periodicEddy){
		for (int i = 0; i < odtl->odtp->Nparticle; i++){
			if(linePos.at(i) >= odtl->ed->leftEdge && linePos.at(i) <= odtl->ed->rightEdge){
            	vvelP.at(i) = vvelPnew.at(i);
            	yPos.at(i) = yPosnew.at(i);

            	uvelPnew.at(i) = 0;
            	vvelPnew.at(i) = 0;
            	wvelPnew.at(i) = 0;
            	yPosnew.at(i) = 0;

                //---------------------Set odtline position 
                linePos.at(i) = yPos.at(i);
        	}
    	}
	}
	else{
		double edCycle = odtl->odtp->domainLength; //to catch all particles on the negative side of the domain
		for (int i = 0; i < odtl->odtp->Nparticle; i++){
			pJumpL = 0.0;
            pJumpU = 0.0;
			if(yPosnew.at(i) > odtl->odtp->domainLength*0.5){
                pJumpL = odtl->odtp->domainLength;
                pJumpU = odtl->uvel->d.at(odtl->ngrd-1) - odtl->uvel->d.at(0);
                //pJumpU = odtl->odtp->probType == "SHEARFLOW" ? odtl->io->initParams["Srate"].as<double>()*odtl->io->params["domainLength"].as<double>() : 0.0;
            }
            else if( yPosnew.at(i)*linePos.at(i) < 0.0){
                    pJumpU = odtl->uvel->d.at(odtl->ngrd-1) - odtl->uvel->d.at(0);
            }
			if(linePos.at(i) >= odtl->ed->leftEdge || linePos.at(i) <= odtl->ed->rightEdge-edCycle){
			   vvelP.at(i) = vvelPnew.at(i);
			   yPos.at(i) = yPosnew.at(i)-pJumpL;
               if(linePos.at(i)*yPos.at(i) < 0.0)
                   uvelP.at(i) = uvelP.at(i) + (yPos.at(i)<0.0 ? -1.0 : 1.0)*pJumpU; 
			}
            uvelPnew.at(i) = 0;
            vvelPnew.at(i) = 0;
            wvelPnew.at(i) = 0;
            yPosnew.at(i) = 0;
            zPosnew.at(i) = 0;
            //---------------------Set odtline position 
            linePos.at(i) = yPos.at(i);
		}
	}
	getDisplacement(y0, yPos);
}
/////////////////////////////////////////////////////////////////////
/** Compute solution of eddy crossing time of particle (Newton process)
* @param velE \input eddy velocity
* @param velG \input gas velocity
* @param velP \input particle velocity
* @param fDivTauP \input quotient of corrector factor and response time
* @param g \input gravity component
* @param partPos \input initial particle position
* @param edgePos \input position of eddy box edge
* @param T1 \input starting time of eddy
* @param T2 \input end time of eddy
*/
double particle::eddyCrossingTime(double velE, double velG, double velP, double fDivTauP, double g, double partPos, double edgePos, double T1, double T2){
	
	double Tcross;
	double f1 = pPosRelToEdge(velE, velG, velP, fDivTauP, g, partPos, edgePos, T1);
   	double f2 = pPosRelToEdge(velE, velG, velP, fDivTauP, g, partPos, edgePos, T2);
   	double Tn, fn, fp;
	
	double T2old = T2;
	double f1old = f1;
	double f2old = f2;

   	for(int it=1; it<=100; it++) {
       	Tn = (T1+T2)*0.5;
       	fn = pPosRelToEdge(velE, velG, velP, fDivTauP, g, partPos, edgePos, Tn);
		if(abs(fn) < 1E-6)
			break;
       	(f1*fn < 0.0) ? (T2 = Tn) : (T1=Tn);
   	}

	if(abs(fn) < 1E-6)
		 Tcross = Tn;
	else{
    	//----------- clean up solution with Newton
   		double Told = Tn;	
		double fold = fn;
   		int maxIt=1000;
		for(int it=1; it<=maxIt; it++) {
			if(isnan(fn))
				cout << endl << "Fail!!!";
			fp = ddt_pPosRelToEdge(velE, velG, velP, fDivTauP, g, Tn);
			Tn = Tn - fn/fp;
			fn = pPosRelToEdge(velE, velG, velP, fDivTauP, g, partPos, edgePos, Tn);
			if(abs(fn) < 1E-6 ) {
        		Tcross = Tn;
    			break;
    	    }
    		if(it==maxIt){
    			cout << endl << "# warning in eddy::getParticleUVWYafterEddy, not converged" << endl;
                exit(0);
            }
		}
	}
return Tcross;
}
//-----------------------------help functions------------------------------------
//analytical solution of particle trajectory
double particle::pPosRelToEdge(double Evel, double Gvel, double Pvel, double fDivTauP, double AG, double pPos0, double edgePos0, double T) {
return pPos0 + (Gvel+AG/fDivTauP)*T + (Pvel-Gvel-AG/fDivTauP)/fDivTauP*(1.0-exp(-T*fDivTauP)) - (edgePos0 + Evel*T);
}
//derivative of analytical solution of particle trajectory
double particle::ddt_pPosRelToEdge(double Evel, double Gvel, double Pvel, double fDivTauP, double AG, double T) {
return (Gvel+AG/fDivTauP) +(Pvel-Gvel-AG/fDivTauP)*exp(-T*fDivTauP) - (Evel);
}
///////////////////////////////////////////////////////////////////
/** Compute response time tauP and corrector f and set quotient
* @param i \index particle index
* @param uvel \input x-component of particle velocities
* @param vvel \input y-component of particle velocities
* @param vvel \input z-component of particle velocities
*/
void particle::set_fDivTauP(int i, double uvelP, double vvelP, double wvelP, double uvel, double vvel, double wvel){
	
	double gasFreeMeanPath; 		// free mean path of particles (maybe in odtp)
	double Cc;                    	// Cunningham slip factor
    double magRelVel;             	// magnitude relative velocity
    double pi = 3.14159265359;

	gasFreeMeanPath = 68E-9;      	// Ambient pressure
	//Cc = 1 + 2*gasFreeMeanPath/diamP.at(i)*(1.257+0.4*exp(-1.1*diamP.at(i)/(2*gasFreeMeanPath)));
	Cc = 1 ; 
    tauP.at(i) = dens0P.at(i)*diamP.at(i)*diamP.at(i)*Cc/18/odtl->odtp->kvisc0;
    magRelVel = sqrt(pow(uvelP - uvel ,2) + pow(vvelP - vvel ,2) + pow(wvelP - wvel ,2));
    double ReP = odtl->odtp->rho0*magRelVel*diamP.at(i)/odtl->odtp->kvisc0;

   	if(ReP <= 100.)
   		f.at(i) = 1.+ 0.15*pow(ReP,0.687);
   	else{
    	cout << "\nWARNING: LARGE REYNOLDS NUMBER > 1000 " << endl;
		f.at(i) = 1.+ 0.15*pow(ReP,0.687);
	}
   	fDivTauP.at(i) = f.at(i)/tauP.at(i);
}
///////////////////////////////////////////////////////////////////
/** Initialize particle phase statistics
*/
void particle::initStats(){

	//@iStat -> odt-solv->iStat
	//@nPart -> nPart

	int nPart = odtl->odtp->Nparticle;
	
	uMean = vector<double>(nPart, 0.0);
	vMean = vector<double>(nPart, 0.0);
	wMean = vector<double>(nPart, 0.0);
	
    // edstat(j,k,l,i): i = particle index
    //                  j = 0.5*del_U_eddy^2, 0.5*del_U_flow^2, 0.5*del_U_all^2, 0.5*del_U_adpt^2
    //                  k = u, v, w 
    //                  l = average period
    edstat      = vector<vector<vector<vector<double> > > >
                    (4, vector<vector<vector<double> > >
                    (3, vector<vector<double> >
                    (odtl->odtp->nStat, vector<double>(nPart, 0.0) ) ) );
    // cstat(m,k,l,i):  i = particle index
    //                  m = u, u^2   
    //                  k = u, v, w, 
    //                  l = average period
    cstat       = vector<vector<vector<vector<double> > > >
                    (2, vector<vector<vector<double> > >
                    (3, vector<vector<double> >
                    (odtl->odtp->nStat, vector<double>(nPart, 0.0) ) ) );

    //@ctime -> odt-solv->ctime;
    // oldVars(k,i)     i = particle index
    // ctime(l):        l = average period
    //                  k = u, v, w
    oldVars     = vector<vector<double> > (3, vector<double>(nPart, 0.0) );

}
////////////////////////////////////////////////////////////////////
/** This function gathers statistics of velocity and its variance in each diffusion time step.
*	@param tStep \input the time step for the gathering
*/
void particle::statsTime(const double tStep){
	int iStat = odtl->solv->iStat;	

	// u velocity
    if(odtl->odtp->probType != "SHEARFLOW"){
        for(int i=0; i<odtl->odtp->Nparticle; i++){
            cstat[0][0][iStat-1][i] += tStep *     uvelP.at(i);
            cstat[1][0][iStat-1][i] += tStep * pow(uvelP.at(i),2);
            cstat[0][1][iStat-1][i] += tStep *     vvelP.at(i);
            cstat[1][1][iStat-1][i] += tStep * pow(vvelP.at(i),2);
            cstat[0][2][iStat-1][i] += tStep *     wvelP.at(i);
            cstat[1][2][iStat-1][i] += tStep * pow(wvelP.at(i),2);
        }
    }
    else{
        double shearRate = odtl->io->initParams["Srate"].as<double>();
        for(int i=0; i<odtl->odtp->Nparticle; i++){
            cstat[0][0][iStat-1][i] += tStep *     (uvelP.at(i)-yPos.at(i)*shearRate);
            cstat[1][0][iStat-1][i] += tStep * pow(uvelP.at(i)-yPos.at(i)*shearRate,2);
            cstat[0][1][iStat-1][i] += tStep *     vvelP.at(i);
            cstat[1][1][iStat-1][i] += tStep * pow(vvelP.at(i),2);
            cstat[0][2][iStat-1][i] += tStep *     wvelP.at(i);
            cstat[1][2][iStat-1][i] += tStep * pow(wvelP.at(i),2);
        }
    }
}
///////////////////////////////////////////////////////////////////////////////
/**
 *  statsSetOld saves the current odtline to be used in 
 *  statsChange to calculate the difference generated through an eddy event or a 
 *  diffusion step. 
 */
void particle::statsSetOld(){

    for (int i=0; i<odtl->odtp->Nparticle; i++)
        oldVars[0][i] = 0.5*pow(uvelP.at(i)-uMean[i], 2);

    for (int i=0; i<odtl->odtp->Nparticle; i++)
        oldVars[1][i] = 0.5*pow(vvelP.at(i)-vMean[i], 2);

    for (int i=0; i<odtl->odtp->Nparticle; i++)
        oldVars[2][i] = 0.5*pow(wvelP.at(i)-wMean[i], 2);
}
///////////////////////////////////////////////////////////////////////////////
/**
 *  statsChange saves and gathers the changes generated through a process (eddy 
 *  event, diffusion, ...) into edstat.
 *
 *  @param jj    \input switch: 0: eddy event  2: diffusion  4: all processes
 */
void particle::statsChange(const int &jj){

	int iStat = odtl->solv->iStat;
    // jj = 0  ==>  differences caused by eddies
    // jj = 1  ==>  differences caused by diffusion
    // jj = 2  ==>  differences caused by eddies and diffusion (all changes) 
    // jj = 3  ==>  differences caused by adaption after diffusion and eddies

    for (int i=0; i<odtl->odtp->Nparticle; i++){
        edstat[0+jj][0][iStat-1][i] += 0.5*pow(uvelP.at(i)-uMean[i], 2) - oldVars[0][i];
    }

    for (int i=0; i<odtl->odtp->Nparticle; i++){
        edstat[0+jj][1][iStat-1][i] += 0.5*pow(vvelP.at(i)-vMean[i], 2) - oldVars[1][i];
    }

    for (int i=0; i<odtl->odtp->Nparticle; i++){
        edstat[0+jj][2][iStat-1][i] += 0.5*pow(wvelP.at(i)-wMean[i], 2) - oldVars[2][i];
    }
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
void particle::statsOutput(){

    cout << "Start subroutine part->statsOutput" << endl;

    int    N  = 0;
	int nPart = odtl->odtp->Nparticle;

    // deklaration and initialisation of needed arrays
    vector< vector< vector< vector<double> > > > eavg;
    vector< vector< vector< vector<double> > > > cavg;

    // eavg(j,k,l,i):   i = particle index
    //                  j = del_0.5*U^2_eddy, del_0.5*U^2_flow, del_0.5*U^2_all,
    //                      del_0.5*U^2_adpt
    //                  k = u, v, w
    //                  l = average periode
    eavg = vector<vector<vector<vector<double> > > > (4,
           vector<vector<vector<double> > > (3,
           vector<vector<double> > (odtl->odtp->nStat,
           vector<double>(nPart ,0.0) ) ) );
    // cavg(m,k,l,i):     i = coordinate
    //                    m = u, u^2
    //                    k = u, v, w
    //                    l = average periode
    cavg = vector<vector<vector<vector<double> > > > (2,
           vector<vector<vector<double> > > (3,
           vector<vector<double> > (odtl->odtp->nStat,
           vector<double>(nPart ,0.0) ) ) );

    vector< vector<double> > tke;
    vector< vector<double> > d_diff;
    vector< vector<double> > d_eddy;
    vector< vector<double> > d_all;
    vector< vector<double> > d_adpt;
    vector< vector<double> > bal;

	// tke(i,l):        i = particle index
    //                  l = average periode
    // tv, d, ta, p have the same as tke
    tke    = vector<vector<double> > (odtl->odtp->nStat, vector<double>(nPart ,0.0) );
    d_diff = vector<vector<double> > (odtl->odtp->nStat, vector<double>(nPart ,0.0) );
    d_eddy = vector<vector<double> > (odtl->odtp->nStat, vector<double>(nPart ,0.0) );
    d_all  = vector<vector<double> > (odtl->odtp->nStat, vector<double>(nPart ,0.0) );
    d_adpt = vector<vector<double> > (odtl->odtp->nStat, vector<double>(nPart ,0.0) );
    bal    = vector<vector<double> > (odtl->odtp->nStat, vector<double>(nPart ,0.0) );

	cout << "particle :: statsOutput :: calculate cavg" << endl;
    for(int i = 0; i < odtl->odtp->nStat; i++){
        for(int k = 0; k < 2; k++){
            for(int l = 0; l < 3; l++){
                for(int j = 0; j < nPart; j++){
                    cavg[k][l][i][j] = cstat[k][l][i][j] / odtl->solv->ctime[i];
                }
            }
        }
    }

    cout << "particle :: statsOutput :: calculate eavg" << endl;
    for(int k = 0; k < odtl->odtp->nStat; k++){
        for(int i = 0; i < 4; i++){
            for(int l = 0; l < 3; l++){
                for(int j = 0; j < nPart; j++){
                    eavg[i][l][k][j] = edstat[i][l][k][j] / odtl->solv->ctime[k];
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
        for (int j = 0; j < nPart; j++){
            tke[k][j]  = 0.5*cavg[1][0][k][j] ; // <u'^2>
            tke[k][j] += 0.5*cavg[1][1][k][j] ; // <v'^2>
            tke[k][j] += 0.5*cavg[1][2][k][j] ; // <w'^2>

            d_eddy[k][j]  = eavg[0][0][k][j] + eavg[0][1][k][j] + eavg[0][2][k][j]; // eddy event
            d_diff[k][j]  = eavg[1][0][k][j] + eavg[1][1][k][j] + eavg[1][2][k][j]; // diffusion
            d_all[k][j]   = eavg[2][0][k][j] + eavg[2][1][k][j] + eavg[2][2][k][j]; // eddy + diffusion
            d_adpt[k][j]  = eavg[3][0][k][j] + eavg[3][1][k][j] + eavg[3][2][k][j]; // mesh adaption
        }
    }

    for (int k = 0; k < odtl->odtp->nStat; k++){
        for (int j = 0; j < nPart; j++){
            bal[k][j] = d_eddy[k][j] + d_diff[k][j] + d_adpt[k][j];
        }
    }
	// start data output

    cout << "particle :: statsOutput :: start writing output" << endl;
    // -- writing ASCII output
    N = nPart;

    for (int k = 0; k < odtl->odtp->nStat; k++){
        stringstream ss1; string s1;
        ss1.clear(); s1.clear();
        ss1 << setfill('0') << setw(5) << k; ss1 >> s1;
        string name = "../data/test/particle_data/data-average-output_" + s1 + ".dat";
        ofstream file_avg(name.c_str(),ios::app);
        assert(file_avg.is_open());
        file_avg << "# grid points = " << N;
        file_avg << "\n# Domain Size = " << odtl->solv->LdomainS;
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
                 <<    " " << setw(14) << odtl->odtp->tEnd * (k+1) / odtl->odtp->nStat;
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
                 <<     "  bal                    ";
        file_avg << scientific;
        file_avg << setprecision(16);
		for (int j = 0; j < N; j++){
            file_avg << endl;
            file_avg << setw(25) << j
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
                     << setw(25) << bal[k][j];
        }
        file_avg.close();
    }
    cout << "particle :: statsOutput :: output in ASCII files finished" << endl;

	return;

}
////////////////////////////////////////////////////////////////////
/**
*  In Progress!!!!!!!!!!!!!
*/
void particle::collision(const vector<double> &dxc,  double dt){

	if(dt == 0)
		return;

	int nCells = dxc.size();
	int nPart = odtl->odtp->Nparticle;
	vector<int> partInCell(nCells+1, 0);
	vector<int> index(nCells+1, 0);
	vector<int> Xref(nPart, 0);
	double pi = 3.14159265359;
	
    /// Build list of number of particles per cell
	for(int i=0; i<nPart; i++){
		if(activeP.at(i) == false)
			partInCell.at(nCells)++;
		else
			partInCell.at(indP.at(i))++;
	}
	/// Build index list as cumulative sum of number of particles in cell
	int m = 0; 
	for(int i=0; i<nCells+1; i++){
		index.at(i) = m;
		m += partInCell.at(i); 
	}
	/// Build cross-reference list
	vector<int> tmp(nCells+1,0);
	for(int i=0; i<nPart; i++){
		if(activeP.at(i) == true){
			Xref.at(index.at(indP.at(i)) + tmp.at(indP.at(i))) = i;
			tmp.at(indP.at(i)) += 1;
		}
		else{
			Xref.at(index.at(nCells) + tmp.at(nCells)) = i;
            tmp.at(nCells) += 1;
		}
	}	
	/// Select and compute collisions per cell
	for(int j=0; j<nCells; j++){
		int number = partInCell.at(j);
		if(number > 1){
			/// Determine number of candidate collision pairs
			if(odtl->odtp->Lspatial){
				double meanU =  0;
				// Compute mean axial velocity of particles in cell
				for(int i=0; i<number;i++){
					meanU += uvelP.at(Xref.at(index.at(j)+i));
				}
				meanU /= number;
				dt = dt/meanU; // Transform spatial coordinate in temporal			
			}
			double collParam = 0.5*NPartPerParcel.at(1)*pi*diamP.at(1)*diamP.at(1)*dt/dxc.at(j); // !Only for homogeneous parcels!
			double select = collParam*number*number*crmax;
			int nsel = floor(select);
			double crm = crmax;
			for(int isel = 0; isel < nsel; isel++){
				int k = floor(odtl->rand->getRand()*number);
				int kk = int(ceil(k + odtl->rand->getRand()*(number-1)))%number;
				int ip1 = Xref.at(k+index.at(j));
				int ip2 = Xref.at(kk+index.at(j));
				/// Compute relative velocity
				double cr = sqrt(pow(uvelP.at(ip1) - uvelP.at(ip2) ,2) + pow(vvelP.at(ip1) - vvelP.at(ip2) ,2) + pow(wvelP.at(ip1) - wvelP.at(ip2) ,2));
				if(cr > crm) crm = cr;
				/// Accept or reject candidate
				if(cr/crmax > odtl->rand->getRand()){
					col += 1;  // Collision counter
					double cos_th = 1-2*odtl->rand->getRand(); // Collision angle theta
					double sin_th = sqrt(1-cos_th*cos_th);
					double sqrt(1-cos_th*cos_th);
					double phi = 1-2*odtl->rand->getRand(); // Collision angle phi
					/// Apply post-collision velocities (hard-sphere)
					uvelP.at(ip1) = 0.5*(uvelP.at(ip1)+uvelP.at(ip2)) + 0.5*cr*cos_th;
					vvelP.at(ip1) = 0.5*(vvelP.at(ip1)+vvelP.at(ip2)) + 0.5*cr*sin_th*cos(phi);
					wvelP.at(ip1) = 0.5*(wvelP.at(ip1)+wvelP.at(ip2)) + 0.5*cr*sin_th*sin(phi);
					uvelP.at(ip2) = 0.5*(uvelP.at(ip1)+uvelP.at(ip2)) - 0.5*cr*cos_th;	
					vvelP.at(ip2) = 0.5*(vvelP.at(ip1)+vvelP.at(ip2)) - 0.5*cr*sin_th*cos(phi);
                    wvelP.at(ip2) = 0.5*(wvelP.at(ip1)+wvelP.at(ip2)) - 0.5*cr*sin_th*sin(phi);
					collisionData(j);
 				}
				crmax = crm;
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////
/** Particle output 
*@param fname \input Name for output file
*@param time \input Time step of Simulation
*(Very depending on test case)
*/ 
void particle::collisionData(int j){
	
	stringstream ss1;
	string s1;
        
	string name = "../data/test/particle_data/particle_collision.dat";
    ofstream write_output(name.c_str(),ios::app);
    assert(write_output.is_open());

    if(col == 1){
    	write_output << "\n# Total Number of particle = "   << odtl->odtp->Nparticle ;
		write_output << "\n# Max. relative velocity = "   << crmax ;
		if(odtl->odtp->Lspatial) write_output << "\n# No." << setw(19) << "axial position" << setw(19) << "radial position";
		else write_output << "\n# No." << setw(19) << "time" << setw(19) << "radial position";
        write_output << "\n#"<< endl;
	}
	write_output << scientific;
    write_output << setprecision(10);
    write_output << col << setw(19) << odtl->mimx->tend << setw(19) << odtl->pos->d.at(j);	
	write_output << "\n";
    write_output.close();
}

////////////////////////////////////////////////////////////////////
/** Destructor
*/ 
particle::~particle() {}

