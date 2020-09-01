
#include "kernel.h"
#include "odtline.h"
#include <cmath>
#include <numeric>   //accumulate
#include <algorithm> // min_element
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>

///////////////////////////////////////////////////////////////////////////////
/** eddy initialization function
 *
 * @param p_odtl  \input set odtline pointer with.
 */

void kernel::init(odtline *p_odtl) {
    odtl = p_odtl;

    nkernels = 0;
    rhoU2.resize(3,0.0);
    DeltaU2TM.resize(3,0.0);
    injectionTime = odtl->odtp->trst;
    meanTime = odtl->odtp->T11*odtl->odtp->L11/odtl->odtp->domainLength;
    sampleKernelTime();
}

///////////////////////////////////////////////////////////////////////////////
/** Sample an eddy size from the eddy size distribution.
 * This could be as simple as a tophat, but this function is
 * more accurate, and more efficient.
 */
void kernel::kernelEnergyInjection(double time) {

    stringstream ss1;
    string       s1;

    sampleKernelTime();
    kernelSize = odtl->odtp->L11;

    sampleKernelPosition();
    ss1.clear();  ss1 << setfill('0') << setw(4) << nkernels; ss1 >> s1;
    //odtl->io->writeDataFile("odt_"+s1+"_preKern.dat", time);
    if(LperiodicKernel) {
        *odtl->io->ostrm << endl << "#   periodic kernel ";
        ss1.clear();  ss1 << setfill('0') << setw(4) << nkernels; ss1 >> s1;
        double cycleDistance = odtl->cyclePeriodicLine(iEnd);
        double bkp_rightEdge = rightEdge;
        rightEdge = leftEdge + kernelSize;
        iStart = odtl->linePositionToIndex(leftEdge,  true, 6);
        iEnd   = odtl->linePositionToIndex(rightEdge, false, 7);
        integrateTKE(0);   //-> -Get TKE
	    tripMap();
	    iStart = odtl->linePositionToIndex(leftEdge,  true, 6);
        iEnd   = odtl->linePositionToIndex(rightEdge, false, 7);
        integrateTKE(1);  //-> Get Delta TKE caused by TM
	    fillKernel();
	    computeKernelCoeff();
	    applyKernels();
	    iEnd   = odtl->ngrd-1;;
        odtl->backCyclePeriodicLine(cycleDistance);
        rightEdge = bkp_rightEdge;
    }
    else{
        ss1.clear();  ss1 << setfill('0') << setw(4) << nkernels; ss1 >> s1;
        integrateTKE(0);  //-> -Get TKE
        tripMap();
	    iStart = odtl->linePositionToIndex(leftEdge,  true, 6);
        iEnd   = odtl->linePositionToIndex(rightEdge, false, 7);
        integrateTKE(1);  //-> Get Delta TKE caused by TM
	    fillKernel();
	    computeKernelCoeff();
	    applyKernels();
        //odtl->io->writeDataFile("odt_"+s1+"_postKern.dat", time);
    }
    odtl->io->writeDataFile("odt_"+s1+"_kernel.dat", time);
    nkernels++;
}

//////////////////////////////////////////////////////////////////////////////
/** Sample an eddy size from the eddy size distribution.
 * This could be as simple as a tophat, but this function is
 * more accurate, and more efficient.
 */
void kernel::sampleKernelTime() {
    
    double r    = odtl->rand->getRand();
    kernelTime  = -1*meanTime*log(1-r); 
    injectionTime += kernelTime;

}

///////////////////////////////////////////////////////////////////////////////
/** Uniformly sample the kernel position on the line.  For periodic domains the
 * position can be anywhere.  For nonperiodic domains, the position
 * is from 0 to the end where end is the domain size - the kernel size
 * (since the kernel has to fit in the domain).  This also means that
 * the kernel position is the left edge of the kernel event.
 * For periodic domains, rightEdge is greater than leftEdge (even for kernels that
 * wrap the domain, and may be outside the range of the base odtline.)
 */
void kernel::sampleKernelPosition() {

    if(!odtl->odtp->Lperiodic) {
        leftEdge  = odtl->rand->getRand() * (odtl->Ldomain() - kernelSize) + odtl->posf->d.at(0);
        rightEdge = leftEdge + kernelSize;
        if(leftEdge  < odtl->posf->d.at(0))          leftEdge  = odtl->posf->d.at(0);
        if(rightEdge > odtl->posf->d.at(odtl->ngrd)) rightEdge = odtl->posf->d.at(odtl->ngrd);
    }
    else {
        LperiodicKernel = false;             // reset the value
        leftEdge  = odtl->rand->getRand() * odtl->Ldomain() + odtl->posf->d.at(0);
        rightEdge = leftEdge + kernelSize;
        if(rightEdge > (odtl->Ldomain() + odtl->posf->d.at(0))){
            LperiodicKernel = true;
            rightEdge = odtl->posf->d.at(0) + rightEdge - odtl->posf->d.at(odtl->ngrd);
        }
    }

    iStart = odtl->linePositionToIndex(leftEdge,  true, 4);
    iEnd = odtl->linePositionToIndex(rightEdge, false, 5);
}

/////////////////////////////////////////////////////////////////////////////
/** Triplet map the line
 *  This version will adjust the TM so that the three
 *      segments are spaced evenly in cylindrical and spherical flows.
 *      (planar = no difference).
 *  This can be converted back to the usual equal volume segments by making fracVleft = 1/3 etc.
 *  @param line \inout either the odtl or the eddl
 *  @param iS   \input starting cell
 *  @param iE   \input ending cell
 *  (Note, this could be done by the lv objects separately.)
 *  Applies to an eddy line, or the full odt line.
 */
void kernel::tripMap() {

    int iS = iStart;
    int iE = iEnd;

    int negrd    = iE-iS+1;
    int negrd2   = 2*negrd;
    int negrd3   = 3*negrd;

    double fracVleft;
    double fracVmidl;
    double fracVrght;
    double pm1;          // +1 or -1 factor to change sign
    double C = odtl->odtp->cCoord;
    double invC = 1.0/C;

    vector<double>   posf0;
    int k, i, ip, im, ie, iw;

    //--------- Set fracVlft, fracVmid, fracVrgt

    double r1 = (2*leftEdge+rightEdge)/3;
    double r2 = (leftEdge+2*rightEdge)/3;

    pm1 =leftEdge*r1 < 0 ? -1.0 : 1.0;             // handles cells that split x=0
    double Vleft = abs(pow(abs(leftEdge), C) - pm1*pow(abs(r1), C));

    pm1 =r1*r2 < 0 ? -1.0 : 1.0;                   // handles cells that split x=0
    double Vmidl = abs(pow(abs(r1), C) - pm1*pow(abs(r2), C));

    pm1 =r2*rightEdge < 0 ? -1.0 : 1.0;            // handles cells that split x=0
    double Vrght = abs(pow(abs(r2), C) - pm1*pow(abs(rightEdge), C));

    double Vtot = Vleft + Vmidl + Vrght;
    fracVleft = Vleft / Vtot;
    fracVmidl = Vmidl / Vtot;
    fracVrght = Vrght / Vtot;

    //---------

    pos0 = odtl->pos->d;
    posf0 = odtl->posf->d;

    //----------- Grab cell "volume" array (delta(x^c)) in the eddy region

    dxc.resize(negrd,0.0);

    i=0; ip=1;
    pm1 = odtl->posf->d.at(iS+ip)*leftEdge < 0 ? -1.0 : 1.0;                     // handles cells that split x=0
    dxc.at(i) = abs(pow(abs(odtl->posf->d.at(iS+ip)), C) - pm1*pow(abs(leftEdge), C));

    for(int i=1, ip=2; i<negrd-1; i++, ip++) {
        pm1 = odtl->posf->d.at(iS+ip)*odtl->posf->d.at(iS+i) < 0 ? -1.0 : 1.0;   // handles cells that split x=0
        dxc.at(i) = abs(pow(abs(odtl->posf->d.at(iS+ip)), C) - pm1*pow(abs(odtl->posf->d.at(iS+i) ), C));
    }

    i=negrd-1;
    pm1 = odtl->posf->d.at(iS+i)*rightEdge < 0 ? -1.0 : 1.0;                     // handles cells that split x=0
    dxc.at(i) = abs(pow(abs(rightEdge), C) - pm1*pow(abs(odtl->posf->d.at(iS+i)), C));

    //---------- insert cells to be filled with values
 
    for(k=0; k<odtl->v.size(); k++)
        odtl->v.at(k)->d.insert(odtl->v.at(k)->d.begin()+iE+1, negrd2, 0.0);

    odtl->ngrd += negrd2;
    odtl->ngrdf = odtl->ngrd+1;
    iE += negrd2;

    //---------- update face positions

    if(C==1 || leftEdge > 0.0) {      // case: planar, or eddy on the right of the centerline
        // posf.at(iS) is the same
        odtl->posf->d.at(iS+1) = pow(pow(leftEdge,C) + fracVleft*dxc.at(0), invC);
        for(ie=2, iw=1, i=1; ie<=negrd; ie++, iw++, i++)
            odtl->posf->d.at(iS+ie) = pow(pow(odtl->posf->d.at(iS+iw),C) + fracVleft*dxc.at(i), invC);
        for(ie=negrd+1,iw=negrd,i=negrd-1; ie<=negrd2; ie++, iw++, i--)
            odtl->posf->d.at(iS+ie) = pow(pow(odtl->posf->d.at(iS+iw),C) + fracVmidl*dxc.at(i), invC);
        for(ie=negrd2+1,iw=negrd2,i=0; ie<negrd3; ie++, iw++, i++)
            odtl->posf->d.at(iS+ie) = pow(pow(odtl->posf->d.at(iS+iw),C) + fracVrght*dxc.at(i), invC);
        // odtl->posf->d.at(iS+negrd3) is the same due to the inserted cells before.
    }
    else if(rightEdge < 0.0) {                       // case: eddy on the left of centerline

        // posf.at(iS) is the same
        odtl->posf->d.at(iS+1) = -pow(pow(abs(leftEdge),C) - fracVleft*dxc.at(0), invC);
        for(ie=2, iw=1, i=1; ie<=negrd; ie++, iw++, i++)
            odtl->posf->d.at(iS+ie) = -pow(pow(abs(odtl->posf->d.at(iS+iw)),C) - fracVleft*dxc.at(i), invC);
        for(ie=negrd+1,iw=negrd,i=negrd-1; ie<=negrd2; ie++, iw++, i--)
            odtl->posf->d.at(iS+ie) = -pow(pow(abs(odtl->posf->d.at(iS+iw)),C) - fracVmidl*dxc.at(i), invC);
        for(ie=negrd2+1,iw=negrd2,i=0; ie<negrd3; ie++, iw++, i++)
            odtl->posf->d.at(iS+ie) = -pow(pow(abs(odtl->posf->d.at(iS+iw)),C) - fracVrght*dxc.at(i), invC);
        // odtl->posf->d.at(iS+negrd3) is the same due to the inserted cells before.

    }
    else {                                           // case: eddy splits the centerline
        // Note, this can directly replace the above two cases!

        double dmb;

        // posf.at(iS) is the same
        pm1 = leftEdge < 0 ? -1 : 1;
        dmb = pow(abs(leftEdge),C) + pm1*fracVleft*dxc.at(0);
        odtl->posf->d.at(iS+1) = dmb < 0.0 ? pow(-dmb, invC) : pm1*pow(dmb, invC);
        for(ie=2, iw=1, i=1; ie<=negrd; ie++, iw++, i++) {
            pm1 = odtl->posf->d.at(iS+iw) < 0 ? -1 : 1;
            dmb = pow(abs(odtl->posf->d.at(iS+iw)),C) + pm1*fracVleft*dxc.at(i);
            odtl->posf->d.at(iS+ie) = dmb < 0.0 ? pow(-dmb, invC) : pm1*pow(dmb, invC);
        }
        for(ie=negrd+1,iw=negrd,i=negrd-1; ie<=negrd2; ie++, iw++, i--) {
            pm1 = odtl->posf->d.at(iS+iw) < 0 ? -1 : 1;
            dmb = pow(abs(odtl->posf->d.at(iS+iw)),C) + pm1*fracVmidl*dxc.at(i);
            odtl->posf->d.at(iS+ie) = dmb < 0.0 ? pow(-dmb, invC) : pm1*pow(dmb, invC);
        }
        for(ie=negrd2+1,iw=negrd2,i=0; ie<negrd3; ie++, iw++, i++) {
            pm1 = odtl->posf->d.at(iS+iw) < 0 ? -1 : 1;
            dmb = pow(abs(odtl->posf->d.at(iS+iw)),C) + pm1*fracVrght*dxc.at(i);
            odtl->posf->d.at(iS+ie) = dmb < 0.0 ? pow(-dmb, invC) : pm1*pow(dmb, invC);
        }
        // odtl->posf->d.at(iS+negrd3) is the same due to the inserted cells before.
    }

	//---------- update center positions
    odtl->pos->setVar();

    for(k=0; k<odtl->v.size(); k++) {
        if(odtl->v.at(k)->var_name=="pos" || odtl->v.at(k)->var_name=="posf") continue;
        vector<double> v0 = odtl->v.at(k)->d;
        int iD;
        iD = iS == 0 ? odtl->ngrd-1 : iS-1;
        odtl->v.at(k)->d.at(iS)   = odtl->v.at(k)->d.at(iD)+(odtl->pos->d.at(iS)-odtl->pos->d.at(iD))
		*(v0.at(iS)-v0.at(iD))/(pos0.at(iS)-odtl->pos->d.at(iD));
        int iw = iS;
	    double ie = iS+1;
	    for(i=0; i<negrd3-1; i++) {
            while(odtl->pos->d.at(iS+i) > pos0.at(ie)){
                if(ie != pos0.size()-1){
                    iw++;
		            ie++;
                }
                else
                    break;
            }
	        odtl->v.at(k)->d.at(iS+i) = v0.at(iw)+(odtl->pos->d.at(iS+i)-pos0.at(iw))
		        *(v0.at(ie)-v0.at(iw))/(pos0.at(ie)-pos0.at(iw)); 
        }
        iD = iE == odtl->ngrd-1 ? 0 : iE+1; 
	    odtl->v.at(k)->d.at(iE) = odtl->v.at(k)->d.at(iE-1)+(odtl->pos->d.at(iE)-odtl->pos->d.at(iE-1))
            *(odtl->v.at(k)->d.at(iD)-odtl->v.at(k)->d.at(iE-1))/(odtl->pos->d.at(iD)-odtl->pos->d.at(iE-1)); 
    }	
}
///////////////////////////////////////////////////////////////////////////////
/** Fill velocity kernel K (used also for \fun{J=|K|})
*/
void kernel::fillKernel() {

    ngrd     = iEnd - iStart + 1;	
    int nseg = ngrd / 3;
    K.resize(ngrd);

    //---------- 1st Segment
    K.at(0) = 0.5*(odtl->posf->d.at(iStart+1) + leftEdge) - pos0.at(iStart);  // pos is 0.5*(faces), but we want pos=0.5*(face+eddy_edge)
    for(int i=1, j=1; i<nseg; i++, j++)                           //    pos0 was already updated this way in tripmap
        K.at(i) = odtl->pos->d.at(iStart+i) - pos0.at(iStart+j);

    //---------- Second Segment
    for(int i=nseg*2-1, j=0; i>=nseg; i--, j++)
        K.at(i) = odtl->pos->d.at(iStart+i) - pos0.at(iStart+j);

    //---------- Third Segment
    for(int i=nseg*2, j=0; i<nseg*3-1; i++, j++)
        K.at(i) = odtl->pos->d.at(iStart+i) - pos0.at(iStart+j);
    int i=nseg*3-1, j=nseg-1;
    K.at(i) = 0.5*(rightEdge + odtl->posf->d.at(iStart+i)) - pos0.at(iStart+j); // pos is 0.5*(faces), but we want pos=0.5*(face+eddy_edge)
}

///////////////////////////////////////////////////////////////////////////////
/** Fill velocity kernel K (used also for \fun{J=|K|})
*/
void kernel::computeKernelCoeff() {
    double         rhoKK=0;
    vector<double> uRhoK(3,0.0);
    vector<double> intRhoKi(ngrd);
  
    vector<double> injTKE(3,0.0);

    //----------- set dxc

    odtl->mesher->setGridDxc(odtl, dxc, odtl->odtp->cCoord);

    //----------------------

    for(int i=0; i<ngrd; i++) {
        intRhoKi.at(i)  = K.at(i)*odtl->rho->d.at(iStart+i)*dxc.at(iStart+i);
        if(!odtl->odtp->Lspatial) {
            rhoKK          += K.at(i)*intRhoKi.at(i);
            uRhoK.at(0)    += intRhoKi.at(i)*odtl->uvel->d.at(iStart+i);
            uRhoK.at(1)    += intRhoKi.at(i)*odtl->vvel->d.at(iStart+i);
            uRhoK.at(2)    += intRhoKi.at(i)*odtl->wvel->d.at(iStart+i);
        }
        else {       // mass flux, not mass
            rhoKK          += odtl->uvel->d.at(iStart+i) * K.at(i)*intRhoKi.at(i);
            uRhoK.at(0)    += odtl->uvel->d.at(iStart+i) * intRhoKi.at(i)*odtl->uvel->d.at(iStart+i);
            uRhoK.at(1)    += odtl->uvel->d.at(iStart+i) * intRhoKi.at(i)*odtl->vvel->d.at(iStart+i);
            uRhoK.at(2)    += odtl->uvel->d.at(iStart+i) * intRhoKi.at(i)*odtl->wvel->d.at(iStart+i);
        }
    }
    
    double Cke = odtl->odtp->Prod*pow(odtl->odtp->T11,3)/pow(kernelSize,2);
    E = 2*odtl->odtp->rho0*pow(kernelSize,11.0/3.0)*pow(odtl->odtp->T11, -3.0)*pow(Cke, 4.0/3.0);
    for(int i=0; i<3; i++) 
        injTKE.at(i) = 1.0/3.0*E + DeltaU2TM.at(i);
   
    dCoef.resize(3,0.0);

    for(int i=0; i<3; i++){
        double dm = uRhoK.at(i) < 0.0 ? -1.0 : 1.0;
        if(0.0 > pow(uRhoK.at(i),2)+2*injTKE.at(i)*rhoKK)
            dCoef.at(i) = 0.0;
        else
            dCoef.at(i) = 1.0/rhoKK*(-uRhoK.at(i)+dm*sqrt(pow(uRhoK.at(i),2)+2*injTKE.at(i)*rhoKK));
    }
}

///////////////////////////////////////////////////////////////////////////////
/** Apply kernels K and J to the velocity profile.
 * This is called after the kernel coefficients is computed in eddyTau
*/
void kernel::applyKernels() {

    ////////// update velocity profiles

    int npts = iEnd-iStart+1;

    if(npts != K.size()){
        cout << endl << "ERROR, WRONG SIZE IN kernel::applyKernels: " << npts << " " << K.size() << endl;
        exit(0);
    }

    if(odtl->odtp->Lspatial) {
        bool Lflag = false;                  // if kernel results in negative (or < umin) uvel for Lspatial, then don't apply kernels.
        for(int i=0; i<npts; i++)            //     a better approach would be to find the limiting alpha and use that.
            if ( odtl->uvel->d.at(i+iStart) + dCoef.at(0)*K.at(i) < odtl->odtp->umin_spatial ) {
                Lflag = true;
                break;
            }
        if(!Lflag) {
            for(int i=0; i<npts; i++) {
                odtl->uvel->d.at(i+iStart) += dCoef.at(0)*K.at(i) ;
                odtl->vvel->d.at(i+iStart) += dCoef.at(1)*K.at(i) ;
                odtl->wvel->d.at(i+iStart) += dCoef.at(2)*K.at(i) ;
            }
        }
    }
    else {
        for(int i=0; i<npts; i++) {
            odtl->uvel->d.at(i+iStart) += dCoef.at(0)*K.at(i);
            odtl->vvel->d.at(i+iStart) += dCoef.at(1)*K.at(i);
            odtl->wvel->d.at(i+iStart) += dCoef.at(2)*K.at(i);
        }
    }
}
///////////////////////////////////////////////////////////////////////////////
/** Fill velocity kernel K (used also for \fun{J=|K|})
*/
void kernel::integrateTKE(int ii) {

    if(ii == 0){
        for(int i=0; i<3; i++){
            rhoU2.at(i) = 0.0;
        }
        double intRho;
        int n = iEnd-iStart+2;  
        //----------- set dxc
        odtl->mesher->setGridDxc(odtl, dxc, odtl->odtp->cCoord);
        //----------------------

        for(int i=0; i<odtl->ngrd; i++) {
            intRho  = 0.5*odtl->rho->d.at(i)*dxc.at(i);
            if(!odtl->odtp->Lspatial) {
                rhoU2.at(0)    += intRho*odtl->uvel->d.at(i)*odtl->uvel->d.at(i);
                rhoU2.at(1)    += intRho*odtl->vvel->d.at(i)*odtl->vvel->d.at(i);
                rhoU2.at(2)    += intRho*odtl->wvel->d.at(i)*odtl->wvel->d.at(i);
            }
            else {       // mass flux, not mass
                rhoU2.at(0)    += odtl->uvel->d.at(i)*intRho*odtl->uvel->d.at(i)*odtl->uvel->d.at(i);
                rhoU2.at(1)    += odtl->uvel->d.at(i)*intRho*odtl->vvel->d.at(i)*odtl->vvel->d.at(i);
                rhoU2.at(2)    += odtl->uvel->d.at(i)*intRho*odtl->wvel->d.at(i)*odtl->wvel->d.at(i);
            }
        }
    }
    else{
        vector<double>      rhoU2_0 = rhoU2;
        double intRho;
        int n = iEnd-iStart+2; 
        for(int i=0; i<3; i++){
            rhoU2.at(i) = 0.0;
            DeltaU2TM.at(i) = 0.0;
        }

        //----------- set dxc
        odtl->mesher->setGridDxc(odtl, dxc, odtl->odtp->cCoord);
        //----------------------

        for(int i=0; i<odtl->ngrd; i++) {
            intRho  = 0.5*odtl->rho->d.at(i)*dxc.at(i);
            if(!odtl->odtp->Lspatial) {
                rhoU2.at(0)    += intRho*odtl->uvel->d.at(i)*odtl->uvel->d.at(i);
                rhoU2.at(1)    += intRho*odtl->vvel->d.at(i)*odtl->vvel->d.at(i);
                rhoU2.at(2)    += intRho*odtl->wvel->d.at(i)*odtl->wvel->d.at(i);
            }
            else {       // mass flux, not mass
                rhoU2.at(0)    += odtl->uvel->d.at(i)*intRho*odtl->uvel->d.at(i)*odtl->uvel->d.at(i);
                rhoU2.at(1)    += odtl->uvel->d.at(i)*intRho*odtl->vvel->d.at(i)*odtl->vvel->d.at(i);
                rhoU2.at(2)    += odtl->uvel->d.at(i)*intRho*odtl->wvel->d.at(i)*odtl->wvel->d.at(i);
            }
        }
        for(int i=0; i<3; i++) DeltaU2TM.at(i) = rhoU2_0.at(i) - rhoU2.at(i);
    }

}
