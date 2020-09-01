/**
 * @file stats.cc
 * Header file for class stats
 */

#include "stats.h"
#include "odtline.h"
#include "stats_Channel.h"
#include <iomanip>
#include <fstream>
#include <cmath>          // fabs

using namespace std;


///////////////////////////////////////////////////////////////////////////////

/**Constructor function (The one to use).
 * Stats has its own uniform grid and does minimal statistics gathering.
 * It computes the mean and mean square values of velocity and mixture fraction.
 * The stats grid goes from 0 to Ldomain, even for periodic domains, so there
 * is some wrapping to do when the domain is periodic.
 * vTrans is just a dummy vector for transfering from the odt grid to the
 * stats grid.
 *
 * @param odtpp \input parameters object to set pointer.
 * @param Ld    \input domain length.
 * @param npts  \input number of evenly spaced stats points along domain.
 */
void stats_Channel::statsConst(odtline *p_odtl, double Ld, int npts, std::string dataPath) {

    odtl = p_odtl;
    
    Ldomain = Ld;
    
    dataDir = dataPath;
    
    ngrd    = npts;
    ngrdf   = ngrd+1;
    pos     = vector<double>(ngrd,  0.0);
    posf    = vector<double>(ngrdf, 0.0);
      
    vTrans  = vector<double>(ngrd, 0.0);
    
    double dx = Ldomain / ngrd;
    pos[0] = (-Ldomain + dx)/2.0;
    for(int i=1; i<ngrd; i++) 
        pos.at(i) = pos.at(i-1) + dx;
    posf[0] = -Ldomain/2.0;
    for(int i=1; i<ngrdf; i++)
        posf[i] = posf[i-1] + dx;
    
    statsTime			= 0.0;
    
    inTime			= false;    
    ngrdTime			= 1000;				//size for interpolating in time. changes as more time pnts are added.
    
    statsTime_Vec.resize	(ngrdTime, 0.0);

}


///////////////////////////////////////////////////////////////////////////////

/**Call this once the odtline has been initialized. This is done so that the 
 * stats are not initially biased by the starting value in the constructor.
 * This could probably be moved into the constructor too, if desired.
 *
 * @param odtl \input odtline object to get stats for.
 */
void stats_Channel::initStats() {
    
    uMean   = vector<double>(ngrd, 0.0);
    vMean   = vector<double>(ngrd, 0.0);
    wMean   = vector<double>(ngrd, 0.0);

    uMsqr   = vector<double>(ngrd, 0.0);
    vMsqr   = vector<double>(ngrd, 0.0);
    wMsqr   = vector<double>(ngrd, 0.0);
    
    uvFlucMean	= vector<double>(ngrd, 0.0);
    
    odtGrd2statGrd(odtl->posf->d, odtl->uvel->d, 0.0);  // modified to omit pJump non-existent input Parameter for periodic domains
    uMean 	= vTrans;
    
    odtGrd2statGrd(odtl->posf->d, odtl->vvel->d, 0.0);  // modified to omit pJump non-existent input Parameter for periodic domains
    vMean 	= vTrans;

    odtGrd2statGrd(odtl->posf->d, odtl->wvel->d, 0.0);  // modified to omit pJump non-existent input Parameter for periodic domains
    wMean 	= vTrans;

    Ldomain 	= odtl->Ldomain();

    for(int i=0; i<ngrd; i++) {
        uMsqr[i]      = 0.0;  // initialization refers to uMsqr as the square of the fluctuations
        vMsqr[i]      = 0.0;
        wMsqr[i]      = 0.0;
	    uvFlucMean[i] = 0.0;
    }
    
}


///////////////////////////////////////////////////////////////////////////////

/**Compute current values of mean variables.
 * You could also pass in an anyline, but then need to get the right position
 * in the props array to get uvel, etc.
 * Currently this will give a Reynolds average, not a Favre Average
 *
 * @param odtl \input odtline for computing stats.
 * @param t \input current realization time.
 * @param dt \input time step size.
 */
void stats_Channel::updateMeans(const double &dt) {
  
    // x-component****************************************************************************************
    odtGrd2statGrd(odtl->posf->d, odtl->uvel->d, 0.0);		// modified to omit pJump non-existent input Parameter for periodic domains
    onlineTimeAvg( uMean, vTrans, statsTime, dt);				// statsTime= current time in statistics
    onlineSqTimeAvg(uMsqr, uMean, vTrans, statsTime, dt);  
    
    vector<double> uVelCopy(vTrans);					// copy of the u velocity vector 
    
    // y-component****************************************************************************************
    odtGrd2statGrd(odtl->posf->d, odtl->vvel->d, 0.0);	// modified to omit pJump non-existent input Parameter for periodic domains
    onlineTimeAvg( vMean, vTrans, statsTime, dt);
    onlineSqTimeAvg(vMsqr, vMean, vTrans, statsTime, dt);
        
    crossVelFlucAvg(uvFlucMean, uMean, vMean, uVelCopy, vTrans, statsTime, dt);
    
    
    // z-component****************************************************************************************
    odtGrd2statGrd(odtl->posf->d, odtl->wvel->d, 0.0);	// modified to omit pJump non-existent input Parameter for periodic domains
    onlineTimeAvg( wMean, vTrans, statsTime, dt);
    onlineSqTimeAvg(wMsqr, wMean, vTrans, statsTime, dt);
    
    
    statsTime +=dt;
    
    uVelCopy.clear();
    
}



///////////////////////////////////////////////////////////////////////////////

/**Adds an instantaneous profile to a mean profile. Weighting each
 * by their respective ages: the current mean has a weight of the
 * simulation age, while the profile to add gets a weight of dt.
 * \cond
 * <x> = ( x1*dt1 + x2*dt2 + x3*dt3 ) / (t=dt1+dt2+dt3)
 * <x> = ( <x>*t + x4*dt4 ) / (t + dt4)
 *
 * NOTE: BELOW IS JUST LATEX STYLE FORMULAS UNTIL PARAMTERS
 * \endcond
 *
 * \f[
 *      \langle x \rangle = \frac{ x_1 d t_1 + x_2 d t_2 + x_3 d t_3 }{(t = d t_1+d t_2+d t_3)}
 * \f]
 * \f[
 *      \langle x \rangle = \frac{ \langle x \rangle t + x_4 d t_4 }{ t + d t_4}
 * \f]
 *
 *
 * @param vecBase  \inout mean variable to update.
 * @param vecToAdd \input variable to augment the mean.
 * @param t        \input current realization time.
 * @param dt       \input time step size.
 */
void stats_Channel::onlineTimeAvg(std::vector<double> &vecBase,
                     std::vector<double> &vecToAdd,
                     const double &t, const double &dt) {


    for(int i=0; i<vecToAdd.size(); i++) {
        vecBase[i] = (vecBase[i] * t + vecToAdd[i]*dt)/(t+dt);
    }

}

///////////////////////////////////////////////////////////////////////////////
/**Like onlineTimeAvg, but squared fluctuation quantities instead of the value itself.
 * @param vecBase  \inout mean variable to update.
 * @param vecMean  \input velocity mean.
 * @param vecToAdd \input instantaneous velocity.
 * @param t        \input current realization time.
 * @param dt       \input time step size.
 */
void stats_Channel::onlineSqTimeAvg(std::vector<double> &vecBase, std::vector<double> &vecMean,
                     std::vector<double> &vecToAdd,
                     const double &t, const double &dt) {


    for(int i=0; i<vecToAdd.size(); i++)
        vecBase[i] = (vecBase[i] * t + (vecToAdd[i]-vecMean[i])*(vecToAdd[i]-vecMean[i])*dt)/(t+dt);

}


///////////////////////////////////////////////////////////////////////////////

/**Cross velocity fluctuation product
 * @param vecBase  \inout mean variable to update.
 * @param vecMean1 \input first mean velocity component in the cross product
 * @param vecMean2 \input second mean velocity component in the cross product
 * @param vecInst1 \input first instantaneous velocity component in the cross product
 * @param vecInst2 \input second instantaneous velocity component in the cross product
 * @param t        \input current realization time.
 * @param dt       \input time step size.
 */
void stats_Channel::crossVelFlucAvg(std::vector<double> &vecBase, std::vector<double> &vecMean1, std::vector<double> &vecMean2,
                     std::vector<double> &vecInst1, std::vector<double> &vecInst2,
                     const double &t, const double &dt) {


    for(int i=0; i<vecInst2.size(); i++) {
        vecBase[i] = (vecBase[i] * t + (vecInst1[i]-vecMean1[i])*(vecInst2[i]-vecMean2[i])*dt)/(t+dt);
    }

}


///////////////////////////////////////////////////////////////////////////////

/** Outputs the properties of the stats.   
 *
 *  @param fname \input output file name.
 */
void stats_Channel::outputProperties(std::string outputStatsFName) {

   string fname = dataDir + outputStatsFName;
   
   ofstream ofile(fname.c_str()); 
   if(!ofile) 
       cout << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;
   
   ofile << "# grid points = "   << ngrd;
   ofile << "\n# Domain Size = " << Ldomain;
   if(fabs(posf[ngrd]-posf[0]- Ldomain) > 1.0E-6)
       ofile << "\n# last posf-first posf != Ldomain, last posf = " << posf[ngrd];
   ofile << "\n# 1_pos             "
         <<     "2_posf            "
         <<     "3_uMean           "
         <<     "4_vMean           " 
         <<     "5_wMean           "
         <<     "6_uMsqr           "
         <<     "7_vMsqr           "
         <<     "8_wMsqr           "
	 <<	"9_uvFlucM	   "
         ;
   ofile << scientific;
   ofile << setprecision(10);
   for(int i=0; i<ngrd; i++) 
       ofile << endl 
           << setw(19) << pos[i] 
           << setw(19) << posf[i]
           
           << setw(19) << uMean[i]
           << setw(19) << vMean[i]
           << setw(19) << wMean[i]
           
           << setw(19) << uMsqr[i]
           << setw(19) << vMsqr[i]
           << setw(19) << wMsqr[i]
           
           << setw(19) << uvFlucMean[i]
           ;

   ofile.close();

}

///////////////////////////////////////////////////////////////////////////////
