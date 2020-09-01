/**
 * @file stats.cc
 * Header file for class stats
 */

#include "odtline.h"
#include "stats.h"
#include "stats_Channel.h"
#include <iomanip>
#include <fstream>
#include <cmath>          // fabs

using namespace std;


///////////////////////////////////////////////////////////////////////////////

/**Used to transfer a vector on the odt grid to the stats grid.
 * Routine integrates odt variable odtVec, filling vTrans data member.
 * Loop over each stat grid cell and for each cell march along the odt grid.
 * Consider the grids: 
 * \vc{
 * stat:  |           |           |           |
 * odt:   |   |   |   :     |     |           |
 *                  pLeft
 * }
 * \c pLeft refers to the overlap.  Note overlap on left of second stat cell.
 * When start second stat cell, go from pLeft to fifth stat face.
 * This will work with nonuniform stat grids.
 *
 * @param odtposf \input vector of face positions on the odtline.
 * @param odtvec  \input vector of an odtline variable.
 * @param pJump   \input periodic jumps.
 */
void stats::odtGrd2statGrd(std::vector<double> odtposf, 
                           std::vector<double> odtvec, double pJump) {

  
    //stat grdPnts in time and space are different
    bool	Order1st= false;		//linear interpolation 1st order. (old ODT method, or falkos code)
    int		interpol= 2;		//switch: 1=linear interpolation 2nd order, 2=quadratic interpolation 1st order
    
    int 	grdPnts	= ngrd;
    vector<double> posVec(grdPnts, 0.0);	
    posVec=posf;
    
    if(inTime){
      grdPnts	= ngrdTime;
      posVec.resize(grdPnts, 0.0);
      posVec	=statsTime_Vec;
    }
    
    
  
    int    odtngrd = odtposf.size()-1;

    //---------- If odtposf[0] != 0 (e.g, periodic) then rebuild the odt vectors
    //---------- by wrapping around the domain

    if(odtl->odtp->Lperiodic && odtposf[0] != 0.0) {
        
        vector<double> vd;

        int nDoff = static_cast<int> (odtposf[0] / Ldomain);   // # of full domains offset (usually 0)
        if(nDoff != 0) 
            for(int i=0; i<odtposf.size(); i++)
                odtposf[i] -= nDoff * Ldomain;                 // shift domain 

        for(int i=odtngrd-1; i>=0; i--)
            if(odtposf[i] < Ldomain) {
                int ipos = i;                                 // cell that splits hi boundary

                //---------- split cell ipos at the domain boundary

                if(odtposf[ipos+1] != Ldomain) {
                    odtposf.insert(odtposf.begin()+ipos+1, Ldomain);
                    odtvec.insert( odtvec.begin() +ipos+1, odtvec[ipos]);
                    odtngrd++;        
                    ipos++;
                }
                
                //---------- now wrap the posf and variable vecs

                vd = vector<double>(odtvec.begin()+ipos, odtvec.end());
                if(pJump != 0.0) 
                    for(int j=0; j<vd.size(); j++)
                        vd[j] -= pJump;
                odtvec.erase(odtvec.begin()+ipos, odtvec.end());
                odtvec.insert(odtvec.begin(), vd.begin(), vd.end());


                vd = vector<double>(odtposf.begin()+ipos, odtposf.end()-1);
                for(int j=0; j<vd.size(); j++)
                    vd[j] -= Ldomain;
                odtposf.erase(odtposf.begin()+ipos, odtposf.end()-1);
                odtposf.insert(odtposf.begin(), vd.begin(), vd.end());
                odtposf[odtngrd] = Ldomain;
                
                break;
            }
    }

    //---------- Now transfer grids

    int    i, ip;                       // for stat grid
    int    j, jp;                       // for odt grid
    double pLeft;                       // for the overlap

    //---------- transfer grids

    j     = 0;
    jp    = 1;
    pLeft = posVec[0];

    vTrans = vector<double>(grdPnts, 0.0);
    
    if(Order1st && !inTime)
    {
      for(i=0, ip=1; i<grdPnts; i++, ip++) {                     // loop over stat grd
        while( odtposf[jp] < posVec[ip] && jp!=odtngrd) {     // loop over odt grd
            vTrans[i] += odtvec[j] * (odtposf[jp]-pLeft);
            j++;
            jp++;
            pLeft = odtposf[j];
        }
        vTrans[i] += odtvec[j] * (posVec[ip] - pLeft);        // get the leftovers
        pLeft = posVec[ip];
      }
    }
    else if(!Order1st && !inTime)
    {
      // calculate a, b and c for u = a*x^2+b*x+c
      // a = 0 if a linear interpolation is used
      double a, b, c;
      double delx, delxp, delxt, xj;

      delxp = (odtposf[jp+1] - odtposf[j])  /2;
      delx  = (odtposf[jp+2] - odtposf[jp]) /2;
      xj    = (odtposf[j]    + odtposf[jp]) /2;

      if (interpol == 2) {
        a =  (odtvec[jp+1]-odtvec[jp]) / (delx*(delxp+delx))
            +(odtvec[j]-odtvec[jp])   / (delxp*(delxp+delx));
      }
      else
      {
        a = 0.0;
      }
      b =  (odtvec[j]-odtvec[jp+1]) * delxp      / (delx*(delxp+delx))
        +(odtvec[jp]-odtvec[j]) * (delxp+delx) / (delxp*delx)
        -2 * a * xj;
      c = odtvec[j] - a * xj * xj - b * xj;

      for(i=0, ip=1; i<grdPnts; i++, ip++){
        while( odtposf[jp] < posVec[ip] && jp != odtngrd) {
            delxt = (pLeft + odtposf[jp]) / 2;
            vTrans[i] +=  a * (pow(odtposf[jp],3) -pow(pLeft,3)) / 3.0
                         +( b * (odtposf[jp] +pLeft) / 2 + c) * (odtposf[jp] -pLeft);
            j++; jp++;
            pLeft = odtposf[j];
            delx = delxp;
            xj   = (odtposf[j] + odtposf[jp]) /2;
            if (jp == odtngrd) {
                delxp = (odtposf[j]-odtposf[j-2])/2;
                if (interpol == 2){
                    a =  (odtvec[j]-odtvec[j-1]) / (delxp*delx)
                        +(odtvec[j-2]-odtvec[j]) / (delxp*(delxp+delx));
                }
                else{
                    a = 0.0;
                }
                b =  (odtvec[j]-odtvec[j-1]) * (delx+delxp) / (delx*delxp)
                    +(odtvec[j-2]-odtvec[j]) * delx / (delxp*(delx+delxp))
                    -2 * a * xj;
            }
            else {
                delxp = (odtposf[jp+1]-odtposf[j])/2;
                if (interpol ==2) {
                    a =  (odtvec[jp]-odtvec[j])  / (delxp*(delxp+delx))
                        +(odtvec[j-1]-odtvec[j]) / (delx*(delxp+delx));
                }
                else {
                    a = 0.0;
                }
                b =  (odtvec[jp]-odtvec[j])  * delx / (delxp*(delx+delxp))
                    +(odtvec[j]-odtvec[j-1]) * delxp / (delx*(delx+delxp))
                    -2 * a * xj;
            }
            c = odtvec[j] - a * xj * xj - b * xj;
        }
        delxt = (pLeft + posVec[ip]) / 2;
        vTrans[i] +=  a * (pow(posVec[ip],3) -pow(pLeft,3)) / 3.0
                     +( b * (posVec[ip] +pLeft) / 2 + c) * (posVec[ip] -pLeft);
        pLeft = posVec[ip];
      }
    }
    else if(inTime)
    {  
      //cycle through posVec
      for(i = 0; i < grdPnts; i++)
      {
	      
	if( posVec[i] <= odtposf[0] )				//left boundary pnts
	{
	  if( posVec[i] == odtposf[0] )
	    vTrans[i]	= odtvec[0];
	  else
	  {
	    double y2	= odtvec[1];	  double y1	= odtvec[0];
	    double x2	= odtposf[1];	  double x1	= odtposf[0];
	    vTrans[i]	= (y2-y1) / (x2-x1) *( posVec[i]-x1 )+y1;
	  }
	}
	else if( posVec[i] >= odtposf[odtposf.size()-1])	//right boundary
	{
	  if( posVec[i] == odtposf[odtposf.size()-1] )
	    vTrans[i]	= odtvec[odtposf.size()-1];
	  else
	  {
	    double y2	= odtvec[odtposf.size()-1];	  double y1	= odtvec[odtposf.size()-2];
	    double x2	= odtposf[odtposf.size()-1];	  double x1	= odtposf[odtposf.size()-2];
	    vTrans[i]	= (y2-y1) / (x2-x1) *( posVec[i]-x1 )+y1;
	  }
	}
	else
	{
	  j		= 1;
	  while( posVec[i] >= odtposf[j] && j < odtposf.size() )
	    j++;
	  
	  double y2	= odtvec[j];	  double y1	= odtvec[j-1];
	  double x2	= odtposf[j];	  double x1	= odtposf[j-1];
	  vTrans[i]	= (y2-y1) / (x2-x1) *( posVec[i]-x1 )+y1;
	}
      }
            
    }
    

    //---------- normalize
    if( !inTime) 
    {
      for(i=0, ip=1; i<grdPnts; i++, ip++) 
        vTrans[i] /= (posVec[ip]-posVec[i]);
    }
      
}



///////////////////////////////////////////////////////////////////////////////
/** Outputs current eddy data
 * 
 *  @param time_implemented \input current time to output the data
 * 
 */
void stats::outputEddyAcceptedData(double &time_implemented){
  
  //file output name and locaction
  string fname = dataDir + "eddyAcceptedSizes.dat";
  ofstream myfile;
  
  if(odtl->solv->neddies==1)  {
    myfile.open ( fname.c_str(), ios::trunc);
    myfile   <<   "1_num eddies accepted "
	   <<     "2_edddy size          "
	   <<     "3_edddy left edge     "
	   <<     "4_edddy right edge     "
	   <<     "5_num. of eddy trials " 
	   <<     "6_eddy tau 		 "
	   <<     "7_Pe of eddy 	 "
	   <<     "8_Ke of eddy 	 "
	   <<	  "9_VP of eddy		 "
	   <<     "10_timeCurrent 	 "
	   << endl;
    myfile.close(); 
  }

  myfile.open ( fname.c_str(), ios::app);		//ios::trunc		//"../data/ATC_Time.dat"
  myfile<<scientific<<setprecision(10)
	  <<setw(19)<<odtl->solv->neddies
	  <<setw(19)<<odtl->ed->eddySize
	  <<setw(19)<<odtl->ed->leftEdge
	  <<setw(19)<<odtl->ed->rightEdge
	  <<setw(19)<<odtl->solv->iEtrials
	  <<setw(19)<<(1.0/odtl->ed->invTauEddy)
	  <<setw(19)<<odtl->ed->ePE
	  <<setw(19)<<odtl->ed->eKE
	  <<setw(19)<<odtl->ed->eVP
	  <<setw(19)<<time_implemented
	  <<endl;
  myfile.close(); 

}
