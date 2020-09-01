/**
 * @file stats.h
 * Header file for classes stats 
 */

#ifndef STATS_H
#define STATS_H

#include <string>
#include <vector>

class odtline;

using namespace std;

///////////////////////////////////////////////////////////////////////////////

/** Class implementing statistics gathering on the odt line such as means and variances.
 *  
 *  @author Juan Medina; based on previous versions of David Lignell and Zoltan Jozefik
 */

class stats {

    ////////////////////// DATA MEMBERS /////////////////////

    public:

	odtline		    *odtl;		///< Pointer to #odtline object of current realization
	
        int                 ngrd;          	///< Number of grid points. Constant grid assumed here
        int                 ngrdf;         	///< Number of grid faces. Constant grid assumed here
        double              Ldomain;       	///< stat domain size (same as odt)
        
        std::string	    dataDir;

        double		    locIndex;	   	/// index that points to actualLdomain start in pos vector
        double		    statsTime;	   	/// time in stats interval. = 0 after each stats period. Global time is not reset
        int		    ngrdTime;      	/// # of stat grid points in time is different than in space
        bool		    inTime;  	   	/// for odtGrd2statGrd, size of grid is different for Time than for space. boolean for space grd or time grd.    
        
        std::vector<double> pos;           	///< vector of cell center positions
        std::vector<double> posf;          	///< vector of cell face  positions
        std::vector<double> vTrans;        	///< a vector to tranfer odt grid to stat grid
        std::vector<double> statsTime_Vec;  	///< contains 'real time' when stats are updated on stats grid.
        
    ////////////////////// MEMBER FUNCTIONS  /////////////////////

	virtual void initStats(){};
	
	virtual void onlineTimeAvg(std::vector<double> &vecBase,
                      std::vector<double> &vecToAdd,
                      const double &t, const double &dt){};
	
	virtual void updateMeans(const double &dt){};
    
        virtual void outputProperties(std::string outputStatsFName){};
 	
	virtual void odtGrd2statGrd(std::vector<double> odtposf, std::vector<double> odtVec, double pJump=0.0);		// default method
	
	virtual void outputEddyAcceptedData(double &time_implemented);		// default method

    private:


    ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public:

	stats(){ odtl=0; }
        virtual void statsConst(odtline *p_odtl, double Ld=1.0, int npts=100, std::string dataPath = "./"){};               // constructor
        virtual ~stats(){ odtl=0; }

};

#endif
