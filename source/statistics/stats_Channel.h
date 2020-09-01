/**
 * @file stats.h
 * Header file for classes stats 
 */

#ifndef STATS_CHANNEL_H
#define STATS_CHANNEL_H

#include <string>
#include <vector>

using namespace std;

///////////////////////////////////////////////////////////////////////////////

/** Class implementing statistics gathering on the odt line such as means and variances.
 *  
 *  @author Juan Medina; based on previous versions of David Lignell and Zoltan Jozefik
 */

class stats_Channel : public stats {

    ////////////////////// DATA MEMBERS /////////////////////

    public:

        std::vector<double> uMean;         ///< mean u velocities avg over periods
        std::vector<double> vMean;         ///< mean v velocities
        std::vector<double> wMean;         ///< mean w velocities

        std::vector<double> uMsqr;         ///< mean squares fluctuation u avg over periods
        std::vector<double> vMsqr;         ///< mean squares v fluctuations
        std::vector<double> wMsqr;         ///< mean squares w fluctuations
	
	std::vector<double> uvFlucMean;				///< mean of the uv fluctuations' product

    ////////////////////// MEMBER FUNCTIONS  /////////////////////

        virtual void initStats();
	
	virtual void onlineTimeAvg(std::vector<double> &vecBase,
                      std::vector<double> &vecToAdd,
                      const double &t, const double &dt);
	
	virtual void updateMeans(const double &dt);
    
        virtual void outputProperties(std::string outputStatsFName);
	
    private:

	void onlineSqTimeAvg(std::vector<double> &vecBase, std::vector<double> &vecMean,
                      std::vector<double> &vecToAdd, const double &t, const double &dt);
      
	void crossVelFlucAvg(std::vector<double> &vecBase, std::vector<double> &vecMean1, std::vector<double> &vecMean2,std::vector<double> &vecInst1,
			     std::vector<double> &vecInst2, const double &t, const double &dt);
      
    ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public:

	stats_Channel(){ odtl=0; }
        virtual void statsConst(odtline *p_odtl, double Ld=1.0, int npts=100, std::string dataPath = "./");               // constructor
        virtual ~stats_Channel(){ odtl=0; }

};

#endif
