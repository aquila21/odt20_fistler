/**
 * @file particle.h
 * Header file for class particle
 */

#ifndef PARTICLE_H
#define PARTICLE_H

#include "odtline.h"
#include <vector>
#include <string>
#include <map>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing particle object
 *  
 *  @author Marco Fistler
 */

class particle {

    public:

    //////////////////// DATA MEMBERS //////////////////////
    odtline        *odtl;
	odtline        *eddl;
	
    //---------- for particle switch 

    bool            particleON;
	double          initTime; /// initializaton time of particle phase for PEI
	vector<double>  yPos;	  /// y coordinate of particle
	vector<double>  zPos;     /// z coordinate of particle
    vector<double>  linePos;  /// odtline position of particle
	vector<double>  eLinePos; /// odtline position of particle during PEI
	vector<double>  indP;     /// index of particle position to odtline 
	vector<double>  uvelP;    /// u velocity of particle
	vector<double>  vvelP;    /// v velocity of particle 
	vector<double>  wvelP;    /// w velocity of particle
    vector<double>  av_uP;    /// average u velocity of particle in cell
    vector<double>  av_vP;    /// average v velocity of particle in cell
    vector<double>  av_wP;    /// average w velocity of particle in cell
    vector<double>  displ;    /// displacement
  	vector<double>  dens0P;   /// initial density of particle 
    vector<double>  diamP;    /// diameter of particle
    vector<double>  massP;    /// mass of particle 
	vector<double>  fDivTauP; /// quotient of corrector factor and response time
	vector<double>  f;        /// corrector function
	vector<double>  tauP;     /// particle response time
	vector<bool>    activeP;  /// active flag (deactivated if reaching domain boundary)
    vector<int>     typeP;    /// tracer(0), inertial(1), ballistic(2)
	vector<double>  gravity;  /// gravity vector	
	vector<double>  pDistr;   /// particle distribution over domain

    int call;				  /// call of output function
 	vector<double> uvelG;     /// u velocity of gas on particle position
    vector<double> vvelG;     /// v velocity of gas on particle position
    vector<double> wvelG;     /// w velocity of gas on particle position	

	// paramter for particle-eddy interaction
	vector<double> uEddy;
	vector<double> vEddy;
	vector<double> wEddy;
	vector<double> deltaLinePos;///  fluid particle displacement due to TM
	vector<double> randomNo;	/// random number for TM

    // momentum source terms for two-way coupling (diffusion eq.)
 	vector<int>    NPartPerParcel; /// number of "real" particles per parcel (particle cloud)
    vector<double> sumWeights;   /// right hand side source term for u momentum equation
     
    // momentum source terms for two-way coupling (eddy interaction)
            
    vector<double> partMomSrc;   	// 
	vector<double> partEnergSrc;    //

	vector<double> yPosnew;		//
	vector<double> zPosnew;		//
	vector<double> linePosnew;  //
	vector<double> uvelPnew; 	// 
    vector<double> vvelPnew; 	// 
    vector<double> wvelPnew;	//

	bool periodicEddy;          // flag for a periodic eddy	
	double pJumpU;              // periodic velocity jump value
    double pJumpL;				// periodic position jump value

	// variable for particle phase statistics

	vector<double> uMean;         ///< mean u velocities
    vector<double> vMean;         ///< mean v velocities
    vector<double> wMean;         ///< mean w velocities
	vector< vector< vector< vector<double> > > >  edstat;
    // data accumulation for calculation of budget terms     
    vector< vector< vector< vector<double> > > >  cstat;
    // second variable for data accumulation                 
    vector< vector<double> >                      oldVars;
    // variable for saving the old state to calculate the 
	// differenc caused by an eddy or a diffusion step

	// variables for collision detection

	int col;	// number of collisions 
	double crmax; // Maximum relative velocity between particles

 //////////////////// MEMBER FUNCTIONS /////////////////

    public:

    particle(){}
    void init(odtline *_odtl); /// in main.cc after eddl.init()
    ~particle();

	void setVelocity();        
    void writeData(string fname, const double time); /// in inputoutput::dumpLineIfNeeded 
	void advanceParticleAndCalcWeights(const vector<double> &dxc, double dt); /// particle advancement (call in solver::diffusionCatchUpIfNeeded after advanceOdt())
	double momRhsSource(const vector<double> &dxc, double dt, double uG, double rhsMix, double rhsSrc); /// computes penalties for mom equation

	void oneWayCoupl(double time);//  governes one-way coupling (call in solver::calculateSolution) 
	bool twoWayCoupl(double time, int iStart, int iEnd); //  governes two-way coupling (call in solver::sampleEddyAndImplementIfAccepted)
 	bool computeEddyVelocities(); /// computes eddy velocities for particle-eddy interaction
	bool timePEIandNewPartProp(double invTauEddy, double time); // computes particle-eddy interaction time and new position and velocities

	bool eddyTauPartSrc(int iStart, int iEnd, const double Z_value); // recomputes eddy time-scale due to particle-energy transfer
    bool set_kernel_coefficients();
	void applyPartProp(); // applies new particle properties after particle-eddy interaction
	
	void initStats();
	void statsTime(const double tStep);
	void statsSetOld();
	void statsChange(const int &jj);
	void statsOutput();

	void collision(const vector<double> &dxc,  double dt); // governes particle collison
	void collisionData(int j); // collision data output

	//////////////help functions //////////////////
	void checkBCandSetLinePos(int i);
	void partLinePosToIndex(int i);
    void getDisplacement(vector<double> y0, vector<double> yNew);
	void getGasVelocity(int k, int cCoord);
	double eddyCrossingTime(double velE, double velG, double velP, double fDivTauP, double g, double partPos, double edgePos, double T1, double T2);
	double pPosRelToEdge(double Evel, double Gvel, double Pvel, double fDivTauP, double AG, double pPos0, double edgePos0, double T);
	double ddt_pPosRelToEdge(double Evel, double Gvel, double Pvel, double fDivTauP, double AG, double T);
	void set_fDivTauP(int i, double uvelP, double vvelP, double wvelP, double uvel, double vvel, double wvel);
};
////////////////////////////////////////////////////////////////////////////////

#endif
