/**
 * @file odtparam.h
 * Header file for class odtparam
 */

#ifndef ODTPARAM_H
#define ODTPARAM_H

#include "inputoutput.h"
#include <string>
#include <cstdlib>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing inputoutput object
 *
 *  @author David O. Lignell
 */

class odtparam {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        odtline *odtl;           ///< pointer to line object
        inputoutput *io;         ///< pointer to io object (has the input file)

        int     seed;            ///<  random number generator seed (negative to randomize it)
        double  tEnd;            ///<  ending time of realization
        double  domainLength;    ///<  length of domain (m)
        int     ngrd0;           ///<  initial grid points
        double  rho0;            ///<  initial uniform density (kg/m^3)
        double  kvisc0;          ///<  initial uniform kinematic viscosity (m^2/s)
        double  sdiff0;          ///<  initial uniform scalar diffusivity (m^2/s)
        double  dPdx;            ///<  initial pressure gradient (Pa/m)
        double  pres;            ///<  initial pressure (Pa)
        string  chemMechFile;    ///<  name of chemical mechanism file
        string  probType;        ///<  problem type: CHANNEL, CHANNEL_SCALAR, JETMIXL_RXN

        double  Z_param;         ///<  Viscous penalty parameter
        double  A_param;         ///<  Energy Distribution parameter alpha
        double  C_param;         ///<  Eddy frequency parameter
        string  LES_type;        ///<  NONE, THIRDS, ELAPSEDTIME, FRACDOMAIN, INTEGRALSCALE
        double  Z_LES;           ///<  large eddy suppression (nonpositive prevents les test)
        double  diffCFL;         ///<  multiplies min diffusion timestep
        double  cvode_atol;      ///<  absolute tolerace atol for cvode
        double  cvode_rtol;      ///<  relative tolerace rtol for cvode
        double  x0virtual;       ///<  LES virtual origin

        bool    LdoDL;           ///<  flag to do the DL energy from the DL instability
        string  radType;         ///<  radiation flag: OPTHIN, TWOFLUX
        bool    Lbuoyant;        ///<  flag to turn on bouyancy (horizontal odt line)
        bool    LPeEddy;         ///<  flag to turn on potential energy for eddies (vertical odt line)
        double  g;               ///<  gravity (default -9.81)
        string  Lsolver;         ///<  EXPLICIT, SEMI-IMPLICIT, or STRANG
        bool    Lperiodic;       ///<  periodic if true
        bool    Lspatial;        ///<  spatial formulation if true
        bool    Llem;            ///<  true if LEM
        bool    LisFlmlt;        ///<  true if solving an unsteady flamelet

        string  bcType;          ///<  OUTFLOW, PERIODIC, WALL, WALL_OUT
        int     cCoord;          ///<  1 = planar, 2 = cylindrical, 3 = spherical
        double  xDomainCenter;   ///<  position of the center of the domain

        double  gDens;           ///<  grid density for mesher
        double  dxmin;           ///<  min grid spacing: = dxmin / domain length
        double  dxmax;           ///<  max grid spacing = dxmax / domain length

        double  Pmax;            ///<  maximum eddy acceptance probability
        double  Pav;             ///<  Average acceptance probability
        double  dtfac;           ///<  maximum factor to increase dtSmean
        int     nDtSmeanWait;    ///<  number of eddy samples before increase dtSmean
        int     eddyMinCells;    ///<  eddy must overlap at least this many cells
        double  DAtimeFac;       ///<  time until catch-up adaption is DAtimeFac * dtCUmax
        double  tdfac;           ///<  factor between dtCUmax and dtCFL for temporal flows; DEFAULT = 1.0
        int     sLastDA;         ///<  size of the lastDA vector for timing adaptmesh after diff
        double  Lp;              ///<  Most probable eddy size frac of domainLength
        double  Lmax;            ///<  Max eddy size frac of domainLength
        double  Lmin;            ///<  Min eddy size frac of domainLength

        int     modDump;         ///<  accepted eddies before output file
        int     modDisp;         ///<  frequency to display results (# eddies)

        bool    LmultiPhase;     ///<  true if line has more than one "line" phase (soot or particles don't count.)
        double  eSurfTens;       ///<  surface tension, J/m2 for liquid phases

        double  uBClo;           ///<  Dirichlet velocity boundary condition.
        double  uBChi;           ///<  Dirichlet velocity boundary condition.
        double  vBClo;           ///<  Dirichlet velocity boundary condition.
        double  vBChi;           ///<  Dirichlet velocity boundary condition.
        double  wBClo;           ///<  Dirichlet velocity boundary condition.
        double  wBChi;           ///<  Dirichlet velocity boundary condition.
        double  sBClo;           ///<  Dirichlet scalar boundary condition.
        double  sBChi;           ///<  Dirichlet scalar boundary condition.
        string  hWallBCtype;     ///<  ADIABATIC or ISOTHERMAL
        double  TBClo;           ///<  Required if hWallBCtype = ISOTHERMAL
        double  TBChi;           ///<  Required if hWallBCtype = ISOTHERMAL

        bool    Lrestart;        ///<  true to restart from file, else false
        string  rstType;         ///<  "single" or "multiple"
        double  trst;            ///<  restart time (from restart file), default is 0.0;

        double  umin_spatial;    ///< min u for spatial flows; used when kernels pull velocity

        //----------------- Soot variables

        bool    Lsoot;                  ///< true for soot, false for no soot
        int     nsvar;                  ///< number of soot variables transported (# soot moments)
        double  b_coag;                 ///< coagulation rate parameter
        int     nsvar_v;                ///< number of soot variables transported (v direction)
        int     nsvar_s;                ///< number of soot variables transported (s direction)
        double  rho_soot;               ///< solid soot density
        int     Cmin;                   ///< minimum number of carbon atoms in a soot particle
        string  PSD_method;             ///< method name for soot PSD: MONO, QMOM, MOMIC
        string  nucleation_mech;        ///< soot nucleation chemistry flag
        string  growth_mech;            ///< soot growth chemistry flag
        string  oxidation_mech;         ///< soot oxidation chemistry flag
        string  coagulation_mech;       ///< soot coagulation mechanism flag

        //----------------- HIPS quantities

        bool    LisHips;         ///<  true if solving hips
        int     nLevels;         ///< number of levels in the tree: 0, 1, 2, ... N-1
        double  Afac;            ///< level lengthscale reduction factor (0.5)
        double  L0;              ///< tree lengthscale
        double  tau0;            ///< integral timescale
        double  fmix;            ///< timescale factor for micromixing

        //----------------- Forced HIT parameter
        double  kernEv;          ///< kernel event flag
	    double  Prod;            ///< TKE injection parameter
      	double  T11;             ///< large eddy turnover time
    	double  L11;             ///< Intergral length scale

        //----------------- HST parameter
      
        double p3;               ///< probability of eddy event between x and z (p1 = p2 = (1-p3)/2)  
      
	    //----------------- Particle parameter 
        int     partCoupl;       ///< flag particle-gas coupling
        double  betaP;           ///< model constant for interaction time 
        bool    typeI;           ///< instantaneous eddy interaction
        bool    typeC;           ///< continuous eddy interction
        int     Nparticle;       ///< number of particle
        double  Dparticle;       ///< diameter of sphere particle
        double  DensiParticle;   ///< density of particle
        int     typeP;           ///< type of particle (tracer(0), inertial(1), ballistic(2))
        double  Gx;              ///< gravity in x-direction
        double  Gy;              ///< gravity in y-direction
        double  Gz;              ///< gravity in z-direction 
        double  crmax;           ///< Maximum relative velocity of particles (collision)
	int 	BCparticle; 	 ///< Boundary condition for particles		
	double  alphaP;          ///< Energy distribution to neighbor cells during  TWC
	double  tparticleON;     ///< stats index when the particle phase is switched on

        int nStat;
	int ngrdStat;
    //////////////////// MEMBER FUNCTIONS /////////////////

    private:

        template <class T>
        T errMsg(const string param) {
            *io->ostrm << endl << "ERROR: missing parameter: " + param << endl;
            exit(0);
            T dummy = static_cast<T> (0);
            return dummy;
        }

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        odtparam(inputoutput *p_io);
        void init(odtline *p_odtl);
        ~odtparam(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif

