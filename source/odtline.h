/**
 * @file odtline.h
 * Header file for class odtline
 */

#ifndef ODTLINE_H
#define ODTLINE_H

#include "lv.h"
#include "odtcase.h"
#include "inputoutput.h"
#include "odtparam.h"
#include "streams.h"
#include "micromixer.h"
#include "eddy.h"
#include "meshManager.h"
#include "solver.h"
#include "particle.h"
#include "kernel.h"
#include "randomGenerator.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"
#include <vector>
#include <string>
#include <map>

using namespace std;
using namespace Cantera;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing odtline object
 *
 *  @author David O. Lignell
 */

class odtline {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        odtline                 *odtl;     ///< (for one line to point to another (eddl))

        int                     ngrd;      ///< number of grid cells
        int                     ngrdf;     ///< number of grid cell faces = ngrd+1

        vector<lv*>             v;         ///< All line variables are stored in here.

        lv*                     pos;       ///< pointers to gas properties
        lv*                     posf;      ///< access as: posf->d[i], or posf->var_name, etc.
        lv*                     rho;
        lv*                     dvisc;
        lv*                     uvel;
        lv*                     vvel;
        lv*                     wvel;
        lv*                     sdiff;
        lv*                     sca;
        lv*                     phase;
        lv*                     enth;
        lv*                     temp;
        lv*                     mixf;
        lv*                     chi;
        lv*                     hr;
        lv*                     aDL;
        vector<lv*>::iterator   ysp;       ///< access as: ysp=v.begin(), (*ysp)->d[i] or (*(ysp+k))->d[i], or ysp[k]->d[i].
        vector<lv*>::iterator   svar;      ///< iterator for increment to go through moments (*(ysp+k))->d[i];)
        vector<lv*>::iterator   eta;       ///< iterator for increment to go through species etc. (*(ysp+k))->d[i];)

        map<string,lv*>         varMap;

        IdealGasMix             *gas;        ///< pointer to cantera thermochemistry object (reaction rates, Cp, etc.)
        Transport               *tran;       ///< pointer to cantera transport object (viscosity, diffusivity, etc.)
        streams                 *strm;       ///< pointer to gas stream properties
        inputoutput             *io;         ///< pointer to input/output object
        odtparam                *odtp;       ///< pointer to the parameters object
        micromixer              *mimx;       ///< pointer to micromixer for diffusion, reaction, line evolution.
        eddy                    *ed;         ///< pointer to object for eddy operations
        odtline                 *eddl;       ///< pointer to eddyline object
        solver                  *solv;       ///< pointer to solver object
        meshManager             *mesher;     ///< pointer to mesh manager object
        particle		*part;       ///< pointer to particle object
        kernel                  *kern;       ///< pointer to kernel event object

        randomGenerator         *rand;

        int                     nTrans;      ///< number of transported variables on the line.

        odtcase                 *odtc;       ///< odt case class: set specific vars...
		double 					pJump;       ///< jump between two periodic boundaries

    //////////////////// MEMBER FUNCTIONS /////////////////

        int    linePositionToIndex(double position, const bool LowSide, int dbg);
        void   setLineFromRegion(const int i1, const int i2);
        double cyclePeriodicLine(const int icycle);
        void   backCyclePeriodicLine(const double backCycleDistance);
        double Ldomain();

    private:

        void initEddyLine();


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        void init(inputoutput *p_io,
                  meshManager *p_mesher,
                  streams     *p_strm,
                  IdealGasMix *p_gas,
                  Transport   *p_tran,
                  micromixer  *p_mimx,
                  eddy        *p_ed,
                  odtline     *p_eddl,
                  solver      *p_solv,
                  particle    *p_part,
		  kernel      *p_kern,
                  randomGenerator *p_rand,
                  bool        LisEddyLine=false);
        odtline(odtline *p_odtl,
                odtparam *p_odtp);
        ~odtline();

};


////////////////////////////////////////////////////////////////////////////////

#endif

