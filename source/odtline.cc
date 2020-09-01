
#include "odtline.h"
#include "processor.h"
#include "lv.h"
#include "lv_uvw.h"
#include "lv_pos.h"
#include "lv_posf.h"
#include "lv_rho.h"
#include "lv_dvisc.h"
#include "lv_sca.h"
#include "lv_aDL.h"
#include "odtcase_channel.h"
#include "odtcase_channelScalar.h"
#include "odtcase_jetMixlRxn.h"
#include "odtcase_jetFlame.h"
#include "odtcase_coldPropaneJet.h"
#include "odtcase_coldJet.h"
#include "odtcase_flmlt.h"
#include "odtcase_hips.h"
#include "odtcase_hips_comb.h"
#include "odtcase_hips_simpleRxn.h"
#include "odtcase_RT.h"
#include "odtcase_shearFlow.h"
#include "odtcase_HIT.h"
#include <cmath>
#include <iomanip>

extern processor proc;

/////////////////////////////////////////////////////////////////////
/** Constructor
 */

odtline::odtline(odtline *p_odtl,
                  odtparam *p_odtp) {

     odtl = p_odtl;
     odtp = p_odtp;

 }

/////////////////////////////////////////////////////////////////////
/** Initializer
 */

void odtline::init(inputoutput *p_io,
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
                   bool        LisEddyLine) {

    //----------------------

    odtc = 0;               // initialize for destruction of eddy lines

    io     = p_io;
    mesher = p_mesher;
    gas    = p_gas;
    tran   = p_tran;
    strm   = p_strm;
    mimx   = p_mimx;
    ed     = p_ed;
    eddl   = p_eddl;
    solv   = p_solv;
    part   = p_part;
    kern   = p_kern;
    rand   = p_rand;

    //----------------------

    ngrd    = odtp->ngrd0;
    ngrdf   = ngrd + 1;

    //----------------------

    if(LisEddyLine) {        // eddy line needs less data
        initEddyLine();
        return;
    }

    //----------------------
    io->init(this);
    odtp->init(this);
    ed->init(this, eddl);
    solv->init(this);
    // mesher is init below in caseinit for phi
    // strm is init below in caseinit  (odtc), (if needed)
    // mimx is init below since it needs v[] set for cvode

    //---------------------- Continue setting up the case using the case_somecase class.
    // Adds to the above variable list, and initializes solution for the run

     if(odtp->probType == "CHANNEL")
         odtc = new odtcase_channel();    // cold channel flow

     else if(odtp->probType == "CHANNEL_SCALAR")
         odtc = new odtcase_channelScalar();  // cold channel flow with passive scalar

     else if(odtp->probType == "JETMIXL_RXN")
         odtc = new odtcase_jetMixlRxn(); // jet, wake, mixing layer with gaseous reaction

     else if(odtp->probType == "COLDPROPANEJET")
         odtc = new odtcase_coldPropaneJet(); // TNF jet

     else if(odtp->probType == "COLDJET")
         odtc = new odtcase_coldJet(); // Hussein 1994

     else if(odtp->probType == "JETFLAME")
         odtc = new odtcase_jetFlame(); // Shaddix jet

     else if(odtp->probType == "FLMLT")
         odtc = new odtcase_flmlt(); // flamelet

     else if(odtp->probType == "HIPS")
         odtc = new odtcase_hips(); // hips

     else if(odtp->probType == "HIPS_COMB")
         odtc = new odtcase_hips_comb(); // hips_comb

     else if(odtp->probType == "HIPS_SIMPLERXN")
         odtc = new odtcase_hips_simpleRxn(); // hips_simpleRxn

     else if(odtp->probType == "RT")
         odtc = new odtcase_RT();      // simple Rayleigh Taylor flow

     else if(odtp->probType == "SHEARFLOW") // homogeneous shear turbulence
         odtc = new odtcase_shearFlow();

     else if(odtp->probType == "HIT")  // homogeneous isotropic turbulence
         odtc = new odtcase_HIT();
     
     else {
         cout << endl << "ERROR, probType UNKNOWN" << endl;
         exit(0);
     }

    if(odtp->kernEv)
        kern->init(this);

    odtc->init(this);

    //----------------------

    for(int k=0; k<v.size(); k++)
        varMap[v.at(k)->var_name] = v.at(k);

    nTrans = 0;
    for(int k=0; k<v.size(); k++)
        if(v.at(k)->L_transported)
            nTrans++;

    //----------------------

    mimx->init(this);

    //----------------------

    if(odtp->Lrestart) {
        io->loadVarsFromRestartFile();
        io->set_iNextDumpTime(odtp->trst);
    }

    //------MF--------------
    part->init(this);
    //----------------------

}

/////////////////////////////////////////////////////////////////////
/** Destructor
 */

odtline::~odtline() {
    for(int k=0; k<v.size(); k++)
        delete v.at(k);
    if(odtc) delete odtc;
}

/////////////////////////////////////////////////////////////////////
/** Compute size of domain based on faces.
 */

double odtline::Ldomain() {
     return posf->d.at(ngrd) - posf->d.at(0);
}

/////////////////////////////////////////////////////////////////////
/** Initialize data members of the eddy line.
 *  Note, none of the other members of this ODTline should be used (like random).
 *  Note, all variables here should have corresponding variables (by var_name) in the
 *     main odtl. This is needed for using the eddl object.
 */

void odtline::initEddyLine() {

    v.push_back(new lv_pos(  this, "pos",   false, true));
    v.push_back(new lv_posf( this, "posf",  false, true));
    v.push_back(new lv_uvw(  this, "uvel",  true,  true));   // last are: L_transported, L_output
    v.push_back(new lv_uvw(  this, "vvel",  true,  true));
    v.push_back(new lv_uvw(  this, "wvel",  true,  true));
    v.push_back(new lv_rho(  this, "rho",   false, false));
    v.push_back(new lv_dvisc(this, "dvisc", false, false));
    if(odtl->odtp->LdoDL)
       v.push_back(new lv_aDL(this, "aDL",   false, false));

    int k = 0;
    pos   = v.at(k++);
    posf  = v.at(k++);
    uvel  = v.at(k++);
    vvel  = v.at(k++);
    wvel  = v.at(k++);
    rho   = v.at(k++);
    dvisc = v.at(k++);
    if(odtl->odtp->LdoDL)
        aDL   = v.at(k++);

}

/////////////////////////////////////////////////////////////////////
/** Set the line from a region of the odtl.  Normally called by eddy line.
 *  @param i1 \input index of starting cell of odtl to build from
 *  @param i2 \input index of ending cell of odtl to build from
 *  If i2 < i1, we have a periodic region (wrap around the domain).
 *     This only happens in planar cases, not cylindrical or sphericial.
 *  nonwrap: |   | * | * | * | * | * | * |   |   |
 *                i1                  i2
 *  new line consists of *'d cells
 *
 *  Wrap: | 4 | 5 |   |   |   |   | 1 | 2 | 3 |
 *             i2                  i1
 *  New line consists of #'d cells: 1 2 3 4 5}
 */

void odtline::setLineFromRegion(const int i1, const int i2) {

	if(!odtl->ed->LperiodicEddy){
		ngrd  = i2-i1+1;
		ngrdf = ngrd+1;
		for(int k=0; k<v.size(); k++)
        	v.at(k)->setLvFromRegion(i1,i2);
	}
	else{
	ngrd = odtl->ngrd-i1+i2+1;
    	ngrdf = ngrd+1;

    	for(int k=0; k<v.size(); k++){
        	v.at(k)->setLvFromRegion(i1,i2);
			if(v.at(k)->var_name=="uvel" && odtl->odtp->probType=="SHEARFLOW"){
				double vmax          = odtl->io->initParams["Srate"].as<double>()*odtl->io->params["domainLength"].as<double>();
				for(int i = odtl->ngrd-i1; i < ngrd; i++)
					v.at(k)->d.at(i) += vmax;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
/** Find index of cell for given position (residing in cell).
 *  Start search assuming a uniform grid,
 *  then search forward or back till hit the cell index.
 *  If position is on cell face j, then if LowSide true, return j, else j-1.         \n
 * For start of eddy region, set LowSide to true                                     \n
 * For end of eddy region, set LowSide to false                                      \n
 * (This is so triplet maps don't overlap cells)
 *                                                                                   <pre><code>
 * e.g., usual:   | { | | | | } |    5 pts, eddy pos between cell faces
 *       okay:    {   | | | |   }    5 pts, eddy pos on cell faces (1 or both)
 *       bad:     |   { | | }   |    5 pts, eddy pos on internal faces (1 or both)
 *                                                                                   </code></pre>
 * @param position \input position to find corresponding index.
 * @param LowSide  \input flag true, then return j if position is on cell face j, else j-1.
 * @return index of position.
 */

int odtline::linePositionToIndex(double position, const bool LowSide, int dbg) {

    if(abs(position-posf->d.at(0)) < 1.0E-14)
        return 0;
    if(abs(position-posf->d.at(ngrd)) < 1.0E-14)
        return ngrd-1;

    if(position < posf->d.at(0))         // for periodic (from eddies only)
        position += Ldomain();
    if(position > posf->d.at(ngrd))
        position -= Ldomain();

    if(position < posf->d.at(0) || position > posf->d.at(ngrd)) {
       *io->ostrm << "\ndbg = " << dbg << endl; //doldb
       *io->ostrm << scientific;
       *io->ostrm << setprecision(14);
       *io->ostrm << "\n ERROR odt_grid::linePositionToIndex position < posf->d.at(0) or > posf->d.at(ngrd) \n"
               " and at processor's id---> " << proc.myid
               <<" Value of position is---> "<<position << " and values of posf->d.at(0) and posf->d.at(ngrd) are "
               <<posf->d.at(0)<< " and "<<posf->d.at(ngrd) <<" respectively "<< endl;
       //io->outputProperties("dbg.dat", 0.0); //doldb
       exit(0);
    }

    int i;
    int ipos = static_cast<int>((position-posf->d.at(0))/Ldomain()*ngrd);

    if(posf->d.at(ipos+1) > position) {      // case 1: grd skewed more pts on right half
        for(i=ipos+1; i>=0; i--)  {
            if(posf->d.at(i) <= position) {
                if(position == posf->d.at(i)) {
                    if(LowSide)
                        return i;
                    else
                        return i-1;
                }
                else
                    return i;
            }
        }
    }

    else  {                           // case 2: grd skewed more pts on left half
        for(i=ipos+1; i<=ngrdf; i++) {
            if(posf->d.at(i) >= position) {
                if(position == posf->d.at(i)) {
                    if(LowSide)
                        return i;
                    else
                        return i-1;
                }
                else
                    return i-1;
            }
        }
    }

    *io->ostrm << "\n\n******** ERROR IN odt_grid::linePositionToIndex "
         << position << '\t' << posf->d.at(0) << '\t' << posf->d.at(ngrd) << '\t' << endl << endl;

    return -1;
}

/////////////////////////////////////////////////////////////////////
/** Cycle line for periodic flows.
 *  @param icycle \input move all cells before and including this one
 *   to the end of the line.
 *  @return the cycle distance (used for backcycling).
 */

double odtline::cyclePeriodicLine(const int icycle) {

    double cycleDistance = posf->d.at(icycle+1)-posf->d.at(0);
    double posFraction = (odtp->domainLength+pos->d.at(0)-pos->d.at(ngrd-2))/(pos->d.at(ngrd-1)-pos->d.at(ngrd-2));
    if(odtp->probType=="SHEARFLOW")  
        pJump = uvel->d.at(ngrd-1) - uvel->d.at(0);
        //pJump = io->initParams["Srate"].as<double>()*io->params["domainLength"].as<double>();     
        //pJump = uvel->d.at(ngrd-2)+posFraction*(uvel->d.at(ngrd-1)-uvel->d.at(ngrd-2))-uvel->d.at(0);     
    for(int k=0; k<v.size(); k++) {
        if (v.at(k)->var_name=="pos" || v.at(k)->var_name=="posf")
            continue;
        v.at(k)->d.insert(v.at(k)->d.end(),   v.at(k)->d.begin(), v.at(k)->d.begin()+icycle+1);
        v.at(k)->d.erase( v.at(k)->d.begin(), v.at(k)->d.begin()+icycle+1);
		if(v.at(k)->var_name=="uvel" && odtp->probType=="SHEARFLOW"){
        	for(int i = ngrd-icycle-1; i < ngrd; i++)
            	v.at(k)->d.at(i) += pJump;
        }
    }

    //---------- now do posf, and pos
    double xend = posf->d.at(ngrd);
    for(int i=1; i<=icycle+1; i++)
        posf->d.push_back(xend+(posf->d.at(i)-posf->d.at(0)));
    posf->d.erase(posf->d.begin(), posf->d.begin()+icycle+1);

    pos->setVar();     // does a little extra work (whole line) but doesn't happen that often
                       //    only when periodic eddies are accepted.
	
    return cycleDistance;
}

/////////////////////////////////////////////////////////////////////
/** Back cycle line for periodic flows. Intended to be called some time
 *  after cyclePeriodicLine is called.
 *  Splits the cell at posf.at(ngrd) - backCycleDistace, then moves end cells
 *     after the split to the beginning of the line.
 *  @param \input distance from the end to split and move the line.
 */

void odtline::backCyclePeriodicLine(const double backCycleDistance) {

    double xend = posf->d.at(ngrd) - backCycleDistance;     // end loc.
    double icycle = linePositionToIndex(xend, true, 1);  // cycle cells greater than this to beginning
	double Ldomain = io->params["domainLength"].as<double>();

    //------------ split the cell where the back cycle happens

    vector<double> interPos(3);
    if(abs(posf->d.at(icycle) - xend) > 1.0e-15) {
        interPos.at(0) = posf->d.at(icycle);
        interPos.at(1) = xend;
        interPos.at(2) = posf->d.at(icycle+1);
        mesher->splitCell(icycle, 1, interPos);
        icycle++;
    }

    //------------ now move the cells

    int nmove = ngrd-icycle;

    for(int k=0; k<v.size(); k++) {
        if (v.at(k)->var_name=="pos" || v.at(k)->var_name=="posf")
            continue;
		for(int i = 0; i<nmove; i++)
			v.at(k)->d.insert(v.at(k)->d.begin(), v.at(k)->d.at(v.at(k)->d.size()-i-1));
       	v.at(k)->d.erase( v.at(k)->d.begin()+icycle+nmove, v.at(k)->d.end() );
		if(v.at(k)->var_name=="uvel" && odtp->probType == "SHEARFLOW"){
            for(int i = 0; i < ngrd-icycle; i++){
         		v.at(k)->d.at(i) -= pJump;
			}
        }
    }

    //---------- now do posf, and pos

    double xstart_orig = posf->d.at(0);

    posf->d.insert(posf->d.begin(), nmove, 0.0);
    icycle += nmove;
    for(int i=0; i<nmove; i++)
        posf->d.at(nmove-1-i) = xstart_orig - (posf->d.at(posf->d.size()-1) - posf->d.at(posf->d.size()-2-i));

    posf->d.erase(posf->d.begin()+icycle+1, posf->d.end());

    pos->setVar();     // does a little extra work (whole line) but doesn't happen that often
                       //    only when periodic eddies are accepted.
	
	ngrd = pos->d.size();
	ngrdf = ngrd+1; 
}

