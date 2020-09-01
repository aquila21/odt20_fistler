/**
 * @file lv_pos.cc
 * Header file for class lv_pos
 */


#include "lv_pos.h"
#include "odtline.h"
#include <cstdlib>

////////////////////////////////////////////////////////////////////////////////
/*! lv_pos  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_pos::lv_pos(odtline    *line,
               const      string s,
               const bool Lt,
               const bool Lo) {

    odtl          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(odtl->ngrd, 0.0);

    double dx = odtl->odtp->domainLength / odtl->ngrd;
    d.at(0) = odtl->odtp->xDomainCenter + 0.5*(-odtl->odtp->domainLength + dx);
    for(int i=1; i<odtl->ngrd; i++)
        d.at(i) = d.at(i-1) + dx;

}

////////////////////////////////////////////////////////////////////////////////
/*! lv_pos splitCell function
 *
 * @param isplt  \input index of cell to split
 * @param nsplt  \input number of cells to split cell into
 * @param cellFaces \input original left edge, new interior faces, orig. right edge.
 *
 * note, the number of cells is changed in the calling function, not here
 * note, this is organized to be independent of whether we split posf before or after.
 * (that is, we often compute pos1 as simply 0.5*(posf1+posf2), but not here.
 */

void lv_pos::splitCell(const int isplt,
                       const int nsplt,
                       const vector<double> &cellFaces) {

    d.insert( d.begin() +isplt+1, nsplt, 0.0);
    for(int i=isplt, j=0; i<=isplt+nsplt; i++,j++)
        d.at(i) = 0.5*(cellFaces.at(j)+cellFaces.at(j+1));

}

////////////////////////////////////////////////////////////////////////////////
/*! lv_pos merge2cells function
 *
 * Function presumes that the variable being merged is a quantity per unit mass.
 * Merging conservatively: (rho*phi*dx)_merged = (rho*phi*dx)_1 + (rho*phi*dx)_2
 * Then solve for phi_merged.
 *
 * NOTE: this should only be called when posf is correct.
 *
 * @param imrg    \input merge cells imrg and imrg+1
 * @param imrg \input merge cells imrg and imrg+1
 * @param m1   \input mass in cell imrg
 * @param m2   \input mass in cell imrg
 * @param LconstVolume \input (for posf, default is false)
 */

void lv_pos::merge2cells(const int    imrg,
                         const double m1,
                         const double m2,
                         const bool   LconstVolume) {

    if(LconstVolume || odtl->odtp->bcType=="WALL") {
        d.erase(d.begin() + imrg+1);
        d.at(imrg) = 0.5*(odtl->posf->d.at(imrg) + odtl->posf->d.at(imrg+1));
    }
    else
        setVar();
}

////////////////////////////////////////////////////////////////////////////////
/*! lv_pos setVar function
 *  @param ipt \input optional point to compute at
 *  Sets the grid position from posf
 *  So, you should make sure posf is consistent with pos before calling.
 *
 *  NOTE: this should only be called when posf is correct.
 */

void lv_pos::setVar(const int ipt){

    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be = -1" << endl;
        exit(0);
    }

    d.resize(odtl->posf->d.size()-1);

    for(int i=0; i<d.size(); i++)
        d.at(i) = 0.5*(odtl->posf->d.at(i) + odtl->posf->d.at(i+1));

}

////////////////////////////////////////////////////////////////////////////////
/*! Set data array from region of odtline.
 *  @param i1 \input index of starting cell of odtl to build from
 *  @param i2 \input index of ending cell of odtl to build from
 *  See odtline::setLineFromRegion for additional details.
 */

void lv_pos::setLvFromRegion(const int i1, const int i2){

    // note, we are owned by the eddyline, so odtl is eddl, so to get odtl data, use odtl->odtl
    const vector<double> &odtl_data = odtl->odtl->varMap.find(var_name)->second->d;

    if(i2 >= i1)
        d.assign(odtl_data.begin()+i1, odtl_data.begin()+i2+1  );
    else {           // wrap around (periodic assignment)
        d.assign(odtl_data.begin()+i1, odtl_data.end());
        int idmb = d.size();
        d.insert(d.end(), odtl_data.begin(), odtl_data.begin()+i2+1 );
        for(int i=idmb; i<d.size(); i++)
            d.at(i)+=odtl->odtl->Ldomain();
    }

}
