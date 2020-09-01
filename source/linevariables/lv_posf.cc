/**
 * @file lv_posf.cc
 * Header file for class lv_posf
 */


#include "lv_posf.h"
#include "odtline.h"
#include <iostream>
#include <cstdlib>
#include <numeric> //accumulate

////////////////////////////////////////////////////////////////////////////////
/*! lv_posf  constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv_posf::lv_posf(odtline    *line,
                 const      string s,
                 const bool Lt,
                 const bool Lo) {

    odtl          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(odtl->ngrdf, 0.0);

    double dx = odtl->odtp->domainLength / odtl->ngrd;
    d.at(0) = odtl->odtp->xDomainCenter - 0.5*odtl->odtp->domainLength;
    for(int i=1; i<odtl->ngrdf; i++)
        d.at(i) = d.at(i-1) + dx;

}

////////////////////////////////////////////////////////////////////////////////
/*! lv_posf splitCell function
 *
 * @param isplt  \input index of cell to split
 * @param nsplt  \input number of cells to split cell into
 * @param cellFaces \input original left edge, new interior faces, orig. right edge.
 *
 * note, the number of cells is changed in the calling function, not here
 */

void lv_posf::splitCell(const int isplt,
                        const int nsplt,
                        const vector<double> &cellFaces) {

    d.insert( d.begin() + isplt+1, nsplt, 0.0 );
    for(int i=isplt+1, j=1; i<=isplt+nsplt; i++, j++)
        d.at(i) = cellFaces.at(j);

}
////////////////////////////////////////////////////////////////////////////////
/*! lv_posf merger2cells function
 *
 * Function presumes that the variable being merged is a quantity per unit mass.
 * Merging conservatively: (rho*phi*dx)_merged = (rho*phi*dx)_1 + (rho*phi*dx)_2
 * Then solve for phi_merged.
 *
 * @param imrg \input merge cells imrg and imrg+1
 * @param m1   \input mass in cell imrg
 * @param m2   \input mass in cell imrg
 * @param LconstVolume \input Thermo says an adiabatic const P mixing will change volume,
 *            but sometimes we want to retain the old volume (e.g. when we merge a small cell
 *            on the edge of the domain when enforcing the boundaries. In that case, we
 *            will chop or extend the cell anyway, so there is no conservation issue.
 */

void lv_posf::merge2cells(const int    imrg,
                          const double m1,
                          const double m2,
                          const bool   LconstVolume) {

    d.erase(d.begin() + imrg+1);

    if(LconstVolume || odtl->odtp->bcType=="WALL") {   //todo: generalize this (works for constant density flows only, and not spatial (due to velocity)).
        odtl->pos->merge2cells(imrg, m1, m2, LconstVolume);
        return;
    }

    double invC = 1.0/odtl->odtp->cCoord;
    double C    = odtl->odtp->cCoord;

    vector<double> dxc;
    odtl->mesher->setGridDxc(odtl, dxc, odtl->odtp->cCoord);
    dxc.at(imrg) = (m1+m2)/odtl->rho->d.at(imrg);
    if(odtl->odtp->Lspatial)
        dxc.at(imrg) /= odtl->uvel->d.at(imrg);


    odtl->mesher->setGridFromDxc(dxc);       // does pos also --> pos is done for each merge and also at the end in meshManager::merge2cells

    ////----------- outflow boundary on both sides.

    //if(odtl->odtp->bcType=="OUTFLOW") {
    //
    //    double V2tot = accumulate(dxc.begin(), dxc.end(), 0.0);
    //    double dmb;
    //    d[0] = -pow(0.5*V2tot, invC);
    //    for(int ie=1, iw=0, i=0; ie<d.size(); ie++, iw++, i++) {
    //        if(d[iw] <= 0.0) {
    //            dmb = pow(abs(d[iw]),C) - dxc[i];
    //            if(dmb >=0) d[ie] = -pow( dmb, invC);
    //            else        d[ie] =  pow(-dmb, invC);
    //        }
    //        else {
    //            d[ie] = pow(pow(d[iw],C) + dxc[i], invC);
    //        }
    //    }
    //}

    ////------------------ Wall on left, outlet on right. Assume all posf values are >= 0.0, expand to the right.

    //else if(odtl->odtp->bcType == "Wall_OUT") {
    //    d[0] = 0.0;
    //    for(int ie=1, iw=0, i=0; ie<odtl->ngrdf; ie++, iw++, i++)
    //        d[ie] = pow(pow(d[iw],C) + dxc[i], invC);
    //}

    ////------------------

    //else {
    //    cout << endl << "ERROR: lv_posf::merge2cells: not setup for given bcType" << endl;
    //    exit(0);
    //}

}

////////////////////////////////////////////////////////////////////////////////
/*! lv_posf merger2cells function
 *
 * Function presumes that the variable being merged is a quantity per unit mass.
 * Merging conservatively: (rho*phi*dx)_merged = (rho*phi*dx)_1 + (rho*phi*dx)_2
 * Then solve for phi_merged.
 *
 * @param imrg \input merge cells imrg and imrg+1
 * @param m1   \input mass in cell imrg
 * @param m2   \input mass in cell imrg
 * @param LconstVolume \input Thermo says an adiabatic const P mixing will change volume,
 *            but sometimes we want to retain the old volume (e.g. when we merge a small cell
 *            on the edge of the domain when enforcing the boundaries. In that case, we
 *            will chop or extend the cell anyway, so there is no conservation issue.

void lv_posf::merge2cells(const int    imrg,
                          const double m1,
                          const double m2,
                          const bool   LconstVolume) {

    d.erase(d.begin() + imrg+1);

    if(LconstVolume || odtl->odtp->bcType=="WALL")    //todo: generalize this (works for constant density flows only, and not spatial (due to velocity)).
        return;

    vector<double> dxc;
    odtl->mesher->setGridDxc(odtl, dxc);

    double invC = 1.0/odtl->odtp->cCoord;
    double C    = odtl->odtp->cCoord;

    //----------- outflow boundary on both sides.

    if(odtl->odtp->bcType=="OUTFLOW") {

        double dmb;
        double pm1;

        //----------- do cell imrg faces

        double xc;                   // position in cell imrg with half cell vol on each side
        if(d.at(imrg) >= 0.0) {         // all on right (before expand)
            xc = pow(pow(d.at(imrg),C)+0.5*dxc.at(imrg),invC);      // xc is found using the original cell vol.
            dxc.at(imrg) = odtl->odtp->Lspatial ? (m1+m2)/odtl->rho->d.at(imrg)/odtl->uvel->d.at(imrg)  : (m1+m2)/odtl->rho->d.at(imrg); // now we update the vol to compute faces from xc
            d.at(imrg+1) = pow(pow(xc,C)+0.5*dxc.at(imrg),invC);
            dmb       = pow(d.at(imrg+1),C)-dxc.at(imrg);
            pm1       = dmb<0.0 ? -1.0 : 1.0;   // (may cross after expand)
            d.at(imrg)   = pm1*pow(pm1*dmb,invC);
        }
        else if(d.at(imrg+1) <= 0.0) {  // all on left (before expand)
            xc = -pow(pow(abs(d.at(imrg+1)),C)+0.5*dxc.at(imrg),invC);
            dxc.at(imrg) = odtl->odtp->Lspatial ? (m1+m2)/odtl->rho->d.at(imrg)/odtl->uvel->d.at(imrg)  : (m1+m2)/odtl->rho->d.at(imrg);
            d.at(imrg)   = -pow(pow(abs(xc),C)+0.5*dxc.at(imrg),invC);
            dmb       = pow(abs(d.at(imrg)),C) - dxc.at(imrg);
            pm1       = dmb<0.0 ? -1.0 : 1.0;   // (may cross after expand)
            d.at(imrg+1) = -pm1*pow(pm1*dmb,invC);

        }
        else {                       // cell splits center (before expand)
            dmb = pow(d.at(imrg+1),C) - 0.5*dxc.at(imrg);
            pm1 = dmb<0.0 ? -1.0 : 1.0;
            xc  = pm1*pow(pm1*dmb,invC);
            dxc.at(imrg) = odtl->odtp->Lspatial ? (m1+m2)/odtl->rho->d.at(imrg)/odtl->uvel->d.at(imrg)  : (m1+m2)/odtl->rho->d.at(imrg);
            if(xc >= 0.0) {
                d.at(imrg+1) = pow(pow(xc,C)+0.5*dxc.at(imrg),invC);
                dmb = pow(d.at(imrg+1),C) - dxc.at(imrg);
                pm1 = dmb<0.0 ? -1.0 : 1.0;
                d.at(imrg) = pm1*pow(pm1*dmb,invC);
            }
            else {
                d.at(imrg) = -pow(pow(abs(xc),C)+0.5*dxc.at(imrg),invC);
                dmb = pow(abs(d.at(imrg)),C)-dxc.at(imrg);
                pm1 = dmb<0.0 ? -1.0 : 1.0;
                d.at(imrg+1) = -pm1*pow(pm1*dmb,invC);
            }
        }

        //----------- work right: imrg+1 to end

        for(int ie=imrg+2, iw=imrg+1, i=imrg+1; ie<d.size(); ie++, iw++, i++) {
            if(d.at(iw) > 0.0)
                d.at(ie) = pow(pow(d.at(iw),C)+dxc.at(i),invC);
            else {
                dmb = pow(abs(d.at(iw)),C)-dxc.at(i);
                pm1 = dmb<0.0 ? -1.0 : 1.0;
                d.at(ie) = -pm1*pow(pm1*dmb,invC);
            }
        }

        //----------- work left: imrg-1 to 0

        for(int iw=imrg-1, i=imrg-1, ie=imrg; iw>=0; iw--, i--, ie--) {
            if(d.at(ie) < 0.0)
                d.at(iw) = -pow(pow(abs(d.at(ie)),C)+dxc.at(i),invC);
            else {
                dmb = pow(d.at(ie),C) - dxc.at(i);
                pm1 = dmb<0.0 ? -1.0 : 1.0;
                d.at(iw) = pm1*pow(pm1*dmb,invC);
            }
        }
    }

    //------------------ Wall on left, outlet on right. Assume all posf values are >= 0.0, expand to the right.

    else if(odtl->odtp->bcType == "Wall_OUT") {
        dxc.at(imrg) = (m1+m2)/odtl->rho->d.at(imrg);
        for(int iw=imrg, ie=imrg+1, i=imrg; ie<d.size(); ie++, iw++, i++)
            d.at(ie) = pow(pow(d.at(iw),C)+dxc.at(imrg),invC);
    }

    //------------------

    else {
        cout << endl << "ERROR: lv_posf::merge2cells: not setup for given bcType" << endl;
        exit(0);
    }

}
 */

////////////////////////////////////////////////////////////////////////////////
/*! Set data array from region of odtline.
 *  @param i1 \input index of starting cell of odtl to build from
 *  @param i2 \input index of ending cell of odtl to build from
 *  See odtline::setLineFromRegion for additional details.
 */

void lv_posf::setLvFromRegion(const int i1, const int i2){

    // note, we are owned by the eddyline, so odtl is eddl, so to get odtl data, use odtl->odtl
    const vector<double> &odtl_data = odtl->odtl->varMap.find(var_name)->second->d;

    if(i2 >= i1)
        d.assign(odtl_data.begin()+i1, odtl_data.begin()+i2+2  );
    else {           // wrap around (periodic assignment)
        d.assign(odtl_data.begin()+i1, odtl_data.end()-1);
        double idmb  = d.size();
        d.insert(d.end(), odtl_data.begin(), odtl_data.begin()+i2+2 );
        for(int i=idmb; i<d.size(); i++)
            d.at(i)+=odtl->odtl->Ldomain();
    }

}

////////////////////////////////////////////////////////////////////////////////
/*! Resize data
 */

void lv_posf::resize() {
    d.resize(odtl->ngrdf);
}













