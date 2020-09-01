/**
 * @file lv.cc
 * Header file for class lv
 */


#include "lv.h"
#include "odtline.h"


////////////////////////////////////////////////////////////////////////////////
/*! lv constructor function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

lv::lv(odtline    *line,
       const      string s,
       const bool Lt,
       const bool Lo) {

    odtl          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(odtl->ngrd, 0.0);

    LagSrc = false;

}

////////////////////////////////////////////////////////////////////////////////
/*! lv splitCell function
 *
 * @param isplt  \input index of cell to split
 * @param nsplt  \input number of cells to split cell into
 * @param cellFaces \input original left edge, new interior faces, orig. right edge.
 */

void lv::splitCell(const int isplt,
                   const int nsplt,
                   const vector<double> &cellFaces) {

    d.insert( d.begin() + isplt, nsplt, d.at(isplt) );

}

////////////////////////////////////////////////////////////////////////////////
/*! lv merger2cells function
 *
 * Function presumes that the variable being merged is a quantity per unit mass.
 * Merging conservatively: (rho*phi*dx)_merged = (rho*phi*dx)_1 + (rho*phi*dx)_2
 * Then solve for phi_merged.
 *
 * @param imrg \input merge cells imrg and imrg+1
 * @param m1   \input mass in cell imrg
 * @param m2   \input mass in cell imrg
 * @param LconstVolume \input (for posf, default is false)
 */

void lv::merge2cells(const int    imrg,
                     const double m1,
                     const double m2,
                     const bool   LconstVolume) {

    d.at(imrg) = (d.at(imrg)*m1 + d.at(imrg+1)*m2 ) / (m1+m2);

    d.erase(d.begin() + imrg+1);

}

////////////////////////////////////////////////////////////////////////////////
/*! interpolate a cell centered variable to a face
 *  by harmonic interpolation which gives roughly an upwind flux
 */
void lv::interpVarToFacesHarmonic(const vector<double> &cvar, vector<double> &fvar){

    // todo: fill this in

    double dx1, dx2, k1, k2;
    int i, im;

    double dfirst;      // store the first face diff till end
    double dlast;       // store the last face diff till end
    double denom;

    //------ do edges

    if (odtl->odtp->Lperiodic) {
        i = 0;
        im = odtl->ngrd - 1;

        dx1 = odtl->posf->d.at(odtl->ngrd) - odtl->pos->d.at(im);
        dx2 = odtl->pos->d.at(i) - odtl->posf->d.at(i);
        k1 = cvar.at(im);
        k2 = cvar.at(i);
        denom = k1*dx2+k2*dx1;
        if(abs(denom)==0.0)
            dfirst = 0.0;
        else
            dfirst = k1 * k2 * (dx1 + dx2) / denom; // first face
        dlast = dfirst; // last face
    }
    else {
        dfirst = cvar.at(0); // first face
        dlast = cvar.at(odtl->ngrd-1); // last face
    }

    //------ do interior faces
    // we assume pos[i] is located right in the middle of posf[i] and posf[i+1]
    // so we take always the left half of cells

    for (i=1, im=0; i < odtl->ngrd; i++, im++) {
        dx1 = odtl->pos->d.at(im) - odtl->posf->d.at(im);
        dx2 = odtl->pos->d.at(i) - odtl->posf->d.at(i);
        k1 = cvar.at(im);
        k2 = cvar.at(i);
        denom = k1*dx2+k2*dx1;
        if(abs(denom)==0.0)
            fvar.at(i) = 0.0;
        else
            fvar.at(i) = k1 * k2 * (dx1 + dx2) / denom;
    }

    fvar.at(0)          = dfirst; // insert the first face flux
    fvar.at(odtl->ngrd) = dlast;  // insert the last face flux

}

///////////////////////////////////////////////////////////////////////////////
/** Interpolate variables to single face
 *
 * @param iface \input face index to interpolate to.
 * @param vec   \input variable being interpolated.
 * @return return interpolated variable at desired face.
 */
double lv::linearInterpToFace(const int &iface, const vector<double> &vec) {

    double x1, x2, y1, y2;

    if (iface == 0) {
        if (odtl->odtp->Lperiodic) {
            x1 = odtl->pos->d.at(odtl->ngrd - 1) - odtl->Ldomain();
            x2 = odtl->pos->d.at(0);
            y1 = vec.at(odtl->ngrd - 1);
            y2 = vec.at(0);
        } else {
            return vec.at(0);
        }
    } else if (iface == odtl->ngrd) {
        if (odtl->odtp->Lperiodic) {
            x1 = odtl->pos->d.at(odtl->ngrd - 1);
            x2 = odtl->pos->d.at(0) + odtl->Ldomain();
            y1 = vec.at(odtl->ngrd - 1);
            y2 = vec.at(0);
        } else {
            return vec.at(odtl->ngrd - 1);
        }
    } else {
        x1 = odtl->pos->d.at(iface - 1);
        x2 = odtl->pos->d.at(iface);
        y1 = vec.at(iface - 1);
        y2 = vec.at(iface);
    }

    return y1 + (y2 - y1) / (x2 - x1) * (odtl->posf->d.at(iface) - x1);
}

////////////////////////////////////////////////////////////////////////////////
/*! Set data array from region of odtline.
 *  @param i1 \input index of starting cell of odtl to build from
 *  @param i2 \input index of ending cell of odtl to build from
 *  See odtline::setLineFromRegion for additional details.
 */

void lv::setLvFromRegion(const int i1, const int i2){

    // note, we are owned by the eddyline, so odtl is eddl, so to get odtl data, use odtl->odtl
    const vector<double> &odtl_data = odtl->odtl->varMap.find(var_name)->second->d;

    if(i2 >= i1)
        d.assign(odtl_data.begin()+i1, odtl_data.begin()+i2+1  );
    else {           // wrap around (periodic assignment)
        d.assign(odtl_data.begin()+i1, odtl_data.end());
        d.insert(d.end(), odtl_data.begin(), odtl_data.begin()+i2+1 );
    }

}

////////////////////////////////////////////////////////////////////////////////
/*! Resize data
 */

void lv::resize() {
    d.resize(odtl->ngrd);
}
