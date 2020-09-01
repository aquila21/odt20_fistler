/**
 * @file meshManager.cc
 * Header file for class meshManager
 */

#include "meshManager.h"
#include "odtline.h"
#include "odtparam.h"
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <numeric> //accumulate

using namespace std;

///////////////////////////////////////////////////////////////////////////////
/** meshManager initialization function
 *
 * @param p_odtl  \input set odtline pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void meshManager::init(odtline *p_odtl, const vector<lv*> p_phi) {

    odtl   = p_odtl;
    phi    = p_phi;

    lastDA = vector<double>(odtl->odtp->sLastDA,  0.0);

}

///////////////////////////////////////////////////////////////////////////////
/** User interface to the mesh adapter.
 *
 * @param iLowerDummy \input lower cell to adapt from (approx).
 * @param iUpperDummy \input upper cell to adapt to (approx).
 */

void meshManager::adaptGrid(int iLowerDummy, int iUpperDummy) {

    for(int i=0; i<phi.size(); i++)
        phi.at(i)->setVar();

    for (int i=0; i < odtl->ngrd; i++) {
        if(odtl->posf->d.at(i) > odtl->posf->d.at(i+1)){
            *odtl->io->ostrm << endl << "ERROR: Non-increasing odtLine!";
            //odtl->io->outputProperties("ERROR_odtline.dat", 0.0);
            exit(0);
        }
    }
    if(!odtl->odtp->LmultiPhase){
        // case is for a single phase
        if(iLowerDummy == iUpperDummy)
            return;

        if(iLowerDummy < iUpperDummy) {
            adaptGrid_details(iLowerDummy, iUpperDummy);
        }

        //---------- for periodic where adaption region wraps, adapt twice: once on each section

        else {

            posUpper = odtl->posf->d.at(iUpperDummy+1);

            adaptGrid_details(iLowerDummy, odtl->ngrd-1);

            iUpperDummy = odtl->linePositionToIndex(posUpper, false, 12);

            adaptGrid_details(0, iUpperDummy);
        }
    }
    else {
        // case is for multiple immiscable phases
        if(iLowerDummy == iUpperDummy)
            return;
        if(iLowerDummy < iUpperDummy){

            double posUpperDummy = odtl->posf->d.at(iUpperDummy+1);
            double posMiddDummy;
            double phase = odtl->phase->d.at(iLowerDummy);

            for (int iDummy = iLowerDummy+1; iDummy <= iUpperDummy; iDummy++){
                if (phase != odtl->phase->d.at(iDummy)){

                    // change to new Phase
                    phase = odtl->phase->d.at(iDummy);

                    // save position of phase change
                    posMiddDummy = odtl->posf->d.at(iDummy);

                    // adapt grid within the old phase
                    adaptGrid_details(iLowerDummy, iDummy-1);

                    // recover indexes
                    iUpperDummy = odtl->linePositionToIndex(posUpperDummy, false, 13);
                    iLowerDummy = odtl->linePositionToIndex(posMiddDummy, true, 14);
                    iDummy = iLowerDummy;
                }
            }

            // adaption of the last phase
            adaptGrid_details(iLowerDummy, iUpperDummy);
        }

        //---------- for periodic where adaption region wraps, adapt twice: once on each section
        else{
            posUpper = odtl->posf->d.at(iUpperDummy+1);

            // call adaptGrid twice, once for the lower part and once for the
            // upper part. adaptGrid is needet to take different phases into
            // account.
            adaptGrid(iLowerDummy, odtl->ngrd-1);

            iUpperDummy = odtl->linePositionToIndex(posUpper, false, 15);

            adaptGrid(0, iUpperDummy);
        }
    }

//	if(odtl->odtp->cCoord != 1) splitCellWithZero();

}


///////////////////////////////////////////////////////////////////////////////
/**
 * Adapt the grid
 *
 * @param iLowerDummy \input lower cell to adapt from (approx).
 * @param iUpperDummy \input upper cell to adapt to (approx).
 */

void meshManager::adaptGrid_details(int iLowerDummy, int iUpperDummy) {

    // set xf and yf in the region to adapt

    iLower   = iLowerDummy;
    iUpper   = iUpperDummy;
    posLower = odtl->posf->d.at(iLower);
    posUpper = odtl->posf->d.at(iUpper+1);


    xf = vector<double>(odtl->posf->d.begin()+iLower,    // face positions
            odtl->posf->d.begin()+iUpper+2);

    yf = vector<vector<double> >(phi.size(), vector<double>(xf.size(), 0.0));

    for(int iProf = 0; iProf<(int)phi.size(); iProf++) {

        if(iLower==0)
            yf.at(iProf).at(0) = phi.at(iProf)->d.at(0);
        else
            yf.at(iProf).at(0) = phi.at(iProf)->d.at(iLower-1) + (odtl->posf->d.at(iLower)-odtl->pos->d.at(iLower-1))/
                (odtl->pos->d.at(iLower)-odtl->pos->d.at(iLower-1))*
                (phi.at(iProf)->d.at(iLower)-phi.at(iProf)->d.at(iLower-1));
        if(iUpper==odtl->ngrd-1)
            yf.at(iProf).at(xf.size()-1) = phi.at(iProf)->d.at(iUpper);
        else
            yf.at(iProf).at(xf.size()-1) = phi.at(iProf)->d.at(iUpper) + (odtl->posf->d.at(iUpper+1)-odtl->pos->d.at(iUpper))/
                (odtl->pos->d.at(iUpper+1)-odtl->pos->d.at(iUpper))*
                (phi.at(iProf)->d.at(iUpper+1)-phi.at(iProf)->d.at(iUpper));
        for(int i=1,j=iLower+i; i<(int)xf.size()-1; i++, j++) {
            yf.at(iProf).at(i) = phi.at(iProf)->d.at(j-1) + (odtl->posf->d.at(j)-odtl->pos->d.at(j-1))/
                (odtl->pos->d.at(j)-odtl->pos->d.at(j-1))*
                (phi.at(iProf)->d.at(j)-phi.at(iProf)->d.at(j-1));
        }
    }
    //---------- Scale xf and yf

    //double xfMax = *max_element(xf.begin(), xf.end());
    //double xfMin = *min_element(xf.begin(), xf.end());
    double xfMax = *max_element(odtl->posf->d.begin(), odtl->posf->d.end());
    double xfMin = *min_element(odtl->posf->d.begin(), odtl->posf->d.end());
    for(int i=0; i<(int)xf.size(); i++)
        xf.at(i) = (xf.at(i)-xfMin)/(xfMax-xfMin);

    for(int iProf = 0; iProf<(int)phi.size(); iProf++) {
        //double yfMax = *max_element(yf.at(iProf).begin(), yf.at(iProf).end());
        //double yfMin = *min_element(yf.at(iProf).begin(), yf.at(iProf).end());
        double yfMax = *max_element(phi.at(iProf)->d.begin(), phi.at(iProf)->d.end());
        double yfMin = *min_element(phi.at(iProf)->d.begin(), phi.at(iProf)->d.end());
        for(int i=0; i<(int)xf.size(); i++)
            yf.at(iProf).at(i) = (yf.at(iProf).at(i)-yfMin)/(yfMax-yfMin);
    }

    //---------- Compute the new grid distance function

    vector<double> s_dist(xf.size(),0.0);                 // distance function for x,y
    vector<double> sn_dist;                               // new distance function for new grid X,Y

    double s_distCum = calcDistance(xf, yf, s_dist);

    int nNewGrd = odtl->odtp->gDens*s_distCum+1;

    sn_dist.resize(nNewGrd+1, 0.0);
    xnf.resize(nNewGrd+1, 0.0);

    double dS = s_distCum / nNewGrd;
    sn_dist.at(0) = 0.0;
    for(int i=1, im=0; i<nNewGrd+1; i++, im++)
        sn_dist.at(i) = sn_dist.at(im) + dS;

    //---------- Get the new grid using the old grid, old dist func, new dist func

    interpVec(s_dist, xf, sn_dist, xnf);                  // interp from xf --> xnf

    for(int i=0; i<(int)xnf.size(); i++)                      // Unscale xnf
        xnf.at(i) = xfMin + xnf.at(i)*(xfMax-xfMin);

    xnf.at(0) = posLower;
    xnf.at(xnf.size()-1) = posUpper;
    for(int i=0; i<(int)xnf.size(); i++)                      // doldb ?
        if(xnf.at(i) < posLower || xnf.at(i) > posUpper) {
            *odtl->io->ostrm << endl << "Error in mesher: X is out of bounds, i= " << i << endl;
            *odtl->io->ostrm << "        posLower, posUpper, xnf.at(i) ,xnf.at(0), xnf.at(end) = " <<
                posLower << ", " << posUpper << ", " << xnf.at(i)<< ", " << xnf.at(0)<< ", " << xnf.at(xnf.size()-1) << endl;
            exit(0);
        }

    for(int i=1; i<(int)xnf.size(); i++)                     // doldb ?
        if(xnf.at(i) <= xnf.at(i-1)) {
            //odtl->io->outputProperties("Line_before_error.dat");
            *odtl->io->ostrm << endl << endl;
            *odtl->io->ostrm << "i = " << i << "  xnf.size = " << xnf.size() << "  xnf = "
                 << xnf.at(i-1) <<" "<< xnf.at(i) << endl;
            *odtl->io->ostrm << endl << "Error in mesher: xnf is nonincreasing " << endl;
            exit(0);
        }

    //--------- Add in the rest of the old grid so that we can adapt desired region of WHOLE line

    if(iLower > 0)
        xnf.insert(xnf.begin(), odtl->posf->d.begin(), odtl->posf->d.begin()+iLower);
    if(iUpper < odtl->ngrd-1)
        xnf.insert(xnf.end(),   odtl->posf->d.begin()+iUpper+2, odtl->posf->d.end());

    /////////////////////////////////////////////////////////////////////////////////////////////

    //---------- Adapt the new grid

    ngrdf = xnf.size();
    ngrd  = ngrdf-1;
    setDxArray();
    set_iLower_iUpper();

    //----------- split large cells

    for (int i=0; i<ngrd; i++) {
        if (dx.at(i) > odtl->odtp->dxmax){
            dx.at(i) *= 0.5;
            dx.insert(dx.begin()+i,dx.at(i));
            ngrd++;
            ngrdf++;
            i--;
        }
    }
    xnf.resize(ngrdf);
    xnf.at(0) = odtl->posf->d.at(0);               // recover the grid
    for(int i=1; i<ngrdf; i++)
        xnf.at(i) = xnf.at(i-1)+dx.at(i-1);

    //----------- Small Cells

    if(odtl->odtp->LmultiPhase)
        mergeSmallCellsMP();      // Function does the same as the old one, but has a phase check
    else
        mergeSmallCells();      // Function does the same as the old one, but has a phase check

    xnf.resize(ngrdf);
    xnf.at(0) = odtl->posf->d.at(0);               // recover the grid
    for(int i=1; i<ngrdf; i++)
        xnf.at(i) = xnf.at(i-1)+dx.at(i-1);

    set_iLower_iUpper();        // recover the bounds

    //----------- 2.5 rule

    impose2point5rule();        // modifies dx, ngrd, ngrdf (but xnf is now inconsisten)

    xnf.resize(ngrdf);
    xnf.at(0) = odtl->posf->d.at(0);               // recover the grid
    for(int i=1; i<ngrdf; i++)
        xnf.at(i) = xnf.at(i-1)+dx.at(i-1);

    //---------- Apply new grid to odt line:

    //---------- split odtl where new grid has faces that don't match odtl.
    // odtl     |          |          |          |          |          |
    // xnf      |     ::   |          |    : :                   :     |
    // odtl new |     ||   |          |    | |   |          |    |     |

    vector<double> whereSplit;
    int iCell;

    double rtol = 1.0E-10*odtl->odtp->domainLength/odtl->ngrd;

    for(int j=1; j<(int)xnf.size()-1; j++) {         // loop over new grid cell faces

        iCell = odtl->linePositionToIndex(xnf.at(j), true, 16);

        if(abs(odtl->posf->d.at(iCell) - xnf.at(j)) > rtol) {
            if(abs(odtl->posf->d.at(iCell+1) - xnf.at(j)) > rtol) { // far away from next face
                // this second if is needed due to the fact that somtimes there is
                // a face in xnf that's exact equal to one face in odtl->posf. In
                // this case a cell with dx=0 is generated which may cause trouble.
                whereSplit.clear();
                whereSplit.push_back( xnf.at(j) );
                j++;
                for( ; j<(int)xnf.size()-1; j++) {

                    if( (abs(odtl->posf->d.at(iCell+1) - xnf.at(j)) > rtol)
                        && (odtl->posf->d.at(iCell+1) > xnf.at(j)) ) {

                        whereSplit.push_back( xnf.at(j) );
                    }
                    else{
                        j--;
                        break;
                    }
                }
                // Append the cell left face and right face to whereSplit
                // We won't split there, but it aids the splitCell functions
                whereSplit.insert(whereSplit.begin(), odtl->posf->d.at(iCell));
                whereSplit.push_back(odtl->posf->d.at(iCell+1));
                splitCell(iCell, whereSplit.size()-2, whereSplit);
            }
        }
    }

    //---------- merge odtl where faces present that are not in new grid
    // odtl new |     ||   |          |    | |   :          :    |     |
    // xnf      |     ||   |          |    | |                   |     |
    // odtl new |     ||   |          |    | |                   |     |
    // now odtl matches xnf

    //---------- note, merging cells may change the cell size for conservation, hence all faces move
    // so mark all, then merge

    mark.clear();

    for(int i=1, j=1; i<odtl->ngrdf-1; i++) {        // loop over interior faces of odtl
        if( abs(odtl->posf->d.at(i)-xnf.at(j)) > rtol)
            mark.push_back(i-1);
        else
            j++;
        }

    for(int i=mark.size()-1; i>=0; i--)         // go backwards so we don't have to update indices in mark
        merge2cells(mark.at(i), false);         //     as we merge cells
                                                // NOTE: flag true here will result in nonconservative mass/mass
                                                //       flow in merging, but elininates cell volumes changing accross
                                                //       the whole domain due to mixing. The meshing should keep the
                                                //       profiles the same and just change grid points, especially
                                                //       if it means the grid moves in unphysical ways (well, its not
                                                //       unphysical since its rigorous mixing when false).
                                                //       We are really trading one numerical error for another:
                                                //       shifting versus mass conservation.
                                                //       Try testing the temporal jet with diffusion bypassed
                                                //       (return at the top of advanceODT).

	enforceDomainSize();
}


/////////////////////////////////////////////////////////////////////////////
/** Finds iLower and iUpper based on the new cell face positions
 *
 */
void meshManager::set_iLower_iUpper() {

    iLower = 0;
    iUpper = ngrd-1;

    if(posLower==xnf.at(ngrd) || posUpper==xnf.at(0)) {
        *odtl->io->ostrm << posLower << " " << posUpper << " " << xnf.at(0) << " " << xnf.at(ngrd) << endl;
        *odtl->io->ostrm << endl << "ERROR in adapMesh posLower/Upper on wrong boundary" << endl;
        exit(0);
    }
    for(int i=0; i<ngrd; i++) {
        if(posLower < xnf.at(i+1)) {
            iLower = i;
            break;
        }
    }
    for(int i=ngrd-1; i>=0; i--) {
        if(posUpper > xnf.at(i)) {
            iUpper = i;
            break;
        }
    }

}

///////////////////////////////////////////////////////////////////////////////
/** Cells that are too small are merged with their neighbors.
 *  Cells here are marked simply as whether they are small or not.
 *  The subsequent merging is with either neighbor.  Elsewhere marked
 *  cells are to be merged with the right neighbor.
 *  Small cells are usually generated by triplet map events that
 *  compress the grid (or by accelerating flows in spatial domains that compress the grid).
 *  Mark the small cells by filling mark (member) with indicies in grid
 *  Merge the small cells in size order.  That is, initial size order, not
 *  accounting for intermediate changes.  This eliminates directional bias too.
 *  So rearrange the small cell indicies into sorted order, scMark (local).
 *  Then loop over each small cell till big enough, merging with its smaller neighbor.
 *  The smaller neighbor may or may not be a small cell.  When merge cells,
 *  need to decrement larger indicies (cells) in the scMark array, and delete cells
 *  from scMark if a small cell is merged with another.
 *  No special treatment for periodic boundaries.  Periodic is the same as nonperiodic.
 *  The periodic boundary position is maintained, and no merging occurs across
 *  a periodic boundary (which would move the boundary from zero).
 */
void meshManager::mergeSmallCells() {

    int i, j;

    //---------- mark the small cells

    mark.clear();
    for(i=0; i<=ngrd-1; i++)
        if(dx.at(i) < odtl->odtp->dxmin)
            mark.push_back(i);

    //---------- sort the small cells to remove any directional bias

    vector<int> ind(mark.size());           // index map for the sort
    for(i=0; i<(int)ind.size(); i++)        // sort "ind" using dx as sort condition
        ind.at(i)=i;
    sort(ind.begin(), ind.end(), *this);    // *this invokes the functor "operator()"
    vector<int> scMark(mark.size());        // reorder the mark array
    for(i=0; i<(int)ind.size(); i++)        // could just use mark.at(ind.at(i)),
        scMark.at(i) = mark.at(ind.at(i));           // but reorder for simplicity


    //---------- merge small cells

    int isc;                                // small cell index
    int iSmallerNB;                         // smaller of 2 neighbors
    int iend = scMark.size();               // may change as merge cells
    bool LcellDone;                         // flag to repeat small cells till big

    //---------------------------------

    for(i=0; i<iend; i++) {                 // loop over small cells

        isc = scMark.at(i);
        isc = scMark.at(i);
        LcellDone = false;

        //---------------------------------

        while(!LcellDone) {                 // keep going till small cell is done

            if(isc == 0)
                iSmallerNB = 1;
            else if(isc == ngrd-1)
                iSmallerNB = ngrd-2;
            else
                iSmallerNB = (dx.at(isc-1) <= dx.at(isc+1)) ? isc-1 : isc+1;

            if(dx.at(iSmallerNB) + dx.at(isc) >= odtl->odtp->dxmin)
                LcellDone = true;           // if new cell big enough then done

            //---------------------------------

            if(iSmallerNB > isc) {

                // merge cells

                ngrdf--;
                ngrd--;
                dx.at(isc) = dx.at(isc) + dx.at(isc+1);
                dx.erase( dx.begin() + isc+1 );

                // delete cell iSmallerNB from scMark array (if present)
                // and decrement all cells > isc in scMark

                for(j=i+1; j<iend; j++)
                    if(scMark.at(j) > isc)
                        scMark.at(j)--;
                for(j=i+1; j<iend; j++)
                    if(scMark.at(j)==isc){
                        scMark.erase(scMark.begin()+j);
                        iend--;
                        break;
                    }

            }         // if(iSmallerNB > isc)

            //---------------------------------

            else {

                // merge cells

                ngrdf--;
                ngrd--;
                dx.at(isc-1) = dx.at(isc-1) + dx.at(isc);
                dx.erase( dx.begin() + isc );

                // decrement isc and scMark.at(i) which is isc
                // delete cell isc-1 in scMark if present
                // and decrement

                scMark.at(i)--;
                isc--;

                for(j=i+1; j<iend; j++)
                    if(scMark.at(j) > isc)
                        scMark.at(j)--;
                for(j=i+1; j<iend; j++)
                    if(scMark.at(j)==isc){
                        scMark.erase(scMark.begin()+j);
                        iend--;
                        break;
                    }

            }     // else (i.e., iSmallerNB < isc)
            //---------------------------------
        }         // end while
        //---------------------------------
    }             // end loop over small cells

    //---------- update the eddy regions

} // end function


//////////////////////////////////////////////////////////////////////////////////
/** This function does the same as mergeSmallCell() excapt that it checks, in
 *  addition, the phase of the merged cell. If there is no opportunity for
 *  merging the cell, it is skipped.
 *
 *  ToDo:
 *  If there is a too small cell and both neigbor cells have different phases,
 *  the cell has to be handeled. For this case, I have no idea how to handle.
 */
void meshManager::mergeSmallCellsMP() {

    double x       = odtl->posf->d.at(0);
    double xp      = 0;                    // midd point of previous cell
    double xc      = 0;                    // midd point of current cell
    double xn      = 0;                    // midd point of next cell
    double phase_p = 0;
    double phase_c = 0;
    double phase_n = 0;


    for (int i=0; i<ngrd; i++){
        if (dx.at(i) >= odtl->odtp->dxmin)
            x += dx.at(i);
        else{
            // the current cell is too small

            // if it is the first cell
            if (i == 0){
                xc = x + dx.at(i)/2;
                xn = x + dx.at(i)+dx.at(i+1)/2;
                phase_c = odtl->phase->d.at( odtl->linePositionToIndex( xc, true, 17 ) );
                phase_n = odtl->phase->d.at( odtl->linePositionToIndex( xn, true, 18 ) );

                if (phase_c == phase_n){
                    ngrdf--;
                    ngrd--;
                    dx.at(i) = dx.at(i) + dx.at(i+1);
                    dx.erase( dx.begin() + i+1 );
                    i--;
                } else{
                    x += dx.at(i);
                }
                continue;
            }

            // if it is the last cell
            if (i == ngrd-1){
                xp = x - dx.at(i-1)/2;
                xc = x + dx.at(i)/2;
                phase_p = odtl->phase->d.at( odtl->linePositionToIndex( xp, true, 19 ) );
                phase_c = odtl->phase->d.at( odtl->linePositionToIndex( xc, true, 20 ) );

                if (phase_p == phase_c){
                    ngrdf--;
                    ngrd--;
                    x -= dx.at(i-1); // subtracting the dx of the previous cell, added in the last loop step
                    dx.at(i-1) = dx.at(i-1) + dx.at(i);
                    dx.erase( dx.begin() + i );
                    i--; i--;
                } else{
                    x += dx.at(i);
                }
                continue;
            }

            // it is a cell in the midd of the domain
            xp = x - dx.at(i-1)/2;         // cell center of previous cell
            xc = x + dx.at(i)/2;           // cell center of current cell
            xn = x + dx.at(i) + dx.at(i+1)/2; // cell center of next cell
            phase_p = odtl->phase->d.at( odtl->linePositionToIndex( xp, true, 21 ) );
            phase_c = odtl->phase->d.at( odtl->linePositionToIndex( xc, true, 22 ) );
            phase_n = odtl->phase->d.at( odtl->linePositionToIndex( xn, true, 23 ) );

            // the case, that the previous cell has another phase
            if (phase_p != phase_c && phase_c == phase_n){
                // the next cell has to be taken for merging cells
                ngrdf--;
                ngrd--;
                dx.at(i) = dx.at(i) + dx.at(i+1);
                dx.erase( dx.begin() + i+1 );
                i--;
                continue;
            }

            // the case, that the next cell has another phase
            if (phase_p == phase_c && phase_c != phase_n){
                // the previous cell has to be taken for merging cells
                ngrdf--;
                ngrd--;
                x -= dx.at(i-1); // subtracting the dx of the previous cell, added in the last loop step
                dx.at(i-1) = dx.at(i-1) + dx.at(i);
                dx.erase( dx.begin() + i );
                i--; i--;
                continue;
            }

            // the case, that both neighbor cells have different phases
            if (phase_p != phase_c && phase_c != phase_n){
                x += dx.at(i);
                continue;
            }

            // the last case, that all three cells have the same phase
            if (phase_p == phase_c && phase_c == phase_n){
                // the smaller one of both neighbor cells has to be taken
                if (dx.at(i-1) < dx.at(i+1)){
                    // the previous neigbor has to be taken
                    ngrdf--;
                    ngrd--;
                    x -= dx.at(i-1); // subtracting the dx of the previous cell, added in the last loop step
                    dx.at(i-1) = dx.at(i-1) + dx.at(i);
                    dx.erase( dx.begin() + i );
                    i--; i--;
                } else{
                    // the next neigbor has to be taken
                    ngrdf--;
                    ngrd--;
                    dx.at(i) = dx.at(i) + dx.at(i+1);
                    dx.erase( dx.begin() + i+1 );
                    i--;
                }
                continue;
            } else{
                *odtl->io->ostrm << "\n\ni = " << i;
                *odtl->io->ostrm <<   "\nphase_p = " << phase_p;
                *odtl->io->ostrm <<   "\nphase_c = " << phase_c;
                *odtl->io->ostrm <<   "\nphase_n = " << phase_n;
                *odtl->io->ostrm <<   "\n(phase_p != phase_c && phase_c == phase_n) = " << (phase_p != phase_c && phase_c == phase_n);
                *odtl->io->ostrm <<   "\n(phase_p != phase_c && phase_c == phase_n) = " << (phase_p == phase_c && phase_c != phase_n);
                *odtl->io->ostrm <<   "\n(phase_p != phase_c && phase_c == phase_n) = " << (phase_p != phase_c && phase_c != phase_n);
                *odtl->io->ostrm <<   "\n(phase_p != phase_c && phase_c == phase_n) = " << (phase_p == phase_c && phase_c == phase_n);
                *odtl->io->ostrm << "\nERROR: In the function meshManager::mergeSmallCells_MP()"
                << endl << "This error might not be there caused by logic one of the"
                << endl << "previous cases have to occur."
                << endl << "There must be a NaN problem.";
                exit(0);
            }
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////
/**
 * Imposes the 2.5 rule for Mesh adaption to the grid.
 */

void meshManager::impose2point5rule(){

    int    i;
    double d1;

    //---------- set the mark array with offending cells

    mark.clear();
    int iHi = (iUpper < ngrd-1) ? iUpper   : iUpper-1;
    int iLo = (iLower > 0)      ? iLower-1 : iLower;
    //int iHi = ngrd-2;
    //int iLo = 0;
    for(i=iLo; i<=iHi; i++) {
        d1 = dx.at(i)/dx.at(i+1);
        if(d1 > 2.5 || d1 < 0.4)
            mark.push_back(i);
    }
    i = ngrd-1;
    if(odtl->odtp->Lperiodic && iUpper==i) {
        d1 = dx.at(i)/dx.at(0);
        if(d1 > 2.5 || d1 < 0.4)
            mark.push_back(i);
    }

    //-------- loop over each marked cell and fix the 2.5 offenders

    for(i=0; i<(int)mark.size(); i++)
        fix2point5offender(mark.at(i), i);

}

///////////////////////////////////////////////////////////////////////////////
/** Marks cells that offend their neighbors with respect to the 2.5 rule
 *
 *           Routine is recursive.  Returns if the test is passed, or if we are at
 *                the ends of the domain.                                                \n
 *           If it starts with no "small cells" routine, it won't create any small cells.\n
 *           The order of operations doesn't matter.                                     \n
 *           Find the larger of mPos and its next neighbor (mPos+1) and split that
 *                cell on the side of the smaller cell.  Split it nsplit times in half.  \n
 *                So, <tt><pre>
 *                      |               ||  ==>  |   *   :    :   || </pre></tt>
 *                when splitting twice.\n
 *           Then go to the opposite side of the split cell (marked with \c * above) and
 *                compare to its neighbor on the other side.  In this example, we would
 *                "mark" the cell before the \c * and call the same routine on it.          \n
 *           Hence, we can traverse through the domain.  Generally you keep traversing
 *                till satisfy the rule or hit the edge.  To keep going in one direction
 *                the cells need to keep increasing (so we won't go too far since motion
 *                implies geometric growth).  You also will stop traversing when you run into
 *                a smaller cell.                                                        \n
 *           Consider:\n <tt><pre>
 *                      ||        |            ||            ||                                  |
 *                       0    1          2     3       4     5                  6
 *           mark offenders: 0 2 3 4 5
 *           Start at 0:
 *           0 vs 1 --> split 1 once
 *                      ||    :   |            ||            ||                                  |
 *                       0  1   2        3     4       5     6                  7
 *           2 vs 3 --> split 3 once
 *                      ||    :   |      :     ||            ||                                  |
 *                       0  1   2     3     4  5       6     7                  8
 *           4 vs 5 --> split 4 once
 *                      ||    :   |      :  :  ||            ||                                  |
 *                       0  1   2     3   4  5 6       7     8                  9
 *           3 vs 4 --> pass --> done
 *           now do next in mark array (which was updated as we went to ( X 3 6 7 8)
 *           Start at 3:
 *           3 vs 4 --> pass (trivially) --> done
 *           Start at 6:
 *           6 vs 7 --> split 7 once
 *                      ||    :   |      :  :  ||      :     ||                                  |
 *                       0  1   2     3   4  5 6   7      8  9                  10
 *           8 vs 9 --> pass (say) --> done:  mark array ( X X X 7 9 )
 *           Start at 7:
 *           7 vs 8 --> pass
 *           Start at 9:
 *           9 vs 10 --> split 4 times
 *                      ||    :   |      :  :  ||      :     || :  :    :       :                |
 *                       0  1   2     3   4  5 6   7      8  9 10 11 12    13           14
 *           14 at the end --> DONE
 *           </pre></tt>
 *
 * @param mPos    \input  marks a cell that offends its next neighbor with respect to the 2.5 rule.
 * @param iglobal \input is the index of which loop it is on in the mark arrary,
 *       for the mark array update. (Could just leave this out and do all of them
 *       but only need to update subsequent ones).
 *
 */
void meshManager::fix2point5offender(const int mPos, const int &iglobal) {


    //--------- Simply return if you pass the test

    int ip = (mPos == ngrd-1) ? 0 : mPos+1;
    double ratio = (dx.at(mPos) < dx.at(ip)) ? dx.at(ip)/dx.at(mPos) : dx.at(mPos)/dx.at(ip);
    if(ratio < 2.5)
        return;

    //--------- Split the larger of the two cells

    int     i;
    int     isplt = (dx.at(mPos) < dx.at(ip)) ? ip : mPos;          // split larger of 2 cells
    int     nsplt = static_cast<int>(log2(ratio / 2.5)) + 1 ; // how many splits to do

    vector<double> icp(nsplt);                 // pos of new intrnl cell fcs (rel to lft edge)

    bool splitCellOnRight = (isplt > mPos);
    if(mPos == ngrd-1 && isplt == 0)           // for periodic
        splitCellOnRight = true;

    if(splitCellOnRight) {                     // | |  :  :    :        |  --> nsplt=3
        icp.at(nsplt-1) = dx.at(isplt)*0.5;
        for(i=nsplt-2; i>=0; i--)
            icp.at(i) = icp.at(i+1)*0.5;

        for(i=nsplt-1; i>=1; i--)              // update dx arr, reuse icp arr
            icp.at(i) -= icp.at(i-1);                // i.e. we inserted cells, so insert dx
        dx.at(isplt) *= 0.5;                      // icp was positions, now are dx's
        dx.insert(dx.begin()+isplt, icp.begin(), icp.end());
    }
    else {                                     // |        :    :  :  | |   --> nsplt=3
        icp.at(0) = dx.at(isplt)*0.5;
        for(i=1; i<nsplt; i++)                 // set icp as dx's : 1/2, 1/4, 1/8 ...
            icp.at(i) = icp.at(i-1)*0.5;
        for(i=1; i<nsplt; i++)                 // now offset to 0.5, 0.75, 0.876 ...
            icp.at(i) += icp.at(i-1);

        for(i=0; i<nsplt-1; i++)               // update dx arr, reuse icp arr
            icp.at(i) = icp.at(i+1) - icp.at(i);
        if(nsplt >=2)
            icp.at(nsplt-1) = icp.at(nsplt-2);
        dx.at(isplt) *= 0.5;
        dx.insert(dx.begin()+isplt+1, icp.begin(), icp.end());
    }

    ngrd  += nsplt;                            // update meshManager's ngrd, ngrdf (odtl's are done)
    ngrdf += nsplt;

    //--------- Update the mark array

    for(i=iglobal+1; i<(int)mark.size(); i++)
        if(mark.at(i) > isplt)
            mark.at(i) += nsplt;

    //--------- Now compare the other half of the split cell with neighbor

    int inext;
    if(splitCellOnRight) {
        inext = isplt + nsplt;
        if(inext == ngrd-1 && !odtl->odtp->Lperiodic)
            return;
    }
    else {
        inext = mPos-1;
        if(inext == -1) {
            if(!odtl->odtp->Lperiodic)
                return;
            else
                inext = ngrd-1;
        }
    }

    fix2point5offender(inext, iglobal);                 // recursive call

}

///////////////////////////////////////////////////////////////////////////////
/** Resizes and sets the dx array to reflect the correct distances.
 *
 */
void meshManager::setDxArray() {

    dx.resize(ngrd);
    for(int i=0, ip=1; i<ngrd; i++, ip++)
        dx.at(i) = xnf.at(ip) - xnf.at(i);

}


////////////////////////////////////////////////////////////////////////////////
/** Given a position, find an index (use with pos)
 *
 * @param x      \input vector of cell positions
 * @param val    \input location which is converted to index
 * @param istart \input start looking at this index
 * @return index of position
 */
int meshManager::findPos(const vector<double> &x, const double val, const int &istart) {

    if(val <= x.at(0))
        return 0;
    if(val >= x.at(x.size()-1))
        return x.size()-2;
    for(int i=istart; i<(int)x.size(); i++)
        if(x.at(i) >= val) {
            return i-1;
        }
    return 0;
}


////////////////////////////////////////////////////////////////////////////////
/** Linear interpolation of a single point, given two vectors \c x and <code>y</code>.
 *
 * @param x      \input  vector of cell positions
 * @param y      \input  vector of cell values
 * @param xval   \input  position to interpolate
 * @param yval   \output value at interpolated position
 * @param istart \output start looking at this index (for findPos)
 */
void meshManager::interp1pt(const vector<double> &x, const vector<double> &y,
                            const double &xval, double &yval, int &istart) {

    int i, ip;

    i = findPos(x, xval, istart);
    ip = i+1;
    if(i>0)
        istart = i-1;

    yval = y.at(i) + (xval-x.at(i))*(y.at(ip)-y.at(i))/(x.at(ip)-x.at(i));

}


////////////////////////////////////////////////////////////////////////////////
/** Linear interpolation of a vector of points (call interp1pt each time)
 *
 * @param x      \input  vector of cell positions
 * @param y      \input  vector of cell values
 * @param xn     \input  vector of positions to interp to
 * @param yn     \output vector of interpolated values
 */
void meshManager::interpVec(const vector<double> &x, const vector<double> &y,
                            const vector<double> &xn, vector<double> &yn) {

    if(x.size() != y.size() || xn.size() != yn.size()) {
        cerr << "\nERROR IN INTERPVEC" << endl;
        exit(0);
    }
    int istart = 0;
    for(int i=0; i<(int)xn.size(); i++)
        interp1pt(x,y,xn.at(i),yn.at(i), istart);
}


////////////////////////////////////////////////////////////////////////////////
/** Computes and returns total distance (length) along a curve specified by vectors x and y
 *
 * @param x      \input  vector of cell positions
 * @param y      \input  vector of vector cell values for profiles (phi's) to compare
 * @param sDist  \input  vector that stores the running distance along the curve
 *
 * @return the total distance(length) along the specified curve.
 */
double meshManager::calcDistance(const vector<double> &x, const vector<vector<double> > &y,
                                 vector<double> &sDist) {

    double dx2, dy2;    // dx^2 and dy^2
    double dmb;

    sDist.at(0) = 0.0;
    for (int i=1; i<(int)x.size(); i++) {

        dx2 = x.at(i)-x.at(i-1);
        dx2 *= dx2;

        dy2 = 0.0;
        for(int iProf=0; iProf<(int)y.size(); iProf++) {
            dmb = y.at(iProf).at(i) - y.at(iProf).at(i-1);
            dmb *= dmb;
            if(dmb > dy2) dy2 = dmb;
        }

        sDist.at(i)= sDist.at(i-1)+sqrt(dx2 +dy2);
    }

    return sDist.at(sDist.size()-1);
}


////////////////////////////////////////////////////////////////////////////////
/** Operates on the odtline.
 *  Splits cells at interior elements of cellFaces vector (conservative).
 *  @param isplt  \input cell to split
 *  @param nsplt  \input number of times to split the cell
 *  @param cellFaces \input locations of original left edge, new interior faces, and right edge
 */

void meshManager::splitCell(const int isplt, const int nsplt, const vector<double> &cellFaces) {

    for(int i=0; i< odtl->v.size(); i++)
        odtl->v.at(i)->splitCell(isplt, nsplt, cellFaces);

    odtl->ngrd  += nsplt;
    odtl->ngrdf += nsplt;

}

////////////////////////////////////////////////////////////////////////////////
/** Operates on the odtline
 *  Merges conservatively.
 *  @param imrg \input merge imrg with imrg+1
 *  @param constVolume \input if true, do the merge without allowing a change in cell volume.
           (default is false. True when we merge small cells at edgees after enforcing domain length
           so that the volume doesn't change, requiring another enforce domain length, etc.)

    Below, the notation is |---1---|----2----|  --> |---------3-------|

    For temporal flows:
        Merge with constant volume before and after.
        * V3 = V1 + V2 (V = volume = dxc = delta(x^c)).
        * Let m1 = rho1*V1 and m2 = rho2*V2
        * (Mass conservation implies: m3 = m1 + m2)
        * For all scalars except rho: s3 = (s1*m1 + s2*m2) / (m1+m2)
        * Then set rho3 based on these scalars (like rho=MP/RT, say)
        * Set posf and pos for V3 (for posf, just delete the middle face).

        This will not preserve mass conservation: V3 is set, and rho3 is computed from
            rho=MP/RT instead of from rho3 = (m1+m2)/V3.
            (Or, if mass is conserved, it implies a pressure/temperature change).

        Version without conserving volume:
        * Let m1 = rho1*V1 and m2 = rho2*V2
        * For all scalars except rho: s3 = (s1*m1 + s2*m2) / (m1+m2)
        * Then set rho3 based on these scalars (like rho=MP/RT, say)
        * Set V3 = (m1+m2)/rho3.
        * Set posf and pos for V3 (domain expands or contracts and the whole profiles shift).


    For spatial flows:
        Merge with constant volume before and after.
        * V3 = V1 + V2
        * Let m1 = rho1*u1*V1 and m2 = rho2*u2*V2
        * (Mass flux conservation implies: m3 = m1 + m2)
        * Set all scalars except rho and u: s3 = (s1*m1 + s2*m2)/(m1+m2).
                (Well, do u, but it will be redundant as overwritten below).
        * Set rho3 based on scalars (like rho=MP/RT, say).
        * Set u3 = (m1+m2)/(rho3*V3).      (This step conserves mass flux).
        * Set posf and pos for V3 (for posf, just delete the middle face).

        This will conserve mass flux and all scalar fluxes, but not momentum flux,
        because we required volume conservation.

        Version without conserving volume:
        * Let m1 = rho1*u1*V1 and m2 = rho2*u2*V2
        * Set all scalars except rho: s3 = (s1*m1 + s2*m2)/(m1+m2), INCLUDING u3.
              (Now momentum is conserved).
        * Set rho3 based on scalars (like rho=MP/RT, say).
        * Set V3 = (m1+m2)/(rho3 * u3)
        * Set posf and pos for V3. Here, expand/contract about the expansion center.


        This conserves everything, but does "bad" things:
            * If you change a cell volume, then you shift the profile throughout the domain.
            * That is especially bad if all you are trying to do is Merge two cells together.
            * Example: velocity profile.
                  (This was done with a toy problem in the code: set prof, merge cells).
                           4
                          / \
                         /   \
                        /     \
                       3       5
                      /         \
                     /           \
              1-----2             6----7----8
                             4
                            / \
                           /   \
                          /     \
                         3       5
                        /         \
                       /           \
                1-----2             6----7

                  - Suppose you merge cell 5 with cell 6.
                  - Let the u=0 in cells 1, 2, 6-8.
                  - Assume constant density.
                  - Then m2=0, so u3 = u1, That is the merged cell has u = u5.
                  - Also, V3 = V1 and we delete cell V2 (that is, cell 6),
                         so there is a net loss of volume on the domain.
                  - We remove half the volume from each side of the domain,
                         which is equivalent to contracting about the contraction center.
                  - In a planar flow, this means we simply shift the whole profile to the right,
              * We would then extend the end cells to fill the domain.
              * But this is clearly weird. By merging cells, we have shifted the whole domain and
                   not changed the profile shape.
              * If you merge cell 1 with two, it moves the other direction.
              * Now, normally you can't have a zero velocity in spatial flow, but with "small" velocity,
              *   similar behavior happens.
              * This shifting of the domain caused profiles to wander unphysically during simulation.
                   A trip map on one side will cause more merging/splitting on that side...

 NOTE: not conserving moment (constant volume formulation) causded major momentum conservation issues in jet flows.
 */

void meshManager::merge2cells(const int imrg, const bool LconstVolume) {

    double dxc1, dxc2;               // cell volumes (imrg, imrg+1, final merged cell) = delta(x^c)
    double pm1;                      // +1 or -1 to change sign

    pm1  = odtl->posf->d.at(imrg+1)*odtl->posf->d.at(imrg) < 0 ? -1.0 : 1.0;              // handles cells that split x=0
    dxc1 = abs( pow(abs(odtl->posf->d.at(imrg+1)),     odtl->odtp->cCoord) -
                pm1*pow(abs(odtl->posf->d.at(imrg) ),  odtl->odtp->cCoord) );
    pm1  = odtl->posf->d.at(imrg+2)*odtl->posf->d.at(imrg+1) < 0 ? -1.0 : 1.0;
    dxc2 = abs( pow(abs(odtl->posf->d.at(imrg+2)),      odtl->odtp->cCoord) -
                pm1*pow(abs(odtl->posf->d.at(imrg+1) ), odtl->odtp->cCoord) );

    double m1 = dxc1 * odtl->rho->d.at(imrg);
    double m2 = dxc2 * odtl->rho->d.at(imrg+1);

    if(odtl->odtp->Lspatial) {
        m1 *= odtl->uvel->d.at(imrg);             // m is now mdot
        m2 *= odtl->uvel->d.at(imrg+1);
    }

    //------- merge all but rho, pos, posf (if spatial and LconstVolume this step is redundant for uvel)

    for(int k=0; k<odtl->v.size(); k++) {
        if(odtl->v.at(k)->var_name=="rho" || odtl->v.at(k)->var_name=="pos" || odtl->v.at(k)->var_name=="posf")
            continue;
        odtl->v.at(k)->merge2cells(imrg, m1, m2, LconstVolume);
    }
    //------- now finish up

    odtl->rho->merge2cells( imrg, m1, m2, LconstVolume);
    if(odtl->odtp->Lspatial && LconstVolume) {
        odtl->uvel->d.at(imrg) = (m1 + m2) / (odtl->rho->d.at(imrg) * (dxc1+dxc2));
    }
    odtl->ngrd--;
    odtl->ngrdf--;
    odtl->posf->merge2cells(imrg, m1, m2, LconstVolume);
    // odtl->pos is done in posf on previous line

}

////////////////////////////////////////////////////////////////////////////////
/**Adapt regions of the mesh depending on age since last adaption.
 * Start at cell index cLastDA in the "lastDA" domain (grid) and search right and left
 * for old cells.  Update these cells (tag them with current time), then adapt them.
 * Then cLastDA (from Update) is the next oldest cell not just adapted.  Search left
 * and right for a new region, and so on until nothing on the grid is too old.
 *
 * @param time    \input current time.
 * @param tLastDA \inout time of last diffusion advancement.
 * @param cLastDA \inout cell index.
 * @param dtCUmax \inout the time increment of eddy trial time advancement before we diffuse
 */

void meshManager::adaptAfterSufficientDiffTime (const double &time,
                                                double       &tLastDA,
                                                int          &cLastDA,
                                                double       &dtCUmax) {


    if(time - tLastDA < odtl->odtp->DAtimeFac * dtCUmax)
        return;

    int iStart, iEnd;

    double youngsters = odtl->odtp->DAtimeFac*dtCUmax;

    while(time - tLastDA >= youngsters) {

        //---------- initialize searches

        int leftSide  = cLastDA;
        int rightSide = (cLastDA == odtl->odtp->sLastDA-1) ? cLastDA : cLastDA+1;

        //---------- search for leftSide

        while(time - lastDA.at(leftSide) >= 0.5 * youngsters){
            if(leftSide == 0)
                break;
            leftSide--;
        }

        //---------- correct for possible overshoot

        if(time - lastDA.at(leftSide) < 0.5 * youngsters)
            leftSide++;

        //---------- search for rightSide

        while(time - lastDA.at(rightSide) >= 0.5 * youngsters){
            if(rightSide == odtl->odtp->sLastDA-1)
                break;
            rightSide++;
        }

        //---------- correct for possible overshoot

        if(time - lastDA.at(rightSide) < 0.5 * youngsters)
            rightSide--;

        //---------- find the corresponding range of indices of the posf array

        iStart = odtl->linePositionToIndex(odtl->posf->d[0]+leftSide *odtl->Ldomain()/odtl->odtp->sLastDA, true, 24);
        iEnd   = odtl->linePositionToIndex(odtl->posf->d[0]+rightSide*odtl->Ldomain()/odtl->odtp->sLastDA, false, 25);

        updateDA(time, tLastDA, cLastDA, iStart, iEnd);

        //---------- adapt the line

        adaptGrid(iStart, iEnd);
    }
}

///////////////////////////////////////////////////////////////////////////////
/**Used for selective mesh adaptation after a sufficient time.
 * Write the current value of time into the cells of lastDA that overlap or are
 * contained in the mesh adaption interval.  If cLastDA is one of those cells, scan the lastDA
 * array to find the smallest time value in the array.  Set tLastDA equal to that
 * smallest time value and set cLastDA equal to the index of the cell containing
 * that smallest value.
 *
 * @param time    \input current time.
 * @param tLastDA \inout time of last diffusion advancement.
 * @param cLastDA \inout cell index.
 * @param iStart  \input starting cell index to check.
 * @param iEnd    \input ending cell index to check.
 */

void meshManager::updateDA(const double &time, double &tLastDA, int &cLastDA,
                           int iStart, int iEnd) {

    //----------- find the range of affected lastDA cells

    int leftCell  = static_cast<int>(odtl->odtp->sLastDA * (odtl->posf->d.at(iStart)-odtl->posf->d.at(0)) / odtl->Ldomain());
    if(leftCell==odtl->odtp->sLastDA) leftCell--;
    int rightCell = static_cast<int>(odtl->odtp->sLastDA * (odtl->posf->d.at(iEnd+1)-odtl->posf->d.at(0)) / odtl->Ldomain());
    if(rightCell == odtl->odtp->sLastDA) rightCell--;

    if(leftCell > rightCell && !odtl->odtp->Lperiodic)
        return;            // no affected cells of LastDA array

    //---------- set the lastDA values in this range equal to the current time

    if(leftCell > rightCell) {                    // region like | * | * |   |   | * | * |
        for(int i=leftCell; i<odtl->odtp->sLastDA; i++)
            lastDA.at(i) = time;
        for(int i=0; i<=rightCell; i++)
            lastDA.at(i) = time;
        if( !(cLastDA >= leftCell) && !(cLastDA <= rightCell) )
            return;                               // oldest cell is unaffected
    }
    else {                                        // region like |   | * | * | * | * |   |
        for(int i=leftCell; i<=rightCell; i++)
            lastDA.at(i)=time;
        if(!( (cLastDA >= leftCell) && (cLastDA <= rightCell) ))
            return;                               // oldest cell is unaffected
    }

    //---------- find oldest cell: cLastDA, and earliest last adapt time

    cLastDA = 0;                                  // initialize cLastDA
    for(int i=1; i<odtl->odtp->sLastDA; i++)
        if(lastDA.at(i) < lastDA.at(cLastDA))
            cLastDA = i;                          // find oldest cell
    tLastDA = lastDA.at(cLastDA);                    // earliest last-adapt time

}

///////////////////////////////////////////////////////////////////////////////
/** Adapt region of mesh and keep track of start/stop indices
 *  @param iStart  \inout starting index of region to adapt.
 *  @param iEnd    \inout ending index of region to adapt.
 *  @param dtCUmax \inout the time increment of eddy trial time advancement before we diffuse
 *  @param time    \input current time.
 *  @param tLastDA \inout time of last diffusion advancement.
 *  @param cLastDA \inout cell index.
 */

void meshManager::adaptEddyRegionOfMesh(const double &time, double &tLastDA, int &cLastDA) {

    int iStart = odtl->linePositionToIndex(odtl->ed->leftEdge,  true, 4);
    int iEnd   = odtl->linePositionToIndex(odtl->ed->rightEdge, false, 5);

    if(iStart > 0)          iStart--;
    if(iEnd < odtl->ngrd-1) iEnd++;

    updateDA(time, tLastDA, cLastDA, iStart, iEnd);

	adaptGrid(iStart, iEnd);
}

///////////////////////////////////////////////////////////////////////////////
/** Enforce the domain size by chopping a domain that is too big,
 *  or expanding a domain that is too small.
 *  This is called in the microMixer where the grid cells expand or contract (updateGrid),
 *  and also after the mesher when we have merged all cells.
 *
 *  Note, this only works for open domains (one or two sides).
 */

void meshManager::enforceDomainSize() {

    if(odtl->odtp->bcType == "WALL" || odtl->odtp->LisFlmlt || odtl->odtp->LisHips)
        return;

    double Ld = odtl->Ldomain();

    //----------------------  OUTFLOW ON BOTH SIDES

    if(odtl->odtp->bcType == "OUTFLOW" || odtl->odtp->bcType == "PERIODIC") {

        double delta_hi;        // right side: positive if too big --> chop; negative if too small --> expand
        double delta_lo;        // left side:  positive if too big --> chop; negative if too small --> expand

        if(odtl->odtp->cCoord==1) {       // planar
            delta_hi = 0.5*(Ld - odtl->odtp->domainLength);
            delta_lo = delta_hi;
        }
        else {                            // cylindrical or spherical
            delta_hi =  odtl->posf->d.at(odtl->ngrd) - 0.5*odtl->odtp->domainLength;
            delta_lo = -odtl->posf->d.at(0)          - 0.5*odtl->odtp->domainLength;
        }

        //----------- Do expansions

        if(delta_hi < 0) {       // too small --> expand last cell
            odtl->posf->d.at(odtl->ngrd) -= delta_hi;
            odtl->pos->d.at(odtl->ngrd-1) = 0.5*(odtl->posf->d.at(odtl->ngrd-1) + odtl->posf->d.at(odtl->ngrd));
        }

        if(delta_lo < 0) {       // too small --> expand first cell
            odtl->posf->d.at(0) += delta_lo;
            odtl->pos->d.at(0) = 0.5*(odtl->posf->d.at(0) + odtl->posf->d.at(1));
        }

        //----------- Do contractions

        vector<double> cellFaces;

        if(delta_hi > 0) {       // too big --> chop domain on right

            //before |----0---|----1---|-----2---------|------3-$----|--4---|---5----|
            //after  |----0---|----1---|-----2---------|---3----$=4==|==5===|===6====|
            // or
            //before |----0---|----1---|-----2---------|------3------$--4---|---5----|
            //after  |----0---|----1---|-----2---------|---3---------$==4===|===5====|
            // works if xSplitHi is at the far right face.

            double xSplitHi = odtl->posf->d.at(odtl->ngrd) - delta_hi;
            double iSplitHi = odtl->linePositionToIndex(xSplitHi, false, 28);
            if(odtl->posf->d.at(iSplitHi+1) != xSplitHi) {
                cellFaces.push_back(odtl->posf->d.at(iSplitHi));
                cellFaces.push_back(xSplitHi);
                cellFaces.push_back(odtl->posf->d.at(iSplitHi+1));
                splitCell(iSplitHi, 1, cellFaces);
            }
            for(int k=0; k<odtl->v.size(); k++)
                odtl->v.at(k)->d.erase(odtl->v.at(k)->d.begin()+iSplitHi+1, odtl->v.at(k)->d.end());
            odtl->posf->d.push_back(xSplitHi);
            // pos is fine
            odtl->ngrd -= (odtl->ngrd - iSplitHi - 1);
            odtl->ngrdf = odtl->ngrd+1;
        }

        if(delta_lo > 0) {       // too big --> chop domain on left

            //before |----0---|----1---|-----$--2------|---3--|---4--|---5--|
            //after  |====0===|====1===|==2==$----3----|---4--|---5--|---6--|
            // or
            //before |----0---|----1---$--------2------|---3--|---4--|---5--|
            //after  |====0===|====1===$--------2------|---3--|---4--|---5--|
            // works if xSplitLo is at the far left face.

            double xSplitLo = odtl->posf->d.at(0) + delta_lo;
            double iSplitLo = odtl->linePositionToIndex(xSplitLo, true, 29);
            if(odtl->posf->d.at(iSplitLo) != xSplitLo) {
                cellFaces.resize(0);
                cellFaces.push_back(odtl->posf->d.at(iSplitLo));
                cellFaces.push_back(xSplitLo);
                cellFaces.push_back(odtl->posf->d.at(iSplitLo+1));
                splitCell(iSplitLo, 1, cellFaces);
            }
            else
                iSplitLo--;
            for(int k=0; k<odtl->v.size(); k++)
                odtl->v.at(k)->d.erase(odtl->v.at(k)->d.begin(), odtl->v.at(k)->d.begin()+iSplitLo+1);
            // pos and posf are fine
            odtl->ngrd -= iSplitLo+1;
            odtl->ngrdf = odtl->ngrd+1;
        }

        //------------- merge first cell and last cell (with neighbors) if small

        if(odtl->posf->d.at(odtl->ngrd) - odtl->posf->d.at(odtl->ngrd-1) < odtl->odtp->dxmin)
            merge2cells(odtl->ngrd-2, false);
        if(odtl->posf->d.at(1) - odtl->posf->d.at(0) < odtl->odtp->dxmin)
            merge2cells(0, false);

    }   // endif OUTFLOW

    //----------------------  Wall on left, outflow on the right

    else if(odtl->odtp->bcType == "WALL_OUT") {

        if(Ld < odtl->odtp->domainLength) {
            odtl->posf->d.at(odtl->ngrd) += odtl->odtp->domainLength - Ld;
            odtl->pos->d.at(odtl->ngrd-1) = 0.5*(odtl->posf->d.at(odtl->ngrd-1) + odtl->posf->d.at(odtl->ngrd));
        }

        else if(Ld > odtl->odtp->domainLength) {

            vector<double> cellFaces;

            //before |----0---|----1---|-----2---------|------3-$----|--4---|---5----|
            //after  |----0---|----1---|-----2---------|---3----$=4==|==5===|===6====|
            // or
            //before |----0---|----1---|-----2---------|------3------$--4---|---5----|
            //after  |----0---|----1---|-----2---------|---3---------$==4===|===5====|
            // works if xSplitHi is at the far right face.
            double xSplitHi = odtl->posf->d.at(0) + odtl->odtp->domainLength;
            double iSplitHi = odtl->linePositionToIndex(xSplitHi, false, 30);
            if(odtl->posf->d.at(iSplitHi+1) != xSplitHi) {
                cellFaces.push_back(odtl->posf->d.at(iSplitHi));
                cellFaces.push_back(xSplitHi);
                cellFaces.push_back(odtl->posf->d.at(iSplitHi+1));
                splitCell(iSplitHi, 1, cellFaces);
            }
            for(int k=0; k<odtl->v.size(); k++)
                odtl->v.at(k)->d.erase(odtl->v.at(k)->d.begin()+iSplitHi+1, odtl->v.at(k)->d.end());
            odtl->posf->d.push_back(xSplitHi);
            // pos is fine
            odtl->ngrd -= (odtl->ngrd - iSplitHi - 1);
            odtl->ngrdf = odtl->ngrd+1;

            //------------- merge last cell if its small

            if(odtl->posf->d.at(odtl->ngrd) - odtl->posf->d.at(odtl->ngrd-1) < odtl->odtp->dxmin) {
                merge2cells(odtl->ngrd-2,false);
            }
        }
    }   // endif WALL_OUT

    else {
        cout << endl << "ERROR: meshManager::enforceDomainSize: not setup for given bcType" << endl;
        exit(0);
    }

}

///////////////////////////////////////////////////////////////////////////////
/**Set the cell volume vectors: dxc = | |x_e|^c - |x_w|^c |.
 * @param line \input pointer to odtline to use (eddyline or odtline)
 * @param dxc \inout vector of cell volumes.
 */

void meshManager::setGridDxc(const odtline *line, vector<double> &dxc, double C) {

    dxc.resize(line->posf->d.size()-1);        // note: posf->d.size()-1 not ngrd due to merge2cells for posf kills face, then calls this func, ngrd is fixed after all vars are merged.
    double pm1;
    for(int i=0, ip=1; i<line->posf->d.size()-1; i++, ip++) {
        pm1 = line->posf->d.at(ip)*line->posf->d.at(i) < 0 ? -1.0 : 1.0;              // handles cells that split x=0
        dxc.at(i) = abs(pow(abs(line->posf->d.at(ip)), C) -
                     pm1*pow(abs(line->posf->d.at(i) ), C));
    }
}

///////////////////////////////////////////////////////////////////////////////
/**Set the cell volume vectors: dx = x_e - x_w
 * @param line \input pointer to odtline to use (eddyline or odtline)
 * @param dx \inout vector of cell sizes.
 */

void meshManager::setGridDx(const odtline *line, vector<double> &dx) {

    dx.resize(odtl->posf->d.size()-1);
    for(int i=0, ip=1; i<dx.size(); i++, ip++)
        dx.at(i)  = abs(odtl->posf->d.at(ip) - odtl->posf->d.at(i));
}

/////////////////////////////////////////////////////////////////////////////////
///** If a face is too close to r=0, remove it by merging.
// *  Only for cylindrical or spherical
// */
//
//void meshManager::removeFaceNearZero() {
//
//    int imid = odtl->linePositionToIndex(0.0, true);
//
//    double dxmid = odtl->posf->d.at(imid+1) - odtl->posf->d.at(imid);
//    if(dxmid <= 3.0*odtl->odtp->dxmin) {
//        if(imid < odtl->ngrd-1)
//            merge2cells(imid,true);
//        if(imid > 0)
//            merge2cells(imid-1,true);
//    }
//    return;
//}

///////////////////////////////////////////////////////////////////////////////
/** If a face is too close to r=0, remove it by merging.
 *  Only for cylindrical or spherical
 */

void meshManager::setGridFromDxc(const vector<double> &dxc2) {

    double C    = odtl->odtp->cCoord;
    double invC = 1.0/odtl->odtp->cCoord;

    //------------- Outflow on both sides

    if(odtl->odtp->bcType=="OUTFLOW" || odtl->odtp->bcType=="PERIODIC" ) {

        //-----------  for planar cases, expand/contract about the expansion/contraction "center"
        if(C==1) {
            double V2tot = accumulate(dxc2.begin(), dxc2.end(), 0.0);
            double dmb;
            odtl->posf->d[0] = -pow(0.5*V2tot, invC);
            for(int ie=1, iw=0, i=0; ie<odtl->ngrdf; ie++, iw++, i++) {
                if(odtl->posf->d[iw] <= 0.0) {
                    dmb = pow(abs(odtl->posf->d[iw]),C) - dxc2[i];
                    if(dmb >=0) odtl->posf->d[ie] = -pow( dmb, invC);
                    else        odtl->posf->d[ie] =  pow(-dmb, invC);
                }
                else {
                    odtl->posf->d[ie] = pow(pow(odtl->posf->d[iw],C) + dxc2[i], invC);
                }
            }
        }

        //-----------  for round cases, expand/contract about the r=0.0 point: FIX zero at zero
        else {
            int i0 = odtl->linePositionToIndex(0.0, true, 31);
            double pm1 = odtl->posf->d.at(i0+1)*odtl->posf->d.at(i0) < 0 ? -1.0 : 1.0;              // handles cells that split x=0
            double dxc0 = abs( pow(abs(odtl->posf->d.at(i0+1)), C) - pm1*pow(abs(odtl->posf->d.at(i0) ), C) );
            double fracVolR =  pow(odtl->posf->d[i0+1],C) / dxc0;
            double fracVolL =  1.0 - fracVolR;

            //-------- Right side of domain

            odtl->posf->d[i0+1] = pow(fracVolR*dxc2[i0], invC);
            for(int i=i0+2, im=i0+1; i<odtl->ngrdf; i++, im++)
                odtl->posf->d[i] = pow(pow(odtl->posf->d[im],C) + dxc2[im], invC);

            //-------- Left side of domain

            odtl->posf->d[i0] = -pow(fracVolL*dxc2[i0], invC);
            for(int i=i0-1, im=i0; i>=0; i--, im--)
                odtl->posf->d[i] = -pow(pow(abs(odtl->posf->d[im]),C) + dxc2[i], invC);
        }
    }

    //------------------ Wall on left, outlet on right. Assume all posf values are >= 0.0, expand to the right.

    else if(odtl->odtp->bcType == "WALL_OUT") {

        odtl->posf->d[0] = 0.0;
        for(int ie=1, iw=0, i=0; ie<odtl->ngrdf; ie++, iw++, i++)
            odtl->posf->d[ie] = pow(pow(odtl->posf->d[iw],C) + dxc2[i], invC);
    }

    //------------------

    else {
		cout << endl << "ERROR: meshManager::setGridFromDxc: not setup for given bcType" << endl;
       	exit(0);
    }

    //------------------

    odtl->pos->setVar();

    //------------------
}

/////////////////////////////////////////////////////////////////////////////////
///** Create cell with r=0 at its center
// *  Only for cylindrical or spherical
// */
//
/*void meshManager::splitCellWithZero() {

    if(odtl->odtp->xDomainCenter == 0.0) {
    int iSplitLo = odtl->linePositionToIndex((-1.0)*odtl->odtp->dxCenter/2.0, true, 65); // in case there is a face position with -dxmax/2, it will be in odtl->posf->d[iSplitLo]
    int iSplitHi = odtl->linePositionToIndex(odtl->odtp->dxCenter/2.0, false, 66);  // in case there is a face position with dxmax/2, it will be in odtl->posf->d[iSplitHi+1]
    int iSplitZero = odtl->linePositionToIndex(0.0, true, 67);  // in case there is a face position with 0, it will be in odtl->posf->d[iSplitZero]

    if(abs(odtl->pos->d.at(iSplitZero)) > 1.0E-10) {

        vector<double> cellFacesLo;
        vector<double> cellFacesHi;

        if( odtl->posf->d.at(iSplitHi+1) != odtl->odtp->dxCenter/2.0 ) {    // in case odtl->posf->d[iSplitHi+1] is not dxmax/2, then split and then odtl->posf->d[iSplitHi+1] will be dxmax/2
        cellFacesHi.push_back(odtl->posf->d.at(iSplitHi));
        cellFacesHi.push_back(odtl->odtp->dxCenter/2.0);
        cellFacesHi.push_back(odtl->posf->d.at(iSplitHi+1));
        splitCell(iSplitHi, 1, cellFacesHi);
        }

        if( odtl->posf->d.at(iSplitLo) != (-1.0)*odtl->odtp->dxCenter/2.0 ) { // in case odtl->posf->d[iSplitLo] is not -dxmax/2, then split and then odtl->posf->d[iSplitLo+1] will be -dxmax/2
        cellFacesLo.push_back(odtl->posf->d.at(iSplitLo));
        cellFacesLo.push_back((-1.0)*odtl->odtp->dxCenter/2.0);
        cellFacesLo.push_back(odtl->posf->d.at(iSplitLo+1));
        splitCell(iSplitLo, 1, cellFacesLo);
        iSplitHi++; // correspondingly, due to the cell added, increase iSplitHi
        iSplitLo++; // due to the cell added, increase iSplitLo, such that odtl->posf->d[iSplitLo] is -dxmax/2
        }

        double endMerge = odtl->posf->d.at(iSplitHi+1);
        while( odtl->posf->d.at(iSplitLo+1) != endMerge ) {     // the merge changes odtl->posf->d[iSplitLo+1], since it involves iSplitLo and iSplitLo+1. odtl->posf->d[iSplitLo] remains unchanged at -dxmax/2
        merge2cells(iSplitLo,true);
        }

        merge2cells(iSplitLo+1, true); // merging iSplitLo+1 with iSplitLo+2 changes odtl->posf->d[iSplitLo+2] and dxmax/2 remains at odtl->posf->d[iSplitLo+1]
        merge2cells(iSplitLo-2, true); // merging iSplitLo-2 with iSplitLo-1 changes odtl->posf->d[iSplitLo-1] and -dxmax/2 changes from odtl->posf->d[iSplitLo] to odtl->posf->d[iSplitLo-1] due to the reduction in the number of cells. Also, dxmax/2 changes from odtl->posf->d[iSplitLo+1] to odtl->posf->d[iSplitLo], thus 0 is at odtl->pos->d.at(iSplitLo-1)

        iNextZero = iSplitLo;
    }
    else
        iNextZero = iSplitZero+1;
    }
    else {
    iNextZero = 0;
    }

}*/
