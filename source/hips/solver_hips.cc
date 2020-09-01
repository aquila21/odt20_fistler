
#include "solver_hips.h"
#include "odtline.h"
#include <cmath>     // pow
#include <algorithm> // copy, sort
#include <sstream>
#include <string>
#include <iomanip>   // precision

///////////////////////////////////////////////////////////////////////////////

bool sortFunc(pair<double,int> &a, pair<double,int> &b){
    return a.first < b.first;
}

///////////////////////////////////////////////////////////////////////////////
/** Initializer
    @param p_odtl \input pointer to odtline object
 */

void solver_hips::init(odtline *p_odtl) {

    odtl = p_odtl;

    //------------------- Set number of parcels, and level lengthscales, timescales, and rates

    odtl->ngrd = static_cast<int>(pow(2.0, odtl->odtp->nLevels-1));

    vector<double> levelLengths(odtl->odtp->nLevels);
    vector<double> levelTaus(odtl->odtp->nLevels);
    levelRates   = vector<double>(odtl->odtp->nLevels);

    for(int i=0; i<odtl->odtp->nLevels; i++) {
        levelLengths[i] = odtl->odtp->domainLength * pow(odtl->odtp->Afac,i);
        levelTaus[i]    = odtl->odtp->tau0 *
                          pow(levelLengths[i]/odtl->odtp->domainLength, 2.0/3.0) /
                          odtl->odtp->C_param;
        levelRates[i]   = 1.0/levelTaus[i] * pow(2.0,i);   // lambda rate at each level (for all nodes) (note, the last two levels are not used for eddy events)
    }

    tMix = odtl->odtp->fmix * levelTaus[odtl->odtp->nLevels-1];

    Nm1 = odtl->odtp->nLevels - 1;
    Nm3 = odtl->odtp->nLevels - 3;

    //------------------- Set the parcel addresses (index array)

    pLoc.resize(odtl->ngrd);
    for(int i=0; i<odtl->ngrd; i++)
        pLoc[i] = i;

    //------------------- set the eddy event times at all levels

    setEddyEventTimes();

    for(int i=0; i<eTL.size(); i++)  //doldb
        cout << endl << "eTL: i, T, L: " << i << " " << eTL[i].first << " " << eTL[i].second;

}
///////////////////////////////////////////////////////////////////////////////
/** Sets arrays s.eTimes and s.eLevels which hold the eddy event times and the corresponding tree base level.
    First make an array of arrays for times at each level.
    These are from a Poisson process with the given rate at each level.
    Then collapse these into a single array.
    Then sort these times, along with the corresponding level array.
 */

void solver_hips::setEddyEventTimes(){

    double dmb;

    for(int i=0; i<levelRates.size()-2; i++) {
        double rate = levelRates[i];
        double time = 0.0;
        for(;;) {
            double r = odtl->rand->getRand();
            time += -log(r)/rate;
            if(time > odtl->odtp->tEnd)
                break;
            eTL.push_back(make_pair(time, i));
        }
    }

    //--------- now sort the list of times and the associated levels

    sort(eTL.begin(), eTL.end(), sortFunc);

}

///////////////////////////////////////////////////////////////////////////////
/** Function performs eddy events: parcel swaps.
    @param iLevel \input  level of the tree for the base of the swap.
    @param Qstart \output starting index for the Q-tree.
    @param Rstart \output starting index for the R-tree.
    @param nPswap \output number of parcels swapped.

    Randomly select a node on iLevel.
    Go down two levels and select nodes 0q and 1r, where q, r are randomly 0 or 1
    Find the starting index of the Q-tree and R-tree to swap and the number of parcels.
    Then swap the cells.

    For a 6 level tree: 0, 1, 2, 3, 4, 5:
    If iLevel = 1, then suppose i=1, 0q = 00 and 1r = 11:
    Then we are swaping 0100** with 0111** or (01|00|**) with (01|11|**)
       or i0qs with i1rs, where i = 01; 0q = 00; 1r = 11; and s = **

    We use bitwise shifts for easy powers of 2.
    The swap is done by adding or subtracting a value (shift),
        which should be equivalent to flipping the swapping the two 0q bits and 1r bits.
                                                                                                              Level
                                                                                                            ---------
                                                    *                                                           0
                                                 /     \
                                              /           \
                                           /                 \
                                        /                       \
                                     /                             \
                                  /                                   \
                               /                                         \
                            /                                               \
                           *                                                (*)  01|0000                        1
                          / \                                               / \
                        /     \                                           /     \
                      /         \                                       /         \
                    /             \                                   /             \
                  /                 \                               /                 \
                /                     \                           /                     \
               *                       *                         *                       *                      2
              / \                     / \                       / \                     / \
            /     \                 /     \                   /     \                 /     \
          /         \             /         \               /         \             /         \
         *           *           *           *            [*] 00|**    *           *          [*] 11|**         3
        / \         / \         / \         / \           / \         / \         / \         / \
       /   \       /   \       /   \       /   \         /   \       /   \       /   \       /   \
      *     *     *     *     *     *     *     *       *     *     *     *     *     *     *     *             4
     / \   / \   / \   / \   / \   / \   / \   / \     / \   / \   / \   / \   / \   / \   / \   / \
    00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15   16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31           5
                                                      ^^^^^^^^^^^                         ^^^^^^^^^^^

*/
void solver_hips::selectAndSwapTwoSubtrees(const int iLevel, int &Qstart, int &Rstart, int &nPswap){

    int iNode = odtl->rand->getRandInt((1 << iLevel)-1);  // base node for the swap
    int zero_q;
    int one_r;
    if(iLevel == Nm3) {                               // if swapping parcels (not nodes)
        zero_q = 0;                                       // it doesn't matter which of the two we choose
        one_r = 2;                                        // so choose the left to make simple mixing easier
    }
    else {
        zero_q = odtl->rand->getRandInt(1);          // 0q where q is 0 or 1
        one_r  = 2+odtl->rand->getRandInt(1);        // 1r where r is 0 or 1
    }

    Qstart = (zero_q << (Nm3-iLevel)) + (iNode << (Nm1-iLevel));  // starting index of Q parcels
    Rstart = (one_r  << (Nm3-iLevel)) + (iNode << (Nm1-iLevel));  // starting index of R parcels
    nPswap = 1 << (Nm3-iLevel);                                   // number of parcels that will be swapped

    int Qend = Qstart + nPswap;                      // inclusive indices are Qstart to Qend-1
    int Rend = Rstart + nPswap;                      // inclusive indices are Rstart to Rend-1

    vector<int> aa(pLoc.begin()+Qstart, pLoc.begin()+Qend);
    copy(pLoc.begin()+Rstart, pLoc.begin()+Rend, pLoc.begin()+Qstart);  // python: pLoc[Qstart:Qend]=pLoc[Rstart:Rend]
    copy(aa.begin(), aa.end(), pLoc.begin()+Rstart);                    // python: pLoc[Rstart:Rend]=aa

}

///////////////////////////////////////////////////////////////////////////////
/** The actual solver
 */

void solver_hips::calculateSolution() {

    odtl->io->writeDataFile("hips_init.dat", 0.0);

    //--------------------

    int QS, RS, nPs;
    double tPrev = 0.0;

    for(int i=0; i<eTL.size(); i++) {

        cout << endl << "eddy #, time, level: "
                     << i << "  " << eTL[i].first
                     << "  " << eTL[i].second;
        cout.flush();
        odtl->mimx->advanceOdt(tPrev, eTL[i].first);
        selectAndSwapTwoSubtrees(eTL[i].second, QS, RS, nPs);
        tPrev = eTL[i].first;

        if(i % odtl->odtp->modDump == 0)  {
            stringstream ss1; string s1;
            ss1.clear();  ss1 << setfill('0') << setw(4) << i; ss1 >> s1;
            odtl->io->writeDataFile("hips_eddy_"+s1+".dat", eTL[i].first);
        }

    }
    odtl->mimx->advanceOdt(tPrev, odtl->odtp->tEnd);

    //--------------------

    odtl->io->writeDataFile("hips_final.dat", 0.0);
}

///////////////////////////////////////////////////////////////////////////////
/** Destructor
 */
solver_hips::~solver_hips(){
    if(ed3)   delete ed3;
    if(eddl3) delete eddl3;
}
