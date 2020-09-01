/**
 * @file lv_pos.h
 * Header file for class lv_pos
 */

#ifndef LV_POS_H
#define LV_POS_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_pos of parent lv object.
 *
 *  @author David O. Lignell
 */

class lv_pos : public lv {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setVar(const int ipt=-1);

        virtual void   merge2cells(const int    imrg,
                                   const double m2,
                                   const double m1,
                                   const bool   LconstVolume=false);

        virtual void splitCell(const int            isplt,
                               const int            nsplt,
                               const vector<double> &cellFaces);

        virtual void setLvFromRegion(const int i1, const int i2);

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_pos(){}
        lv_pos(odtline      *line,
               const string s,
               const bool   Lt,
               const bool   Lo=true);

        virtual ~lv_pos(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
