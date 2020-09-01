/**
 * @file lv_posf.h
 * Header file for class lv_posf
 */

#ifndef LV_POSF_H
#define LV_POSF_H

#include "lv.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_posf of parent lv object.
 *
 *  @author David O. Lignell
 */

class lv_posf : public lv {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void   merge2cells(const int    imrg,
                                   const double m2,
                                   const double m1,
                                   const bool   LconstVolume=false);

        virtual void splitCell(const int            isplt,
                               const int            nsplt,
                               const vector<double> &cellFaces);

        virtual void setLvFromRegion(const int i1, const int i2);
        virtual void resize();

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_posf(){}
        lv_posf(odtline      *line,
                const string s,
                const bool   Lt,
                const bool   Lo=true);

        virtual ~lv_posf(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
