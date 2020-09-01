/**
 * @file lv.h
 * Header file for class lv
 */

#ifndef LV_H
#define LV_H

#include <string>
#include <vector>

class odtline;

using namespace std;


////////////////////////////////////////////////////////////////////////////////

/** Class implementing base lv object.
 *  Derived lv will be u,v,w,Yi,etc.
 *
 *  @author David O. Lignell
 */

class lv {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        string                        var_name;               ///< name of variable
        vector<double>                d;                      ///< the data
        bool                          L_transported;          ///< flag true if var is transported
        bool                          L_output;               ///< flag true if included in output
        bool                          LagSrc;                 ///< flag to lag source term in implicit solve (initially put in for enthalpy radiation)

        odtline                       *odtl;                  ///< pointer to line object (parent)

        vector<double>                rhsSrc;                 ///< the data
        vector<double>                rhsMix;                 ///< the data

        vector<double>                flux;

        //---------- for flmlt interface (inherited)

        double                        chi0;                   ///< scalar dissipation rate


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void   setVar(const int ipt=-1){}

        virtual void   merge2cells(const int    imrg,
                                   const double m2,
                                   const double m1,
                                   const bool   LconstVolume=false);

        virtual void   splitCell(const int            isplt,
                                 const int            nsplt,
                                 const vector<double> &cellFaces);

        virtual void   getRhsSrc(const int ipt=-1){if(!L_transported) return;}
        virtual void   getRhsMix(const vector<double> &gf,
                                 const vector<double> &dxc){if(!L_transported) return;}

        virtual void   interpVarToFacesHarmonic(const vector<double> &cvar, vector<double> &fvar);
        virtual double linearInterpToFace(const int &iface, const vector<double> &vec);
        virtual void   setLvFromRegion(const int i1, const int i2);
        virtual void   resize();

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv(){}
        lv(odtline      *line,
           const string s,
           const bool   Lt,
           const bool   Lo=true);

        virtual ~lv(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
