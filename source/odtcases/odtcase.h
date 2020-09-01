/**
 * @file odtcase.h
 * Header file for class odtcase
 */

#ifndef ODTCASE_H
#define ODTCASE_H

class odtline;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing base odtcase object.
 *  Derived odtcase will be u,v,w,Yi,etc.
 *
 *  @author David O. Lignell
 */

class odtcase {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        odtline                       *odtl;                  ///< pointer to line object (parent)

    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setGasStateAtPt(const int &ipt){}
        virtual void setCaseSpecificVars(){}
        virtual void setCaseSpecificVars_cvode(const int &ipt){}
        virtual void enforceMassFractions();
        virtual void enforceSootMom();

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        odtcase(){}
        virtual void init(odtline *p_odtl) { odtl = p_odtl; }
        virtual ~odtcase(){}

};

////////////////////////////////////////////////////////////////////////////////

#endif
