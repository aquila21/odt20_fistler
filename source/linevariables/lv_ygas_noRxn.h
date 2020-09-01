/**
 * @file lv_ygas_noRxn.h
 * Header file for class lv_ygas_noRxn
 */

#ifndef LV_YGAS_NORXN_H
#define LV_YGAS_NORXN_H

#include "lv_ygas.h"
#include <string>
#include <vector>

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child lv_ygas_Rsoot of parent lv_ygas object.
 *  Inherit everything, but add soot rates to the gas source term.
 *
 *  @author David O. Lignell
 */

class lv_ygas_noRxn : public lv_ygas {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    public:

    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void getRhsSrc(const int ipt=-1);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        lv_ygas_noRxn(){}
        lv_ygas_noRxn(odtline     *line,
                      const string s,
                      const bool   Lt,
                      const bool   Lo=true) : lv_ygas(line, s, Lt, Lo) {}

        virtual ~lv_ygas_noRxn(){}

};


////////////////////////////////////////////////////////////////////////////////

#endif
