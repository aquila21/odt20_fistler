/**
 *  @file pdgen.h
 *	Header file for class pda
 */

#ifndef PDA_H
#define PDA_H

#include<vector>

using namespace std;

class pda {

    public:

        void pdAlg(int nm, int np,std::vector<double> &mu, std::vector<double> &wts, std::vector<double> &absc ) ;

        void imtql2(int                               nm, //nmd
                    int                               n,  //np
                    std::vector<double>               &d,  //absc
                    std::vector<double>               &e,  //b
                    std::vector<std::vector<double> > &z,  //P
                    int                               &ierr);

        double pythag(double a, double b);

        pda(){}
};

#endif
