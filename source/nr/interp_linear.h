#ifndef INTERP_LINEAR_H
#define INTERP_LINEAR_H

#include "nr3.h"

struct Base_interp
{
	Int n, mm, jsav, cor, dj;
	const Doub *xx, *yy;
	Base_interp(VecDoub_I &x, const Doub *y, Int m)
		: n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y) {
		dj = MAX(1,(int)pow((Doub)n,0.25));
	}

	Doub interp(Doub x) {
		Int jlo = cor ? hunt(x) : locate(x);
		return rawinterp(jlo,x);
	}

	Int locate(const Doub x);
	Int hunt(const Doub x);

	Doub virtual rawinterp(Int jlo, Doub x) = 0;

};

struct Linear_interp : Base_interp
{
	Linear_interp(VecDoub_I &xv, VecDoub_I &yv)
		: Base_interp(xv,&yv[0],2)  {}
	Doub rawinterp(Int j, Doub x) {
		if (xx[j]==xx[j+1]) return yy[j];
		else return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
	}
};

#endif
