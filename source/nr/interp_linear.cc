#include "interp_linear.h"

Int Base_interp::locate(const Doub x)
{
	Int ju,jm,jl;
	if (n < 2 || mm < 2 || mm > n) throw("locate size error");
	Bool ascnd=(xx[n-1] >= xx[0]);
	jl=0;
	ju=n-1;
	while (ju-jl > 1) {
		jm = (ju+jl) >> 1;
		if (x >= xx[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	cor = abs(jl-jsav) > dj ? 0 : 1;
	jsav = jl;
	return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
}
Int Base_interp::hunt(const Doub x)
{
	Int jl=jsav, jm, ju, inc=1;
	if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
	Bool ascnd=(xx[n-1] >= xx[0]);
	if (jl < 0 || jl > n-1) {
		jl=0;
		ju=n-1;
	} else {
		if (x >= xx[jl] == ascnd) {
			for (;;) {
				ju = jl + inc;
				if (ju >= n-1) { ju = n-1; break;}
				else if (x < xx[ju] == ascnd) break;
				else {
					jl = ju;
					inc += inc;
				}
			}
		} else {
			ju = jl;
			for (;;) {
				jl = jl - inc;
				if (jl <= 0) { jl = 0; break;}
				else if (x >= xx[jl] == ascnd) break;
				else {
					ju = jl;
					inc += inc;
				}
			}
		}
	}
	while (ju-jl > 1) {
		jm = (ju+jl) >> 1;
		if (x >= xx[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	cor = abs(jl-jsav) > dj ? 0 : 1;
	jsav = jl;
	return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
}

