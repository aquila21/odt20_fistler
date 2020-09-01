
#include "odtline.h"
#include "streams.h"
#include "processor.h"
#include "odtparam.h"
#include "micromixer.h"
#include "micromixer_flmlt.h"
#include "micromixer_hips.h"
#include "meshManager.h"
#include "eddy.h"
#include "solver.h"
#include "solver_flmlt.h"
#include "solver_hips.h"
#include "randomGenerator.h"
#include "particle.h"
#include "kernel.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"

#include <iostream>
#include <string>
#include <ctime>
#include <sstream>

using namespace std;
using namespace Cantera;

//////////////////////////////////////////////////////////////

processor proc;

//////////////////////////////////////////////////////////////

int main(int argc, char*argv[]) {


    if(argc<3) {
        cout << endl << "ERROR: code needs caseName and shift arguments" << endl;
        return 1;
    }
    string caseName= argv[1];            // example: temporalJet (../input/temporalJet, without the ../input/)

    int nShiftFileNumbers = 0;
    stringstream ss1;
    string       s1;
    ss1.clear(); ss1 << argv[2];
    ss1 >> nShiftFileNumbers;

    inputoutput io(caseName, nShiftFileNumbers);
    odtparam    odtp(&io);
    streams     strm;
    IdealGasMix gas("../data/"+caseName+"/input/"+odtp.chemMechFile);
    Transport   *tran = newTransportMgr("Mix", &gas);
    eddy        ed;
    meshManager mesher;
    solver      *solv;
    micromixer  *mimx;
    particle	part;
    kernel      kern;
    if(odtp.LisFlmlt) {
        solv = new solver_flmlt();
        mimx = new micromixer_flmlt();
    }
    else if(odtp.LisHips) {
        solv = new solver_hips();
        mimx = new micromixer_hips();
    }
    else {
        solv = new solver();
        mimx = new micromixer();
    }

    odtline odtl(NULL,  &odtp);
    odtline eddl(&odtl, &odtp);

    randomGenerator rand(odtp.seed);

    odtl.init(&io,  &mesher, &strm, &gas, tran, mimx, &ed, &eddl, solv, &part, &kern, &rand);
    eddl.init(NULL, NULL, NULL,  NULL, NULL, NULL,  NULL, NULL, NULL, NULL,  NULL, NULL, true);
    //
    //-------------------

    time_t mytimeStart, mytimeEnd;
    mytimeStart = time(0);
    *io.ostrm << endl << "#################################################################";
    *io.ostrm << endl << "#  Start Time = " << ctime(&mytimeStart);
    *io.ostrm << endl << "#################################################################";


    //-------------------

    odtl.solv->calculateSolution();

    //odtl.io->outputProperties("../data/init.dat", 0.0); //doldb
    //odtl.mimx->advanceOdt(0.0, odtl.odtp->tEnd);        //doldb

    //-------------------

    //delete mimx;
    //delete solv;

    mytimeEnd = time(0);
    *io.ostrm << endl << "#################################################################";
    *io.ostrm << endl << "#  Start Time = " << ctime(&mytimeStart);
    *io.ostrm << endl << "#  End Time   = " << ctime(&mytimeEnd);
    *io.ostrm << endl << "#################################################################";
    *io.ostrm << endl;

    return 0;


}
