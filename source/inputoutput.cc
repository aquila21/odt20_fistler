
#include "inputoutput.h"
#include "odtline.h"
#include "processor.h"
#include <sys/stat.h>             // for mkdir
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>

extern processor proc;

///////////////////////////////////////////////////////////////////////////////
/** inputoutput initialization function
 *
 * @param p_odtl  \input set odtline pointer with.
 */

void inputoutput::init(odtline *p_odtl) {
    odtl    = p_odtl;
    LdoDump = false;
}

///////////////////////////////////////////////////////////////////////////////
/** inputoutput constructor function
 *
 * @param p_caseName \input directory with input files
 * @param p_nShift   \input shift the file numbers by this amount (used for multiple sets of parallel runs).
 */

inputoutput::inputoutput(const string p_caseName, const int p_nShift){

    caseName     = p_caseName;
    inputFileDir = "../data/"+caseName+"/input/";

    inputFile   = YAML::LoadFile(inputFileDir+"odt_input.yaml");     ///< use these "nodes" to access parameters as needed
    params      = inputFile["params"];
    sootParams  = inputFile["sootParams"];
    streamProps = inputFile["streamProps"];
    initParams  = inputFile["initParams"];
	initParticleParams = inputFile["initParticleParams"];
    dTimes      = inputFile["dumpTimes"];
    bcCond      = inputFile["bcCond"];

	int N = params["tEnd"].as<double>()/dTimes[0].as<double>();

	for(int i=1; i<N+1; i++)
        dumpTimes.push_back(i*dTimes[0].as<double>());
    dumpTimes.push_back(1.0E200);                       ///< add an extra "infinity" for easy handling of the last time
    iNextDumpTime = 0;

//    for(int i=0; i<dTimes.size(); i++)
//        dumpTimes.push_back(dTimes[i].as<double>());
//    dumpTimes.push_back(1.0E200);                       ///< add an extra "infinity" for easy handling of the last time
//    iNextDumpTime = 0;

    //----------- set the data directory and runtime file

    string fname;
    stringstream ss1;
    string       s1;

    ss1.clear(); ss1 << setfill('0') << setw(5) << proc.myid + p_nShift;
    s1 = ss1.str();
    dataDir = "../data/"+caseName+"/data/data_" + s1 + "/";   // e.g., "../data_00001", etc.

    int iflag = mkdir(dataDir.c_str(), 0755);
    if(iflag != 0) {
        cout << "\n********** Error, process " << proc.myid << "failed to create "
             << dataDir << ", or it was already there" << endl;
        exit(0);
    }

    fname = "../data/"+caseName+"/runtime/runtime_" + s1;
    ostrm = new ofstream(fname.c_str());

    //----------- set gnuplot file

    fname = dataDir + "plot_odt.gnu";
    gnufile.open(fname.c_str());
    if(!gnufile) {
        cout << endl << "ERROR OPENING FILE: " << dataDir+"plot_odt.gnu" << endl;
        exit(0);
    }

}

///////////////////////////////////////////////////////////////////////////////
/** Destructor function
 */

inputoutput::~inputoutput() {

    delete ostrm;
    gnufile.close();

}

///////////////////////////////////////////////////////////////////////////////
/** Writes a data file of the line properties (in order)
 *
 * @param fname \input name of the file to write including path
 *  @param time \input time of the output
 */

void inputoutput::outputProperties(const string fname, const double time) {

    string       s1;
    stringstream ss1;

    //--------------------------

    for(int i=0; i<odtl->v.size(); i++)
        odtl->v.at(i)->setVar();             // make sure the variable is up to date

    //--------------------------

    ofstream ofile(fname.c_str());
    if(!ofile) {
        *ostrm << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;
        exit(0);
    }

    *ostrm << endl << "# Writing outputfile: " << fname;
    ostrm->flush();

    //--------------------------

    ofile << "# time = "   << time;

    ofile << "\n# Grid points = "   << odtl->ngrd;

    if(!odtl->odtp->LisHips)
        ofile << "\n# Domain Size = " << odtl->Ldomain();

    ofile << "\n# Pressure (Pa) = " << odtl->odtp->pres << endl;

    ofile << "#";
    for(int i=0,j=1; i<odtl->v.size(); i++){
        if(odtl->v.at(i)->L_output)
            ofile << setw(14) << j++ << "_" << odtl->v.at(i)->var_name;
    }

    int iploc;
    ofile << scientific;
    ofile << setprecision(10);
    for(int i=0; i<odtl->ngrd; i++) {
        iploc = (odtl->odtp->LisHips) ? odtl->solv->pLoc[i] : i;    // HiPS uses an index array
        ofile << endl;
        for(int k=0; k<odtl->v.size(); k++)
            if(odtl->v.at(k)->L_output)
                ofile << setw(19) << odtl->v.at(k)->d.at(iploc);
    }

    ofile.close();

}

///////////////////////////////////////////////////////////////////////////////
/** Set iNextDumpTime from time. Used for restarts.
 * @param time \input time to use to set iNextDumpTime.
 */

void inputoutput::set_iNextDumpTime(double time) {

    for(int i=0; i<dumpTimes.size(); i++)
        if(dumpTimes[i] >= time) {
            iNextDumpTime = i;
            break;
        }
}

///////////////////////////////////////////////////////////////////////////////
/** Dumps a line, sets flag, increments next dump
 */

void inputoutput::dumpLineIfNeeded(){

    if(!LdoDump) return;

    stringstream ss;
    ss << setfill('0') << setw(5) << iNextDumpTime;
    string fname = "dmp_" + ss.str() + ".dat";
	//outputProperties(fname, dumpTimes.at(iNextDumpTime));

	writeDataFile(fname, dumpTimes.at(iNextDumpTime));
	if(odtl->part->particleON) 
		odtl->part->writeData("particle_"+fname, dumpTimes.at(iNextDumpTime));

    iNextDumpTime++;
    LdoDump = false;
}

///////////////////////////////////////////////////////////////////////////////
/** Dumps a line, sets flag, increments next dump
 *  @param fnameRaw \input file name without the path (just the name).
 *  @param time \input time of the output
 */

void inputoutput::writeDataFile(const string fnameRaw, const double time) {

    string fname =dataDir + fnameRaw;
    outputProperties(fname, time);
    gnufile << "plot '" << fnameRaw << "' us 1:5; pause -1;" << endl;

}

///////////////////////////////////////////////////////////////////////////////
/**Output title of properties displayed to screen. */

void inputoutput::outputHeader() {

    *ostrm << endl << "#--------------------------------------------------"
                   << "--------------------------------------------------------------------";
    *ostrm << endl;
    *ostrm << setw(5) << "# EE,"
        << setw(12) << "time,"
        << setw(12) << "t-t0,"
        << setw(10) << "nEtry,"
        << setw(6)  << "ngrd,"
        << setw(12) << "edSize,"
        << setw(12) << "edPos,"
        << setw(12) << "edPa,"
        << setw(12) << "nEposs,"
        << setw(12) << "PaAvg,"
        << setw(12) << "invTauEddy"
        ;

    *ostrm << endl << "#--------------------------------------------------"
                   << "--------------------------------------------------------------------";
}

///////////////////////////////////////////////////////////////////////////////
/**Outputs the data corresponding to outputHeader.
 * After a given number of accepted eddies, output this info.
 *
 */

void inputoutput::outputProgress() {

    double dmb = 0.5*(odtl->ed->leftEdge + odtl->ed->rightEdge);
    if(dmb > odtl->posf->d.at(odtl->ngrd))
        dmb = dmb-odtl->Ldomain();

    *ostrm << scientific << setprecision(3) << endl;
    *ostrm << setw(5)  << odtl->solv->neddies                    //  1: EE
           << setw(12) << odtl->solv->time                       //  2: time
           << setw(12) << odtl->solv->time-odtl->solv->t0        //  3: t-t0
           << setw(10) << odtl->solv->iEtrials                   //  4: nEtry
           << setw(6)  << odtl->ngrd                             //  5: ngrd
           << setw(12) << odtl->ed->eddySize                     //  6: edSize
           << setw(12) << dmb                                    //  7: edPos
           << setw(12) << odtl->ed->Pa                           //  8: edPa
           << setw(12) << odtl->solv->nPaSumC                    //  9: nEposs
           << setw(12) << odtl->solv->PaSumC/odtl->solv->nPaSumC // 10: PaAvg
           << setw(12) << odtl->ed->invTauEddy                   // 11: invTauEddy
        ;
    ostrm->flush();
}

///////////////////////////////////////////////////////////////////////////////
/** Restart
 *  The number of columns in the restart file should match the number and order of linevariables
 *      that are output to a data file. (That is, the order of the linevariables in odtline).
 */

void inputoutput::loadVarsFromRestartFile() {

    string fname;
    stringstream ss1;
    string       s1;

    for(int k=0; k<odtl->v.size(); k++) {
        if(odtl->v[k]->L_transported && !odtl->v[k]->L_output) {
            cout << endl << "ERROR: to restart, all transported variables need to be in the restart file" << endl;
            exit(0);
        }
    }

    if(odtl->odtp->rstType == "multiple") {
        ss1.clear(); ss1 << setfill('0') << setw(5) << proc.myid;
        fname = inputFileDir + "restart/restart_" + ss1.str() + ".dat";
    }
    else
        fname = inputFileDir + "/restart.dat";

    ifstream ifile(fname.c_str());
    if(!ifile) {
        cout << endl << "ERROR: reading restart file " << fname << endl;
        exit(0);
    }

    //------------- Get file header information

    getline(ifile, s1);                        // read line "# time = 1.1" (this is the restart time
    ss1.clear();
    ss1.str(s1);
    ss1 >> s1 >> s1 >> s1 >> odtl->odtp->trst;

    getline(ifile, s1);                        // read line "# Grid points = 100"
    ss1.clear();
    ss1.str(s1);
    ss1 >> s1 >> s1 >> s1 >> s1 >> odtl->ngrd;
    odtl->ngrdf = odtl->ngrd+1;

    getline(ifile, s1);                        // read line "# Domain Size = 2" (don't use)
    getline(ifile, s1);                        // read line "# Pressure (Pa) = 101325
    getline(ifile, s1);                        // read line "# column headers

    //------------- Get file data columns

    for(int k=0; k<odtl->v.size(); k++)
       odtl->v[k]->d.resize(odtl->ngrd);
    odtl->posf->d.resize(odtl->ngrdf);

    for(int i=0; i<odtl->ngrd; i++) {
        for(int k=0; k<odtl->v.size(); k++) {
            if(!odtl->v[k]->L_output)
                continue;
            ifile >> odtl->v[k]->d[i];
        }
    }

    odtl->posf->d[odtl->ngrd] = odtl->posf->d[0] + odtl->odtp->domainLength;

    //------------- Set the variables

    for(int k=0; k<odtl->v.size(); k++)
        odtl->v[k]->setVar();

}

///////////////////////////////////////////////////////////////////////////////
/**Output title of properties displayed to screen. */

void inputoutput::outputKernHeader() {

    *ostrm << endl << "#--------------------------------------------------"
                   << "--------------------------------------------------------------------";
    *ostrm << endl;
    *ostrm << setw(5) << "# K,"
        << setw(12) << "time,"
        << setw(12) << "injTime,"
        << setw(12) << "kernSize,"
        << setw(12)  << "kernLE,"
        << setw(12) << "injTKE,"
        << setw(12) << "deltaE,"
        ;

    *ostrm << endl << "#--------------------------------------------------"
                   << "--------------------------------------------------------------------";
}

///////////////////////////////////////////////////////////////////////////////
/**Outputs the data corresponding to outputHeader.
 * After a given number of accepted eddies, output this info.
 *
 */
void inputoutput::outputKernProgress() {

    double dmb = 0.0;
    if(odtl->kern->rightEdge > odtl->kern->leftEdge){
        dmb = 0.5*(odtl->kern->leftEdge + odtl->kern->rightEdge);
    }
    else if(abs(odtl->kern->rightEdge) > odtl->kern->leftEdge){
        dmb = 0.5*(odtl->kern->rightEdge + odtl->kern->leftEdge - odtl->Ldomain());
    }
    else{
        dmb = 0.5*(odtl->kern->rightEdge + odtl->kern->leftEdge + odtl->Ldomain());
    }
  
    *ostrm << scientific << setprecision(3) << endl;
    *ostrm << setw(5)  << odtl->kern->nkernels                   //  1: K
           << setw(12) << odtl->solv->time                       //  2: time
           << setw(12) << odtl->kern->kernelTime                 //  3: inj. time
           << setw(12) << odtl->kern->kernelSize                 //  4: kernel size
           << setw(12) << odtl->kern->leftEdge                   //  5: kernel pos
           << setw(12) << odtl->kern->E                          //  6: inj TKE
           << setw(12) << odtl->kern->DeltaE                     //  7: actual energy diff
        ;
    ostrm->flush();
}


