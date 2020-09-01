/**
 * @file inputoutput.h
 * Header file for class inputoutput
 */

#ifndef INPUTOUTPUT_H
#define INPUTOUTPUT_H

#include <vector>
#include <string>
#include <ostream>
#include <fstream>
#include "yaml-cpp/yaml.h"

class odtline;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing inputoutput object
 *
 *  @author David O. Lignell
 */

class inputoutput {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        odtline                  *odtl;          ///< pointer to line object

        ostream                  *ostrm;         ///< Runtime: points to cout or to a file

        string                   caseName;       ///< input file directory
        string                   inputFileDir;   ///< input file directory
        string                   dataDir;        ///< data directory (output)

        YAML::Node               inputFile;      ///< yaml input file object base node
        YAML::Node               params;         ///< yaml sub node
        YAML::Node               sootParams;     ///< yaml sub node
        YAML::Node               streamProps;    ///< yaml sub node
        YAML::Node               initParams;     ///< yaml sub node
	YAML::Node               initParticleParams; ///< yaml sub node
        YAML::Node               dTimes;         ///< yaml sub node
        YAML::Node               bcCond;         ///< yaml sub node

        vector<double>           dumpTimes;      ///< vector of dump times
        int                      iNextDumpTime;  ///< index of next dump time
        bool                     LdoDump;        ///< flag for whether we are dumping a file

        ofstream                 gnufile;        ///< gnuplot script file

    //////////////////// MEMBER FUNCTIONS /////////////////

    void outputProperties(const string fname,
                          const double time);    ///< actually write the data file
    void dumpLineIfNeeded();                     ///< calls outputProperties for dumpTimes
    void writeDataFile(const string fnameRaw,
                       const double time);       ///< writes the gnuplot file and calls outputProperties
    void outputHeader();                         ///< output header info during odt solution
    void outputProgress();                       ///< output data going with the header info
    void outputKernHeader();                     ///< output header info during odt solution (Kernel events)
    void outputKernProgress();                   ///< output data going with the header info (Kernel events)

    void loadVarsFromRestartFile();

    void set_iNextDumpTime(double time);

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        inputoutput(const string p_inputFileDir, const int nShift);
        void init(odtline *p_odtl);
        ~inputoutput();

};


////////////////////////////////////////////////////////////////////////////////

#endif

