#include "odtparam.h"
#include "odtline.h"

///////////////////////////////////////////////////////////////////////////////
/** odtparam initialization function
 * @param p_odtl  \input set odtline pointer with.
 */

void odtparam::init(odtline *p_odtl) {
    odtl = p_odtl;
}

///////////////////////////////////////////////////////////////////////////////
/** odtparam constructor function
 *
 * @param p_io  \input set inputoutput pointer with.
 */

odtparam::odtparam(inputoutput *p_io) {

    io = p_io;

    seed           = io->params["seed"]           ? io->params["seed"].as<int>()             : errMsg<int>("seed");
    tEnd           = io->params["tEnd"]           ? io->params["tEnd"].as<double>()          : errMsg<double>("tEnd");
    domainLength   = io->params["domainLength"]   ? io->params["domainLength"].as<double>()  : errMsg<double>("domainLength");
    ngrd0          = io->params["ngrd0"]          ? io->params["ngrd0"].as<int>()            : 1000;     //errMsg<int>("ngrd0");
    rho0           = io->params["rho0"]           ? io->params["rho0"].as<double>()          : 1.0;      //errMsg<double>("rho0");
    kvisc0         = io->params["kvisc0"]         ? io->params["kvisc0"].as<double>()        : 0.001694; //errMsg<double>("kvisc0");
    sdiff0         = io->params["sdiff0"]         ? io->params["sdiff0"].as<double>()        : 0.001694; //errMsg<double>("sdiff0");
    dPdx           = io->params["dPdx"]           ? io->params["dPdx"].as<double>()          : 0.0;
    pres           = io->params["pres"]           ? io->params["pres"].as<double>()          : 101325.0;
    chemMechFile   = io->params["chemMechFile"]   ? io->params["chemMechFile"].as<string>()  : errMsg<string>("chemMechFile");
    probType       = io->params["probType"]       ? io->params["probType"].as<string>()      : errMsg<string>("probType");

    Z_param        = io->params["Z_param"]        ? io->params["Z_param"].as<double>()       : 400.0;    //errMsg<double>("Z_param");
    A_param        = io->params["A_param"]        ? io->params["A_param"].as<double>()       : 0.666667; //errMsg<double>("A_param");
    C_param        = io->params["C_param"]        ? io->params["C_param"].as<double>()       : 5.0;      //errMsg<double>("C_param");
    LES_type       = io->params["LES_type"]       ? io->params["LES_type"].as<string>()      : "NONE";   //errMsg<string>("LES_type");
    Z_LES          = io->params["Z_LES"]          ? io->params["Z_LES"].as<double>()         : 0.0;
    x0virtual      = io->params["x0virtual"]      ? io->params["x0virtual"].as<double>()     : 0.0;
    diffCFL        = io->params["diffCFL"]        ? io->params["diffCFL"].as<double>()       : errMsg<double>("diffCFL");
    cvode_atol     = io->params["cvode_atol"]     ? io->params["cvode_atol"].as<double>()    : 1.0E-10;
    cvode_rtol     = io->params["cvode_rtol"]     ? io->params["cvode_rtol"].as<double>()    : 1.0E-4;

    radType        = io->params["radType"]        ? io->params["radType"].as<string>()       : "NONE";
    Lbuoyant       = io->params["Lbuoyant"]       ? io->params["Lbuoyant"].as<bool>()        : false;
    LPeEddy        = io->params["LPeEddy"]        ? io->params["LPeEddy"].as<bool>()         : false;
    g              = io->params["g"]              ? io->params["g"].as<double>()             : -9.81;
    LdoDL          = io->params["LdoDL"]          ? io->params["LdoDL"].as<bool>()           : false;
    Lsolver        = io->params["Lsolver"]        ? io->params["Lsolver"].as<string>()       : errMsg<string>("Lsolver");
    Lperiodic      = io->params["Lperiodic"]      ? io->params["Lperiodic"].as<bool>()       : false;
    Lspatial       = io->params["Lspatial"]       ? io->params["Lspatial"].as<bool>()        : false;
    Llem           = io->params["Llem"]           ? io->params["Llem"].as<bool>()            : false;
    LisFlmlt       = io->params["LisFlmlt"]       ? io->params["LisFlmlt"].as<bool>()        : false;

    bcType         = io->params["bcType"]         ? io->params["bcType"].as<string>()        : errMsg<string>("bcType");
    cCoord         = io->params["cCoord"]         ? io->params["cCoord"].as<int>()           : 1;
    xDomainCenter  = io->params["xDomainCenter"]  ? io->params["xDomainCenter"].as<double>() : 0.0;


    gDens          = io->params["gDens"]          ? io->params["gDens"].as<double>()         : 30;  //errMsg<double>("gDens");
    dxmin          = io->params["dxmin"]          ? io->params["dxmin"].as<double>()         : 0.0; //errMsg<double>("dxmin");
    dxmax          = io->params["dxmax"]          ? io->params["dxmax"].as<double>()         : 0.2;

    Pmax           = io->params["Pmax"]           ? io->params["Pmax"].as<double>()          : 0.4;    //errMsg<double>("Pmax");
    Pav            = io->params["Pav"]            ? io->params["Pav"].as<double>()           : 0.02;   //errMsg<double>("Pav");
    dtfac          = io->params["dtfac"]          ? io->params["dtfac"].as<double>()         : 2.0;    //errMsg<double>("dtfac");
    nDtSmeanWait   = io->params["nDtSmeanWait"]   ? io->params["nDtSmeanWait"].as<int>()     : 100000; //errMsg<int>("nDtSmeanWait");
    eddyMinCells   = io->params["eddyMinCells"]   ? io->params["eddyMinCells"].as<int>()     : 3;      //errMsg<int>("eddyMinCells");
    DAtimeFac      = io->params["DAtimeFac"]      ? io->params["DAtimeFac"].as<double>()     : 10.0;   //errMsg<double>("DAtimeFac");
    tdfac          = io->params["tdfac"]          ? io->params["tdfac"].as<double>()         : 1.0;
    sLastDA        = io->params["sLastDA"]        ? io->params["sLastDA"].as<int>()          : 100;    //errMsg<int>("sLastDA");
    Lp             = io->params["Lp"]             ? io->params["Lp"].as<double>()            : 0.01;   //errMsg<double>("Lp");
    Lmax           = io->params["Lmax"]           ? io->params["Lmax"].as<double>()          : 1.0;    //errMsg<double>("Lmax");
    Lmin           = io->params["Lmin"]           ? io->params["Lmin"].as<double>()          : dxmin*eddyMinCells; //errMsg<double>("Lmin");

    modDump        = io->params["modDump"]        ? io->params["modDump"].as<int>()          : 1000000; //errMsg<int>("modDump");
    modDisp        = io->params["modDisp"]        ? io->params["modDisp"].as<int>()          : 1;

    LmultiPhase    = io->params["LmultiPhase"]    ? io->params["LmultiPhase"].as<bool>()     : false;
    eSurfTens      = io->params["eSurfTens"]      ? io->params["eSurfTens"].as<double>()     : 0.0;

    Lrestart       = io->params["Lrestart"]       ? io->params["Lrestart"].as<bool>()        : false;
    rstType        = io->params["rstType"]        ? io->params["rstType"].as<string>()       : "single";    // "single" or "multiple"
    trst = 0.0; // (dont read this in, it comes from the restart file

    umin_spatial   = io->params["umin_spatial"]   ? io->params["umin_spatial"].as<double>()  : 0.5;

    // HIPS variables ---------------------

    LisHips        = io->params["LisHips"]        ? io->params["LisHips"].as<bool>()         : false;
    nLevels        = io->params["nLevels"]        ? io->params["nLevels"].as<int>()          : 0;
    Afac           = io->params["Afac"]           ? io->params["Afac"].as<double>()          : 0.5;
    tau0           = io->params["tau0"]           ? io->params["tau0"].as<double>()          : 0.0;
    fmix           = io->params["fmix"]           ? io->params["fmix"].as<double>()          : 0.0;

    // Soot variables ---------------------

    Lsoot             = io->params["Lsoot"]                ? io->params["Lsoot"].as<bool>()                     : false;
    nsvar             = io->sootParams["nsvar"]            ? io->sootParams["nsvar"].as<int>()                  : 0;
    b_coag            = io->sootParams["b_coag"]           ? io->sootParams["b_coag"].as<double>()              : 1.0;                  ///< coagulation constant                  
    nsvar_v           = io->sootParams["nsvar_v"]          ? io->sootParams["nsvar_v"].as<int>()                : 0;
    nsvar_s           = io->sootParams["nsvar_s"]          ? io->sootParams["nsvar_s"].as<int>()                : 0;
    rho_soot          = io->sootParams["rho_s"]            ? io->sootParams["rho_s"].as<double>()               : 2000.0;              ///< solid soot density, kg/m3
    Cmin              = io->sootParams["Cmin"]             ? io->sootParams["Cmin"].as<int>()                   : 10;                  ///< minimum number of carbon atoms in a soot particle
    nsvar             = io->sootParams["nsvar"]            ? io->sootParams["nsvar"].as<int>()                  : 0;                 ///< number of soot variables (i.e. moments, bins)           
    PSD_method        = io->sootParams["PSD_method"]       ? io->sootParams["PSD_method"].as<string>()          : "NONE";            
    nucleation_mech   = io->sootParams["nucleation_mech"]  ? io->sootParams["nucleation_mech"].as<string>()     : "NONE";    
    growth_mech       = io->sootParams["growth_mech"]      ? io->sootParams["growth_mech"].as<string>()         : "NONE";        
    oxidation_mech    = io->sootParams["oxidation_mech"]   ? io->sootParams["oxidation_mech"].as<string>()      : "NONE";        
    coagulation_mech  = io->sootParams["coagulation_mech"] ? io->sootParams["coagulation_mech"].as<string>()    : "NONE";     

    //----------forced HIT parameter-----------
    kernEv         = io->params["kernEv"]    ? io->params["kernEv"].as<bool>()                : false;
    Prod           = io->params["Prod"]      ? io->params["Prod"].as<double>()                : 0.0;
    T11            = io->params["T11"]       ? io->params["T11"].as<double>()                 : 0.0;
    L11            = io->params["L11"]       ? io->params["L11"].as<double>()                 : 0.0;
    //---------- HST parameter-----------

    p3             = io->params["p3"]        ? io->params["p3"].as<double>()                  : 0.0;
 
    //----------Particle parameters (MF) -----------
    partCoupl      = io->params["partCoupl"]       ? io->params["partCoupl"].as<int>()        : 0;
    tparticleON    = io->params["tparticleON"]     ? io->params["tparticleON"].as<double>()   : 0.0;
    betaP          = io->params["betaP"]           ? io->params["betaP"].as<double>()         : 1.0;
    Nparticle      = io->params["Nparticle"]       ? io->params["Nparticle"].as<int>()        : 0;
    Dparticle      = io->params["Dparticle"]       ? io->params["Dparticle"].as<double>()     : 0.0;
    DensiParticle  = io->params["DensiParticle"]   ? io->params["DensiParticle"].as<double>() : 0.0;
    typeP          = io->params["typeP"]           ? io->params["typeP"].as<int>()            : 1;
    Gx             = io->params["Gx"]              ? io->params["Gx"].as<double>()            : 0.0;
    Gy             = io->params["Gy"]              ? io->params["Gy"].as<double>()            : 0.0;
    Gz             = io->params["Gz"]              ? io->params["Gz"].as<double>()            : 0.0;
    typeI          = io->params["typeI"]           ? io->params["typeI"].as<bool>()           : true;
    typeC          = io->params["typeC"]           ? io->params["typeC"].as<bool>()           : false;
    crmax          = io->params["crmax"]           ? io->params["crmax"].as<double>()         : 10.0;
    BCparticle     = io->params["BCparticle"]      ? io->params["BCparticle"].as<int>()       : 0.0;
    alphaP         = io->params["alphaP"]          ? io->params["alphaP"].as<double>()        : 1.0;

    nStat          = io->params["nStat"]           ? io->params["nStat"].as<int>()            : 2.0;
    ngrdStat       = io->params["ngrdStat"]        ? io->params["ngrdStat"].as<int>()         : 200.0;

    // Dirichlet velocity BCs (no-slip)
    uBClo       = io->bcCond["uBClo"]       ? io->bcCond["uBClo"].as<double>()       : 0.0;
    uBChi       = io->bcCond["uBChi"]       ? io->bcCond["uBChi"].as<double>()       : 0.0;
    vBClo       = io->bcCond["vBClo"]       ? io->bcCond["vBClo"].as<double>()       : 0.0;
    vBChi       = io->bcCond["vBChi"]       ? io->bcCond["vBChi"].as<double>()       : 0.0;
    wBClo       = io->bcCond["wBClo"]       ? io->bcCond["wBClo"].as<double>()       : 0.0;
    wBChi       = io->bcCond["wBChi"]       ? io->bcCond["wBChi"].as<double>()       : 0.0;

    // Dirichlet scalar BCs
    sBClo       = io->bcCond["sBClo"]       ? io->bcCond["sBClo"].as<double>()       : 0.0; //errMsg<double>("sBClo");
    sBChi       = io->bcCond["sBChi"]       ? io->bcCond["sBChi"].as<double>()       : 0.0; //errMsg<double>("sBChi");

    // Dirichlet temperature BCs
    hWallBCtype = io->bcCond["hWallBCtype"] ? io->bcCond["hWallBCtype"].as<string>() : "ADIABATIC";
    if(hWallBCtype == "ISOTHERMAL") {
        TBClo   = io->bcCond["TBClo"]       ? io->bcCond["TBClo"].as<double>()       : errMsg<double>("TBClo");
        TBChi   = io->bcCond["TBChi"]       ? io->bcCond["TBChi"].as<double>()       : errMsg<double>("TBChi");
    }
    if(radType != "NONE") {
        TBClo   = io->bcCond["TBClo"]       ? io->bcCond["TBClo"].as<double>()       : errMsg<double>("TBClo");
        TBChi   = io->bcCond["TBChi"]       ? io->bcCond["TBChi"].as<double>()       : errMsg<double>("TBChi");
    }

    //--------------------- un-normalize

    dxmin *= domainLength;
    dxmax *= domainLength;
    Lmax  *= domainLength;
    Lmin  *= domainLength;
    Lp    *= domainLength;

    //--------------------- sanity checks

    if( (cCoord == 2) && (xDomainCenter != 0.0) && ( abs(xDomainCenter) < 0.5*domainLength) ) {
        cout << endl << "ERROR: for cylindrical, set xDomainCenter to be 0.0 (e.g, a pipe) or such as gives an annulus";
        cout << endl << "       That is, xDomainCenter should give a domain, that when rotated about the origin, doesnt overlap itself";
        exit(0);
    }
    if( (cCoord == 3) && (xDomainCenter != 0.0) )  {
        cout << endl << "ERROR: spherical case requires xDomainCenter to be zero";
        exit(0);
    }

}

