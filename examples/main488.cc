
// Authors: Torbjorn Sjostrand <torbjorn.sjostrand@fysik.lu.se>
//          Chloé Gaudu <gaudu@uni-wuppertal.de>

// Keywords: switch beam; switch collision energy; reuse MPI initialization

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  vector<int> idA_list = {2212, -2212, 2112, -2112, 111, 211, -211, 311, 321, -321, 130, 310, 3122, 3212, 3222, 3112, 3322};
  vector<int> idB_list = {2212, 1000060120, 1000070140, 1000080160, 1000180400};
  double idA = idA_list[0];
  double idB = idB_list[3];
  double pA_max = 1e12;
  double pA_min = 1e2; 
  int nEvent = 10;
  int nPts = 11;

  const bool doLog = true;
  std::string out = "main488_"+std::to_string(nEvent);
  ofstream logBuf;
  std::streambuf* oldCout;
  if(doLog) {
    oldCout = std::cout.rdbuf(logBuf.rdbuf());
    logBuf.open((out == "" ? "pythia.log" : out + ".log"));
  }

  Pythia pythia;

  //pythia.readString("HeavyIon:mode = 2");

  pythia.readString("Beams:allowIDASwitch = on");
  pythia.settings.mode("Beams:idA", idA);
  pythia.settings.mode("Beams:idB", idB);
  
  double mA = pythia.particleData.m0(idA);
  double mB = pythia.particleData.m0(2212);
  double eA_max = sqrt(pow2(pA_max) + pow2(pythia.particleData.m0(3322))); // Xi0 heaviest 
  double eA_min = sqrt(pow2(pA_min) + pow2(pythia.particleData.m0(111))); // pi0 lightest
  double eCM_max = sqrt(2.* eA_max * mB);
  double eCM_min = sqrt(2.* eA_min * mB);

  pythia.readString("Beams:allowVariableEnergy = on");
  pythia.settings.parm("HeavyIon:varECMMin", eCM_min);
  pythia.settings.parm("HeavyIon:varECMMax", eCM_max);
  pythia.settings.parm("HeavyIon:varECMSigFitNPts", nPts);  
  pythia.readString("Beams:frameType = 1");
  pythia.settings.parm("Beams:eCM", eCM_max);

  pythia.readString("SoftQCD:all = on");

  // If you set reuseInit = 2, the old files will be replaced overwritten they already exist.
  pythia.readString("MultipartonInteractions:reuseInit = 3");
  pythia.readString("MultipartonInteractions:initFile = main488.mpi");
  pythia.readString("HeavyIon:SasdMpiReuseInit = 3");
  pythia.readString("HeavyIon:SasdMpiInitFile = main488.sasd.mpi");
  pythia.readString("HeavyIon:SigFitReuseInit = 3");
  pythia.readString("HeavyIon:SigFitInitFile = main488.sigfit");

  if (!pythia.init()) {
    cout << " Pythia failed to initialize." << endl;
    return -1;
  }

  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    
    double idA_switch = idA_list[iEvent%17];
    double idB_switch = idB_list[iEvent%5];
    pythia.setBeamIDs(idA_switch, idB_switch);

    double pA_switch = pA_min + pythia.rndm.flat() * (pA_max - pA_min);
    double mA_switch = pythia.particleData.m0(idA_switch);
    double eA_switch = sqrt(pow2(pA_switch) + pow2(mA_switch));
    double eCM_switch = sqrt(2.* eA_switch * mB);
    pythia.setKinematics(eCM_switch);

    /*std::cout << "Event: " << iEvent << std::endl;
    pythia.event.list();*/
    
    pythia.next();
  }

  pythia.stat();

  if (doLog) std::cout.rdbuf(oldCout);

  return 0;
}