
// Authors: Chloé Gaudu <gaudu@uni-wuppertal.de>

// Keywords: switch beam; switch collision energy; reuse MPI initialization

// Attempt to replicate the Corsika 8 Interaction.inl file version from the
// 620-upgrade-to-pythia-8-310 branch on 29 July 2024 to test out the reuse 
// files.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  vector<int> idA_list = {2212, 211, -211, 1000070140};
  vector<int> idB_list = {1000070140, 1000080160, 1000180400}; 
  int nEvent = 1;

  const bool doLog = true;
  std::string out = "test488_c8_replica_"+std::to_string(nEvent);
  ofstream logBuf;
  std::streambuf* oldCout;
  if(doLog) {
    oldCout = std::cout.rdbuf(logBuf.rdbuf());
    logBuf.open((out == "" ? "pythia.log" : out + ".log"));
  }

  Pythia pythia;

  pythia.readString("Beams:allowIDASwitch = on");
  double idA = idA_list[0]; // pythia_.settings.mode("Beams:idA", 2212)
  double idB = idB_list[0]; // pythia_.settings.mode("Beams:idB", 1000070140);
  pythia.settings.mode("Beams:idA", idA);
  pythia.settings.mode("Beams:idB", idB);
  
  double pA_max = 1e5;
  double mA_max = pythia.particleData.m0(3322);
  double eA_max = sqrt(pow2(pA_max) + pow2(mA_max));
  double mB = pythia.particleData.m0(2212);
  double eCM_max = sqrt(pow2(mA_max) + pow2(mB) + 2.*eA_max*mB);

  pythia.readString("Beams:allowVariableEnergy = on");
  pythia.readString("Beams:frameType = 1");
  pythia.settings.parm("Beams:eCM", 10000.);

  pythia.readString("SoftQCD:all = on");

  pythia.readString("MultipartonInteractions:reuseInit = 3");
  pythia.readString("MultipartonInteractions:initFile = main424.mpi");
  pythia.readString("HeavyIon:SasdMpiReuseInit = 3");
  pythia.readString("HeavyIon:SasdMpiInitFile = main424.sasd.mpi");
  pythia.readString("HeavyIon:SigFitReuseInit = 3");
  pythia.readString("HeavyIon:SigFitInitFile = main424.sigfit");

  // replicate exact Corsika 8 Interaction.inl which were not used in my Pythia tests before
  pythia.readString("HadronLevel:Decay = off");
  pythia.readString("Check:epTolErr = 0.1");
  pythia.readString("Check:epTolWarn = 0.0001");
  pythia.readString("Check:mTolErr = 0.01");
  pythia.readString("Stat:showProcessLevel = off");
  pythia.readString("Stat:showPartonLevel = off");

  if (!pythia.init()) {
    throw std::runtime_error("Pythia::Interaction: Initialization failed!");
    return -1;
  }

  //cout << "starting event loop..." << endl; 
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    std::cout << "Event: " << iEvent << std::endl;  

    /*double idA_switch = idA_list[iEvent%4];
    double idB_switch = idB_list[iEvent%3];
    std::cout << "switching beam " << std::endl;
    pythia.setBeamIDs(idA_switch, idB_switch);*/

    /*double pA_switch = pA_min + pythia.rndm.flat() * (pA_max - pA_min);
    double mA_switch = pythia.particleData.m0(idA_switch);
    double eA_switch = sqrt(pow2(pA_switch) + pow2(mA_switch));
    double eCM_switch = sqrt(2.* eA_switch * mB);
    std::cout << "setting kinematics " << std::endl;
    pythia.setKinematics(eCM_switch);*/
    
    //pythia.next();
    pythia.event.list();
  }

  pythia.stat();

  if (doLog) std::cout.rdbuf(oldCout);

  return 0;
}
