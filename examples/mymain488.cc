
// Authors: Chloé Gaudu <gaudu@uni-wuppertal.de>

// Keywords: switch beam; switch collision energy; reuse MPI initialization

// Attempt to replicate the Corsika 8 Interaction.inl file version from the
// 620-upgrade-to-pythia-8-310 branch on 29 July 2024 to test out the reuse 
// files. Models the idA and idB used in testPythia8Interaction.inl. 

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Code const projectile = GENERATE(Code::Proton, Code::AntiProton, Code::Neutron, Code::PiPlus, Code::KPlus);
  vector<int> idA_list = {2212, -2212, 2112, 211, 321};
  // Code const target = GENERATE(Code::Proton, Code::Nitrogen, Code::Oxygen, Code::Argon); 
  vector<int> idB_list = {2212, 1000070140, 1000080160, 1000180400}; 

  const bool doLog = true;
  std::string out = "mymain488_c8_Kp_p_proj400TeV";
  ofstream logBuf;
  std::streambuf* oldCout;
  if(doLog) {
    oldCout = std::cout.rdbuf(logBuf.rdbuf());
    logBuf.open((out == "" ? "pythia.log" : out + ".log"));
  }

  Pythia pythia;

  pythia.readString("Beams:allowIDASwitch = on");
  double idA = idA_list[0];
  double idB = idB_list[1];
  pythia.settings.mode("Beams:idA", idA);
  pythia.settings.mode("Beams:idB", idB);

  pythia.readString("Beams:allowVariableEnergy = on");
  pythia.readString("Beams:frameType = 1");
  pythia.settings.parm("Beams:eCM", 10000.);

  pythia.readString("SoftQCD:all = on");

  pythia.readString("MultipartonInteractions:reuseInit = 2");
  pythia.readString("MultipartonInteractions:initFile = main424.mpi");
  pythia.readString("HeavyIon:SasdMpiReuseInit = 2");
  pythia.readString("HeavyIon:SasdMpiInitFile = main424.sasd.mpi");
  pythia.readString("HeavyIon:SigFitReuseInit = 2");
  pythia.readString("HeavyIon:SigFitInitFile = main424.sigfit");

  // replicate exact Corsika 8 Interaction.inl which were not used in my Pythia tests before
  pythia.readString("HadronLevel:Decay = off");
  pythia.readString("Check:epTolErr = 0.1");
  pythia.readString("Check:epTolWarn = 0.0001");
  pythia.readString("Check:mTolErr = 0.01");
  pythia.readString("Stat:showProcessLevel = off");
  pythia.readString("Stat:showPartonLevel = off");

  if (!pythia.init()) {
    if (doLog) std::cout.rdbuf(oldCout);  
    throw std::runtime_error("Pythia::Interaction: Initialization failed!");
    return -1;
  }

  /*if (!pythia.next()) {
      if (doLog) std::cout.rdbuf(oldCout);
      throw std::runtime_error("Pythia::Interaction: Interaction failed!");
  }*/
  
  double idA_pythia = 321;
  double idB_pythia = 2212;
  double eCM_pythia = sqrt(pow2(pythia.particleData.m0(idA_pythia)) + pow2(pythia.particleData.m0(idB_pythia))+ 2.*400000.*pythia.particleData.m0(idB_pythia));
  pythia.setBeamIDs(idA_pythia, idB_pythia);
  std::cout << "eCM_pythia = " << eCM_pythia << std::endl; 
  pythia.setKinematics(eCM_pythia);

  if (!pythia.next()) {
      if (doLog) std::cout.rdbuf(oldCout);
      throw std::runtime_error("Pythia::Interaction: Interaction failed!");
  }

  pythia.event.list();
  pythia.stat();

  if (doLog) std::cout.rdbuf(oldCout);

  return 0;
}
