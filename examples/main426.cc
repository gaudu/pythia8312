// main426.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Marius Utheim <marius.m.utheim@jyu.fi>;
//          Chloé Gaudu <gaudu@uni-wuppertal.de>

// Keywords: hadron-ion collisions, optimization

// The purpose of this example is to generate initialization files that can
// be used to speed up initialization in hadron-hadron or hadron-ion runs.
// It produces data for energies from 10^2 to 10^12 GeV. 
// All hadron-nucleon and hadron-ion interactions are possible. 
// Initialization data is saved in all.mpi, all.sasd.mpi, and all.sigfit.

// After initializing, it is possible to change energy and beam types on
// a per-event basis.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  double pA_max = 1e12;
  double pA_min = 1e2; 
  int nPts = 11;

  vector<int> idA_list = {111, 211, -211, 311, 321, -321, 130, 310,
                          2212, -2212, 2112, -2112,
                          3122, 3212, 3222, 3112, 3322, 
                          1000020040, 1000070140, 1000260560
  };

  const bool doLog = true;
  std::string out = "main426";
  ofstream logBuf;
  std::streambuf* oldCout;
  if(doLog) {
    oldCout = std::cout.rdbuf(logBuf.rdbuf());
    logBuf.open((out == "" ? "pythia.log" : out + ".log"));
  }

  Pythia pythia;
  // Use Angantyr even when initializing with pp.
  pythia.readString("HeavyIon:mode = 2");

  double mA_max = pythia.particleData.m0(3322); // Xi0 heaviest
  double mA_min = pythia.particleData.m0(111); // pi0 lightest
  double mB = pythia.particleData.m0(2212);
  double eA_max = sqrt(pow2(pA_max) + pow2(mA_max)); 
  double eA_min = sqrt(pow2(pA_min) + pow2(mA_min)); 
  double eCM_max = sqrt(pow2(mA_max) + pow2(mB) + 2.* eA_max * mB);
  double eCM_min = sqrt(pow2(mA_min) + pow2(mB) + 2.* eA_min * mB);

  // Variable energy parameters.
  pythia.readString("Beams:allowVariableEnergy = on");
  pythia.settings.parm("HeavyIon:varECMMax", eCM_max);
  pythia.settings.parm("HeavyIon:varECMMin", eCM_min);
  pythia.settings.parm("HeavyIon:varECMSigFitNPts", nPts);

  // Variable beam parameters.
  pythia.readString("Beams:allowIDASwitch = on");
  pythia.settings.mvec("Beams:idAList", idA_list);

  // Specify where to save. If you set reuseInit = 2, the old files will be
  // replaced overwritten they already exist.
  pythia.readString("MultipartonInteractions:reuseInit = 3");
  pythia.readString("MultipartonInteractions:initFile = main426.mpi");
  pythia.readString("HeavyIon:SasdMpiReuseInit = 3");
  pythia.readString("HeavyIon:SasdMpiInitFile = main426.sasd.mpi");
  pythia.readString("HeavyIon:SigFitReuseInit = 3");
  pythia.readString("HeavyIon:SigFitInitFile = main426.sigfit");

  if (!pythia.init()) {
    cout << " Pythia failed to initialize." << endl;
    return -1;
  }
  
  if (doLog) std::cout.rdbuf(oldCout);

  // After initializing, we're done.
  return 0;
}
