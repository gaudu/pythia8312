// Authors: Chloé Gaudu <gaudu@uni-wuppertal.de>

// Keyword: analysis; hA collisions; AA collisions; cross-section

// Usage: ./main1011 Beams:idA Beams:idB p_lab(projectile)
// can be run using ./run_pleiades_xsec_script_8312 on Pleiades

#include "Pythia8/HeavyIons.h"
#include "Pythia8/Pythia.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

template <typename K, typename V>
std::vector<K> get_map_keys(const std::map<K, V>& input_map) {
    std::vector<K> keys;
    for (const auto& pair : input_map) {
        keys.push_back(pair.first);
    }
    return keys;
}

std::string to_scientific_string(double input) { 
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(0) << input;
    std::string result = stream.str();
    size_t dotPos = result.find('.');
    if (dotPos != std::string::npos) {
        result.erase(dotPos, 1); 
    }
    size_t pos = result.find("e+");
    if (pos != std::string::npos) {
        result.erase(pos + 1, 1); 
    }

    return result;
}

using namespace Pythia8;

int main(int argc, char *argv[]) {
#include <stdio.h>
#include <stdlib.h>
  if (argc<=1) {
    std::cout << "argc<=1" << '\n';
    if (argc==5) {
      std::cout << "argc==5" << '\n'; 
    }
    return 1; 
  }

  float nEvents = 1e5;
  int idA = atoi(argv[1]);
  int idB = atoi(argv[2]);

  std::map<std::string, std::string> idA_map = {
    {"2212", "p"}, 
    {"-2212", "ap"}, 
    {"2112", "n"}, 
    {"-2112", "an"},
    {"111", "pi0"}, 
    {"211", "pip"}, 
    {"-211", "pim"},
    {"311", "k0"},
    {"321", "kp"}, 
    {"-321", "km"}, 
    {"130", "k0l"}, 
    {"310", "k0s"},
    {"3122", "lambda0"}, 
    {"3212", "sigma0"},
    {"3222", "sigmap"}, 
    {"3112", "sigmam"}, 
    {"3322", "xi0"},
    {"1000020040", "4He"}, 
    {"1000070140", "14N"}, 
    {"1000260560", "56Fe"}
    };
  std::map<std::string, std::string> idB_map = {
    {"2212", "1H"}, 
    {"1000060120", "12C"}, 
    {"1000070140", "14N"}, 
    {"1000080160", "16O"}, 
    {"1000180400", "40Ar"}
    };
  std::map<int, std::string> subprocess_map = {
    {101, "nd"}, 
    {102, "el"}, 
    {103, "XB"}, 
    {104, "AX"}, 
    {105, "XX"}, 
    {106, "AXB"}
    };

  // logfile initialization
  const bool doLog = true;
  std::string out = "main1011_sNN_v2_"+idA_map[std::string(argv[1])]+"_"+idB_map[std::string(argv[2])]+"_pNow"+to_scientific_string(std::atof(argv[3]))+"_nEvents"+to_scientific_string(nEvents);
  ofstream logBuf;
  std::streambuf* oldCout;
  if(doLog) {
    oldCout = std::cout.rdbuf(logBuf.rdbuf());
    logBuf.open((out == "" ? "pythia.log" : out + ".log"));
  }
  
  // Pythia configuration
  Pythia pythia; 

  // "effective" nucleon mass m_N = m_A / A
  std::map<std::string, double> mB_map = {
    {"2212", pythia.particleData.m0(2212)}, 
    {"1000060120", pythia.particleData.m0(1000060120)/12}, 
    {"1000070140", pythia.particleData.m0(1000070140)/14}, 
    {"1000080160", pythia.particleData.m0(1000080160)/16}, 
    {"1000180400", pythia.particleData.m0(1000180400)/40}
    };

  pythia.settings.mode("Beams:idA", idA);
  pythia.settings.mode("Beams:idB", idB);
  
  double pNow = std::stod(argv[3]); // GeV
  double eCMNow = ( Vec4(0., 0., pNow, pNow * sqrt(1 + pow2(mB_map[std::string(argv[2])] / pNow)))
                 + Vec4(0., 0., 0., mB_map[std::string(argv[2])]) ).mCalc(); // GeV
  pythia.settings.parm("Beams:eCM", eCMNow);
  pythia.readString("SoftQCD:all = on");

  // use Angantyr for minimum-bias pp collisions
  //pythia.readString("HeavyIon:mode = 2");

  if (!pythia.init()) {
    std::cout << "Pythia failed to initialize." << std::endl;
    return -1;
  }

  std::ofstream xz{out + ".dat"};

  // loop over events
  for ( int iEvent = 0; iEvent < nEvents; ++iEvent ) { 
    /*pythia.event.list();
    pythia.info.list();
    pythia.process.list();*/
    if (!pythia.next()) continue;
  }
  
  std::vector<int> sigmaGen_keys = get_map_keys(pythia.info.sigGenM);

  pythia.stat();

  if (!xz.is_open()) {
    std::cerr << "Error: failed to open file for writing." << std::endl;
    return 1; 
  }

  xz << "plab" << '\t' << argv[3] << '\n';
  xz << "s_NN" << '\t' << eCMNow << '\n';
  xz << "proj_id" << '\t' << argv[1] << '\n';
  xz << "targ_id" << '\t' << argv[2] << '\n';
  xz << "sig_tot" << '\t' << pythia.info.sigmaGen(0) << '\n';
  xz << "err_tot" << '\t' << pythia.info.sigmaErr(0) << '\n';
  for (const auto& key : sigmaGen_keys) {
        std::string label_sig = "sig_" + subprocess_map[key];
        std::string label_err = "err_" + subprocess_map[key];
        xz << label_sig << '\t' << pythia.info.sigmaGen(key) << '\n';
        xz << label_err << '\t' << pythia.info.sigmaErr(key) << '\n';
    }
  bool has102 = std::find(sigmaGen_keys.begin(), sigmaGen_keys.end(), 102) != sigmaGen_keys.end();
  if (has102) {
  xz << "sig_inel" << '\t' << pythia.info.sigmaGen(0) - pythia.info.sigmaGen(102) << '\n';
  xz << "err_inel" << '\t' << sqrt(pow(pythia.info.sigmaErr(0),2) + pow(pythia.info.sigmaErr(102),2)) << '\n';
  }

  if (!xz) {
    std::cerr << "Error: failed to write data to the file." << std::endl;
    return 1; 
  }
  
  xz.close();
  std::cout << "Data written to output file successfully." << std::endl;
  
  // logfile output
  if (doLog) std::cout.rdbuf(oldCout);

  return 0;

}
