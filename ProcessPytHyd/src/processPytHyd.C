//Author: Chris McGinn (2019.10.07)

//cpp
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

//FASTJET
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "fastjet/contrib/SoftKiller.hh"

#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

//ROOT
#include "TDatime.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TNamed.h"
#include "TTree.h"

//JetMethods
#include "MakePytHyd/include/checkMakeDir.h"
#include "MakePytHyd/include/returnRootFileContentsList.h"
#include "ProcessPytHyd/include/etaPhiFunc.h"

int processPytHyd(const std::string inFileName)
{
  if(!checkFileExt(inFileName, "root")) return 1;

  std::string outFileName = inFileName.substr(0, inFileName.find(".root"));
  while(outFileName.find("/") != std::string::npos){
    outFileName.replace(0, outFileName.find("/")+1, "");
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  outFileName = "output/" + dateStr + "/" + outFileName + "_PROCESSED_" + dateStr + ".root";

  const double minJtPt = 30.;
  const double maxJtAbsEta = 2.;
  const double maxGlobalAbsEta = 5.;
  const double rParam = 0.4;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, rParam, fastjet::E_scheme);
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* processedTree_p = new TTree("processedTree_p", "");

  Float_t sumHFEt_;
  Int_t cent_;

  const Int_t nMaxJets = 500;
  Int_t njtSig_;
  Float_t jtptSig_[nMaxJets];
  Float_t jtetaSig_[nMaxJets];
  Float_t jtphiSig_[nMaxJets];
  Bool_t hasSigBkgdMatchSig_[nMaxJets];
  Bool_t hasSigBkgdAMatchSig_[nMaxJets];
  Bool_t hasSigBkgdSKMatchSig_[nMaxJets];

  const Int_t nMaxParticles = 100000;
  Int_t nconstSig_;
  Float_t constptSig_[nMaxParticles];
  Float_t constetaSig_[nMaxParticles];
  Float_t constphiSig_[nMaxParticles];
  Int_t constjtposSig_[nMaxParticles];
  
  Int_t njtSigBkgd_;
  Float_t jtptSigBkgd_[nMaxJets];
  Float_t jtetaSigBkgd_[nMaxJets];
  Float_t jtphiSigBkgd_[nMaxJets];
  Int_t sigposSigBkgd_[nMaxJets];
  
  Int_t njtSigBkgdA_;
  Float_t jtptSigBkgdA_[nMaxJets];
  Float_t jtetaSigBkgdA_[nMaxJets];
  Float_t jtphiSigBkgdA_[nMaxJets];
  Int_t sigposSigBkgdA_[nMaxJets];

  Int_t njtSigBkgdSK_;
  Float_t jtptSigBkgdSK_[nMaxJets];
  Float_t jtetaSigBkgdSK_[nMaxJets];
  Float_t jtphiSigBkgdSK_[nMaxJets];
  Int_t sigposSigBkgdSK_[nMaxJets];

  processedTree_p->Branch("sumHFEt", &sumHFEt_, "sumHFEt/F");
  processedTree_p->Branch("cent", &cent_, "cent/I");

  processedTree_p->Branch("njtSig", &njtSig_, "njtSig/I");
  processedTree_p->Branch("jtptSig", jtptSig_, "jtptSig[njtSig]/F");
  processedTree_p->Branch("jtetaSig", jtetaSig_, "jtetaSig[njtSig]/F");
  processedTree_p->Branch("jtphiSig", jtphiSig_, "jtphiSig[njtSig]/F");
  processedTree_p->Branch("hasSigBkgdMatchSig", hasSigBkgdMatchSig_, "hasSigBkgdMatchSig[njtSig]/O");
  processedTree_p->Branch("hasSigBkgdAMatchSig", hasSigBkgdAMatchSig_, "hasSigBkgdAMatchSig[njtSig]/O");
  processedTree_p->Branch("hasSigBkgdSKMatchSig", hasSigBkgdSKMatchSig_, "hasSigBkgdSKMatchSig[njtSig]/O");

  processedTree_p->Branch("nconstSig", &nconstSig_, "nconstSig/I");
  processedTree_p->Branch("constptSig", constptSig_, "constptSig[nconstSig]/F");
  processedTree_p->Branch("constphiSig", constphiSig_, "constphiSig[nconstSig]/F");
  processedTree_p->Branch("constetaSig", constetaSig_, "constetaSig[nconstSig]/F");
  processedTree_p->Branch("constjtposSig", constjtposSig_, "constjtposSig[nconstSig]/I");
  
  processedTree_p->Branch("njtSigBkgd", &njtSigBkgd_, "njtSigBkgd/I");
  processedTree_p->Branch("jtptSigBkgd", jtptSigBkgd_, "jtptSigBkgd[njtSigBkgd]/F");
  processedTree_p->Branch("jtetaSigBkgd", jtetaSigBkgd_, "jtetaSigBkgd[njtSigBkgd]/F");
  processedTree_p->Branch("jtphiSigBkgd", jtphiSigBkgd_, "jtphiSigBkgd[njtSigBkgd]/F");
  processedTree_p->Branch("sigposSigBkgd", sigposSigBkgd_, "sigposSigBkgd[njtSigBkgd]/I");

  processedTree_p->Branch("njtSigBkgdA", &njtSigBkgdA_, "njtSigBkgdA/I");
  processedTree_p->Branch("jtptSigBkgdA", jtptSigBkgdA_, "jtptSigBkgdA[njtSigBkgdA]/F");
  processedTree_p->Branch("jtetaSigBkgdA", jtetaSigBkgdA_, "jtetaSigBkgdA[njtSigBkgdA]/F");
  processedTree_p->Branch("jtphiSigBkgdA", jtphiSigBkgdA_, "jtphiSigBkgdA[njtSigBkgdA]/F");
  processedTree_p->Branch("sigposSigBkgdA", sigposSigBkgdA_, "sigposSigBkgdA[njtSigBkgdA]/I");

  processedTree_p->Branch("njtSigBkgdSK", &njtSigBkgdSK_, "njtSigBkgdSK/I");
  processedTree_p->Branch("jtptSigBkgdSK", jtptSigBkgdSK_, "jtptSigBkgdSK[njtSigBkgdSK]/F");
  processedTree_p->Branch("jtetaSigBkgdSK", jtetaSigBkgdSK_, "jtetaSigBkgdSK[njtSigBkgdSK]/F");
  processedTree_p->Branch("jtphiSigBkgdSK", jtphiSigBkgdSK_, "jtphiSigBkgdSK[njtSigBkgdSK]/F");
  processedTree_p->Branch("sigposSigBkgdSK", sigposSigBkgdSK_, "sigposSigBkgdSK[njtSigBkgdSK]/I");

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> inNameList = returnRootFileContentsList(inFile_p, "TNamed");
  std::vector<TNamed> outNameList;
  for(unsigned int nI = 0; nI < inNameList.size(); ++nI){
    std::string name = inNameList[nI].substr(inNameList[nI].find("/")+1, inNameList[nI].size());
    std::string title = ((TNamed*)inFile_p->Get(inNameList[nI].c_str()))->GetTitle();
    std::cout << name << ", " << title << std::endl;
    outNameList.push_back(TNamed(name.c_str(), title.c_str()));
  }
  
  std::vector<std::string> treeList = returnRootFileContentsList(inFile_p, "TTree");
  if(treeList.size() != 1){
    std::cout << "Given input \'" << inFileName << "\' contains more than or less than 1 ttree. return 1" << std::endl;
    return 1;
  }
  TTree* inTree_p = (TTree*)inFile_p->Get(treeList[0].c_str());

  Int_t nPart_;
  Float_t pt_[nMaxParticles];
  Float_t eta_[nMaxParticles];
  Float_t phi_[nMaxParticles];
  Float_t mass_[nMaxParticles];
  Bool_t isSig_[nMaxParticles];

  inTree_p->SetBranchStatus("*", 0);
  inTree_p->SetBranchStatus("nPart", 1);
  inTree_p->SetBranchStatus("pt", 1);
  inTree_p->SetBranchStatus("eta", 1);
  //  inTree_p->SetBranchStatus("phi", 1);
  //  inTree_p->SetBranchStatus("mass", 1);
  //  inTree_p->SetBranchStatus("isSig", 1);

  inTree_p->SetBranchAddress("nPart", &nPart_);
  inTree_p->SetBranchAddress("pt", pt_);
  inTree_p->SetBranchAddress("eta", eta_);
  //  inTree_p->SetBranchAddress("phi", phi_);
  //  inTree_p->SetBranchAddress("mass", mass_);
  //  inTree_p->SetBranchAddress("isSig", isSig_);

  const Int_t nEntries = TMath::Min(1000000, (Int_t)inTree_p->GetEntries());
  const Int_t nDiv = TMath::Max(1, nEntries/20);

  std::vector<float> sumHFEtVect;
  
  std::cout << "Pre-Processing " << nEntries << " events..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;
    inTree_p->GetEntry(entry);

    sumHFEt_ = 0.0;

    for(Int_t pI = 0; pI < nPart_; ++pI){
      if(TMath::Abs(eta_[pI]) > 3.) sumHFEt_ += pt_[pI];
    }

    sumHFEtVect.push_back(sumHFEt_);
  }

  std::sort(std::begin(sumHFEtVect), std::end(sumHFEtVect));
  const Int_t nCent = 200;
  double interval = sumHFEtVect.size()/nCent;
  std::vector<float> centVect = {0.0};
  for(Int_t cI = 1; cI < nCent; ++cI){
    centVect.push_back(sumHFEtVect[(int)(cI*interval)]);
  }
  centVect.push_back(sumHFEtVect[sumHFEtVect.size()-1] + 10.);

  std::cout << "Centrality: ";
  for(unsigned int cI = 0; cI < centVect.size()-1; ++cI){
    std::cout << centVect[cI] << ", ";
  }
  std::cout << centVect[centVect.size()-1] << "." << std::endl;
    
    
  inTree_p->SetBranchStatus("phi", 1);
  inTree_p->SetBranchStatus("mass", 1);
  inTree_p->SetBranchStatus("isSig", 1);

  inTree_p->SetBranchAddress("phi", phi_);
  inTree_p->SetBranchAddress("mass", mass_);
  inTree_p->SetBranchAddress("isSig", isSig_);
  
  std::cout << "Processing " << nEntries << " events..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;
    inTree_p->GetEntry(entry);

    sumHFEt_ = 0.0;

    std::vector<fastjet::PseudoJet> sigParticles;
    std::vector<fastjet::PseudoJet> particles;
    TLorentzVector tL;
    
    for(Int_t pI = 0; pI < nPart_; ++pI){
      if(TMath::Abs(eta_[pI]) > 3.) sumHFEt_ += pt_[pI];

      tL.SetPtEtaPhiM(pt_[pI], eta_[pI], phi_[pI], mass_[pI]);
      fastjet::PseudoJet part(tL.Px(), tL.Py(), tL.Pz(), tL.E());
      
      if(isSig_[pI]) sigParticles.push_back(part);
      particles.push_back(part);
    }

    cent_ = -1;
    for(unsigned int cI = 1; cI < centVect.size(); ++cI){
      if(sumHFEt_ < centVect[cI]){
	cent_ = 200 - cI;
	break;
      }
    }

    //area based stuff
    fastjet::GridMedianBackgroundEstimator bge(maxGlobalAbsEta, 0.5); //see FASTJET manual
    bge.set_compute_rho_m(true);
    bge.set_particles(particles);
    
    fastjet::Subtractor subtractor(&bge);
    subtractor.set_use_rho_m(true);
    subtractor.set_safe_mass(true);

    fastjet::contrib::SoftKiller soft_killer(5.0, 0.4);
    std::vector<fastjet::PseudoJet> soft_killed_event = soft_killer(particles);

    
    fastjet::ClusterSequence csSig(sigParticles, jet_def);
    std::vector<fastjet::PseudoJet> jetsSig = fastjet::sorted_by_pt(csSig.inclusive_jets());


    njtSig_ = 0;
    nconstSig_ = 0;
    for(unsigned int jI = 0; jI < jetsSig.size(); ++jI){
      if(jetsSig[jI].pt() < minJtPt) continue;
      if(TMath::Abs(jetsSig[jI].eta()) >= maxJtAbsEta) continue;

      jtptSig_[njtSig_] = jetsSig[jI].pt();
      jtetaSig_[njtSig_] = jetsSig[jI].eta();
      jtphiSig_[njtSig_] = jetsSig[jI].phi_std();
      hasSigBkgdMatchSig_[njtSig_] = false;
      hasSigBkgdAMatchSig_[njtSig_] = false;
      hasSigBkgdSKMatchSig_[njtSig_] = false;
      
      std::vector<fastjet::PseudoJet> constit = jetsSig[jI].constituents();
      for(unsigned int cI = 0; cI < constit.size(); ++cI){
	constptSig_[nconstSig_] = constit[cI].pt();
	constetaSig_[nconstSig_] = constit[cI].eta();
	constphiSig_[nconstSig_] = constit[cI].phi_std();
	constjtposSig_[nconstSig_] = njtSig_;
	
	++nconstSig_;
	if(nconstSig_ >= nMaxParticles){
	  std::cout << "WARNING YOU FUCKED UP YOU FUCKED UP YOU FUCKED UP" << std::endl;
	}
      }

      ++njtSig_;
    }

    fastjet::ClusterSequence csSigBkgd(particles, jet_def);
    std::vector<fastjet::PseudoJet> jetsSigBkgd = fastjet::sorted_by_pt(csSigBkgd.inclusive_jets());

    njtSigBkgd_ = 0;
    for(unsigned int jI = 0; jI < jetsSigBkgd.size(); ++jI){
      if(jetsSigBkgd[jI].pt() < minJtPt) continue;
      if(TMath::Abs(jetsSigBkgd[jI].eta()) >= maxJtAbsEta) continue;

      jtptSigBkgd_[njtSigBkgd_] = jetsSigBkgd[jI].pt();
      jtetaSigBkgd_[njtSigBkgd_] = jetsSigBkgd[jI].eta();
      jtphiSigBkgd_[njtSigBkgd_] = jetsSigBkgd[jI].phi_std();
      sigposSigBkgd_[njtSigBkgd_] = -1;

      for(Int_t jI2 = 0; jI2 < njtSig_; ++jI2){
	if(hasSigBkgdMatchSig_[jI2]) continue;

	if(getDR(jtetaSig_[jI2], jtphiSig_[jI2], jtetaSigBkgd_[jI], jtphiSigBkgd_[jI]) < rParam/2.){

	  sigposSigBkgd_[njtSigBkgd_] = jI2;
	  hasSigBkgdMatchSig_[jI2] = true;
	  break;
	}     
      }
      
      ++njtSigBkgd_;
    }

    double ghost_area = 0.01;
    int active_area_repeats = 1;
    fastjet::GhostedAreaSpec ghost_spec(maxGlobalAbsEta, active_area_repeats, ghost_area);
    fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, ghost_spec);

    fastjet::ClusterSequenceArea csSigBkgdA(particles, jet_def, area_def);
    std::vector<fastjet::PseudoJet> jetsSigBkgdA = fastjet::sorted_by_pt(subtractor(csSigBkgdA.inclusive_jets()));

    njtSigBkgdA_ = 0;
    for(unsigned int jI = 0; jI < jetsSigBkgdA.size(); ++jI){
      if(jetsSigBkgdA[jI].pt() < minJtPt) continue;
      if(TMath::Abs(jetsSigBkgdA[jI].eta()) >= maxJtAbsEta) continue;

      jtptSigBkgdA_[njtSigBkgdA_] = jetsSigBkgdA[jI].pt();
      jtetaSigBkgdA_[njtSigBkgdA_] = jetsSigBkgdA[jI].eta();
      jtphiSigBkgdA_[njtSigBkgdA_] = jetsSigBkgdA[jI].phi_std();
      sigposSigBkgdA_[njtSigBkgdA_] = -1;

      for(Int_t jI2 = 0; jI2 < njtSig_; ++jI2){
	if(hasSigBkgdAMatchSig_[jI2]) continue;

	if(getDR(jtetaSig_[jI2], jtphiSig_[jI2], jtetaSigBkgdA_[jI], jtphiSigBkgdA_[jI]) < rParam/2.){

	  sigposSigBkgdA_[njtSigBkgdA_] = jI2;
	  hasSigBkgdAMatchSig_[jI2] = true;
	  break;
	}     
      }
      
      ++njtSigBkgdA_;
    }


    fastjet::ClusterSequenceArea csSigBkgdSK(soft_killed_event, jet_def, area_def);
    std::vector<fastjet::PseudoJet> jetsSigBkgdSK = fastjet::sorted_by_pt(csSigBkgdSK.inclusive_jets());

    njtSigBkgdSK_ = 0;
    for(unsigned int jI = 0; jI < jetsSigBkgdSK.size(); ++jI){
      if(jetsSigBkgdSK[jI].pt() < minJtPt) continue;
      if(TMath::Abs(jetsSigBkgdSK[jI].eta()) >= maxJtAbsEta) continue;

      jtptSigBkgdSK_[njtSigBkgdSK_] = jetsSigBkgdSK[jI].pt();
      jtetaSigBkgdSK_[njtSigBkgdSK_] = jetsSigBkgdSK[jI].eta();
      jtphiSigBkgdSK_[njtSigBkgdSK_] = jetsSigBkgdSK[jI].phi_std();
      sigposSigBkgdSK_[njtSigBkgdSK_] = -1;

      for(Int_t jI2 = 0; jI2 < njtSig_; ++jI2){
	if(hasSigBkgdSKMatchSig_[jI2]) continue;

	if(getDR(jtetaSig_[jI2], jtphiSig_[jI2], jtetaSigBkgdSK_[jI], jtphiSigBkgdSK_[jI]) < rParam/2.){

	  sigposSigBkgdSK_[njtSigBkgdSK_] = jI2;
	  hasSigBkgdSKMatchSig_[jI2] = true;
	  break;
	}     
      }
      
      ++njtSigBkgdSK_;
    }

    processedTree_p->Fill();
  }  

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  processedTree_p->Write("", TObject::kOverwrite);
  delete processedTree_p;

  TDirectoryFile* params_p = (TDirectoryFile*)outFile_p->mkdir("params");
  params_p->cd();
  
  for(unsigned int oI = 0; oI < outNameList.size(); ++oI){
    outNameList[oI].Write("", TObject::kOverwrite);
  }

  params_p->Close();
  delete params_p;
    
  outFile_p->Close();
  delete outFile_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/processPytHyd.exe <inFileName>. return 1" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += processPytHyd(argv[1]);
  return retVal;
}
