//Author: Chris McGinn (2019.10.04)

//cpp
#include <iostream>
#include <string>

//ROOT
#include "TDatime.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TNamed.h"
#include "TTree.h"

//Local
#include "MakePytHyd/include/checkMakeDir.h"
#include "MakePytHyd/include/configParser.h"
#include "MakePytHyd/include/generalTreeHandler.h"
#include "MakePytHyd/include/returnRootFileContentsList.h"

int combinePythiaWithBkgd(const std::string configName)
{
  if(!checkFile(configName)){
    std::cout << "Given configName \'" << configName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  configParser config;
  if(!config.Init(configName)) return 1;

  bool doBkgd = true;
  
  const std::string inPytFileName = config.GetParamValByName("INPYTFILENAME");
  const std::string inBkgdFileName = config.GetParamValByName("INBKGDFILENAME");  
  if(!checkFileExt(inPytFileName, "root")) return 1;
  if(!checkFileExt(inBkgdFileName, "root")) doBkgd = false;  // we need pure pythia mode

  const std::string inPytTreeName = config.GetParamValByName("INPYTTREENAME");
  const std::string inBkgdTreeName = config.GetParamValByName("INBKGDTREENAME");  

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  std::string outFileName = "output/" + dateStr + "/combinedPythiaWithBkgd_" + dateStr + ".root";
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* outTree_p = new TTree("particleTree", "");

  const Int_t nMaxQGParticles = 100000;
  Int_t nQGPart_ = 0;
  Float_t qgpt_[nMaxQGParticles];
  Float_t qgeta_[nMaxQGParticles];
  Float_t qgphi_[nMaxQGParticles];
  Int_t qgpdgid_[nMaxQGParticles];

  const Int_t nMaxParticles = 100000;
  Int_t nPart_ = 0;
  Float_t pt_[nMaxParticles];
  Float_t eta_[nMaxParticles];
  Float_t phi_[nMaxParticles];
  Float_t mass_[nMaxParticles];
  Int_t pdgid_[nMaxParticles];
  Short_t chg_[nMaxParticles];
  Bool_t isSig_[nMaxParticles];

  outTree_p->Branch("nQGPart", &nQGPart_, "nQGPart/I");
  outTree_p->Branch("qgpt", qgpt_, "qgpt[nPart]/F");
  outTree_p->Branch("qgphi", qgphi_, "qgphi[nPart]/F");
  outTree_p->Branch("qgeta", qgeta_, "qgeta[nPart]/F");
  outTree_p->Branch("qgpdgid", qgpdgid_, "qgpdgid[nPart]/I");

  outTree_p->Branch("nPart", &nPart_, "nPart/I");
  outTree_p->Branch("pt", pt_, "pt[nPart]/F");
  outTree_p->Branch("phi", phi_, "phi[nPart]/F");
  outTree_p->Branch("eta", eta_, "eta[nPart]/F");
  outTree_p->Branch("mass", mass_, "mass[nPart]/F");
  outTree_p->Branch("pdgid", pdgid_, "pdgid[nPart]/I");
  outTree_p->Branch("chg", chg_, "chg[nPart]/S");
  outTree_p->Branch("isSig", isSig_, "isSig[nPart]/O");

  TFile* inPytFile_p = new TFile(inPytFileName.c_str(), "READ");
  std::vector<std::string> pytTreeList = returnRootFileContentsList(inPytFile_p, "TTree");
  if(!vectContainsStr(inPytTreeName, &pytTreeList)){
    std::cout << "Requested treeName \'" << inPytTreeName << "\' is not found in \'" << inPytFileName << "\'. return 1" << std::endl;
    return 1;
  }
  TTree* pytTree_p = (TTree*)inPytFile_p->Get(inPytTreeName.c_str());
  generalTreeHandler pytHandler;
  if(!pytHandler.Init(config, true, pytTree_p)) return 1;

  TFile* inBkgdFile_p = nullptr;
  TTree* bkgdTree_p = nullptr;
  generalTreeHandler bkgdHandler;
  
  if(doBkgd){
    inBkgdFile_p = new TFile(inBkgdFileName.c_str(), "READ");
    std::vector<std::string> bkgdTreeList = returnRootFileContentsList(inBkgdFile_p, "TTree");
    if(!vectContainsStr(inBkgdTreeName, &bkgdTreeList)){
      std::cout << "Requested treeName \'" << inBkgdTreeName << "\' is not found in \'" << inBkgdFileName << "\'. return 1" << std::endl;
      return 1;
    }
    bkgdTree_p = (TTree*)inBkgdFile_p->Get(inBkgdTreeName.c_str());
    if(!bkgdHandler.Init(config, false, bkgdTree_p)) return 1;
  }
  
  const Int_t nEntriesPyt = pytTree_p->GetEntries();
  Int_t nEntriesBkgd = pytTree_p->GetEntries();
  if(doBkgd) nEntriesBkgd = bkgdTree_p->GetEntries();
  const Int_t nEntries = TMath::Min(10000000, TMath::Max(nEntriesPyt, nEntriesBkgd));
  const Int_t nEntriesMin = TMath::Min(nEntriesPyt, nEntriesBkgd);
  const Int_t nDiv = TMath::Max(1, nEntries/200);

  std::cout << "Processing " << nEntries << " events, signal+bkgd" << std::endl;
  std::cout << " Mismatch between signal and background will mean we reuse events " << ((double)nEntries)/((double)nEntriesMin) << " times" << std::endl;

  std::vector<float> qgptVect, qgetaVect, qgphiVect, ptVect, etaVect, phiVect, massVect;
  std::vector<int> qgpdgVect, pdgVect;
  std::vector<short> chgVect;
  
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;

    pytTree_p->GetEntry(entry%nEntriesPyt);
    if(doBkgd) bkgdTree_p->GetEntry(entry%nEntriesBkgd);

    pytHandler.Update();
    if(doBkgd) bkgdHandler.Update();

    qgptVect = pytHandler.GetQGPtVect();
    qgetaVect = pytHandler.GetQGEtaVect();
    qgphiVect = pytHandler.GetQGPhiVect();
    qgpdgVect = pytHandler.GetQGPdgVect();

    ptVect = pytHandler.GetPtVect();
    etaVect = pytHandler.GetEtaVect();
    phiVect = pytHandler.GetPhiVect();
    massVect = pytHandler.GetMassVect();
    pdgVect = pytHandler.GetPdgVect();
    chgVect = pytHandler.GetChgVect();

    nQGPart_ = 0;
    for(unsigned int pI = 0; pI < qgptVect.size(); ++pI){
      qgpt_[nPart_] = qgptVect[pI];
      qgeta_[nPart_] = qgetaVect[pI];
      qgphi_[nPart_] = qgphiVect[pI];
      qgpdgid_[nPart_] = qgpdgVect[pI];
      
      ++nPart_;
    }

    nPart_ = 0;
    for(unsigned int pI = 0; pI < ptVect.size(); ++pI){
      pt_[nPart_] = ptVect[pI];
      eta_[nPart_] = etaVect[pI];
      phi_[nPart_] = phiVect[pI];
      mass_[nPart_] = massVect[pI];
      pdgid_[nPart_] = pdgVect[pI];
      chg_[nPart_] = chgVect[pI];
      isSig_[nPart_] = true;
      
      ++nPart_;
    }

    if(doBkgd){
      ptVect = bkgdHandler.GetPtVect();
      etaVect = bkgdHandler.GetEtaVect();
      phiVect = bkgdHandler.GetPhiVect();
      massVect = bkgdHandler.GetMassVect();
      pdgVect = bkgdHandler.GetPdgVect();
      chgVect = bkgdHandler.GetChgVect();
      
      for(unsigned int pI = 0; pI < ptVect.size(); ++pI){
	pt_[nPart_] = ptVect[pI];
	eta_[nPart_] = etaVect[pI];
	phi_[nPart_] = phiVect[pI];
	mass_[nPart_] = massVect[pI];
	pdgid_[nPart_] = pdgVect[pI];
	chg_[nPart_] = chgVect[pI];
	isSig_[nPart_] = false;
	
	++nPart_;
      }
    }

    outTree_p->Fill();
  }

  pytHandler.Clean();

  if(doBkgd){
    bkgdHandler.Clean();
    inBkgdFile_p->Close();
    delete inBkgdFile_p;
  }
  
  inPytFile_p->Close();
  delete inPytFile_p;
  
  outFile_p->cd();

  outTree_p->Write("", TObject::kOverwrite);
  delete outTree_p;

  TDirectoryFile* paramDir_p = (TDirectoryFile*)outFile_p->mkdir("params");
  paramDir_p->cd();
  
  config.WriteToFile(outFile_p, paramDir_p);
  
  outFile_p->Close();
  delete outFile_p;

  config.Clean();
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/combinePythiaWithBkgd.exe <inConfigFileName>. return 1" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += combinePythiaWithBkgd(argv[1]);
  return retVal;
}
