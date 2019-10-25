//Author: Chris McGinn (2019.10.16)

//cpp
#include <iostream>
#include <string>

//ROOT
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TTree.h"

//JetMethods
#include "EvalPytHyd/include/getLinBins.h"

#include "MakePytHyd/include/checkMakeDir.h"

#include "ProcessPytHyd/include/etaPhiFunc.h"

int makeVacuumAreaHist(std::string inFileName)
{
  if(!checkFileExt(inFileName, "root")) return 1;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  std::string outFileName = inFileName.substr(0, inFileName.rfind(".root"));
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  outFileName = "output/" + dateStr + "/" + outFileName + "_VacuumArea_" + dateStr + ".root";  
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1D* jetArea_p = new TH1D("jetArea_h", ";Area;Counts", 50, 0.0, 1.0);
  
  //we will adopt a (rough) CMS Geometry, which means 0.087x0.087 deltaEtaxdeltaPhi
  const Int_t nPhiBins = 72;
  Double_t phiBins[nPhiBins+1];
  getLinBins(-TMath::Pi(), TMath::Pi(), nPhiBins, phiBins);

  const Int_t nEtaBins = 44;
  Double_t etaBins[nEtaBins+1];
  const Double_t deltaEta = 0.087;
  getLinBins(-(nEtaBins*deltaEta)/2., (nEtaBins*deltaEta)/2., nEtaBins, etaBins);

  std::cout << "PhiBins: " << std::endl;
  for(Int_t pI = 0; pI < nPhiBins; ++pI){
    std::cout << " " << phiBins[pI] << ",";
  }
  std::cout << " " << phiBins[nPhiBins] << "." << std::endl;

  std::cout << "EtaBins: " << std::endl;
  for(Int_t pI = 0; pI < nEtaBins; ++pI){
    std::cout << " " << etaBins[pI] << ",";
  }
  std::cout << " " << etaBins[nEtaBins] << "." << std::endl;

  const Double_t deltaRJet = 0.4;
  const Double_t maxJtAbsEta = etaBins[nEtaBins] - deltaRJet;
  const Double_t minJtPt = 80;
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("processedTree_p");

  const Int_t nMaxParticles = 100000;
  Int_t nconstSig_;
  Float_t constptSig_[nMaxParticles];
  Float_t constphiSig_[nMaxParticles];
  Float_t constetaSig_[nMaxParticles];
  Int_t constjtposSig_[nMaxParticles];

  const Int_t nMaxJets = 500;
  Int_t njtSig_;
  Float_t jtptSig_[nMaxJets];
  Float_t jtphiSig_[nMaxJets];
  Float_t jtetaSig_[nMaxJets];

  inTree_p->SetBranchStatus("*", 0);
  inTree_p->SetBranchStatus("nconstSig", 1);
  inTree_p->SetBranchStatus("constptSig", 1);
  inTree_p->SetBranchStatus("constphiSig", 1);
  inTree_p->SetBranchStatus("constetaSig", 1);
  inTree_p->SetBranchStatus("constjtposSig", 1);

  inTree_p->SetBranchAddress("nconstSig", &nconstSig_);
  inTree_p->SetBranchAddress("constptSig", constptSig_);
  inTree_p->SetBranchAddress("constphiSig", constphiSig_);
  inTree_p->SetBranchAddress("constetaSig", constetaSig_);
  inTree_p->SetBranchAddress("constjtposSig", constjtposSig_);
  
  inTree_p->SetBranchStatus("njtSig", 1);
  inTree_p->SetBranchStatus("jtptSig", 1);
  inTree_p->SetBranchStatus("jtphiSig", 1);
  inTree_p->SetBranchStatus("jtetaSig", 1);

  inTree_p->SetBranchAddress("njtSig", &njtSig_);
  inTree_p->SetBranchAddress("jtptSig", jtptSig_);
  inTree_p->SetBranchAddress("jtphiSig", jtphiSig_);
  inTree_p->SetBranchAddress("jtetaSig", jtetaSig_);

  const Int_t nEntries = inTree_p->GetEntries();
  const Int_t nDiv = TMath::Max(1, nEntries/20);

  std::cout << "Processing " << nEntries << " events..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;   
    inTree_p->GetEntry(entry);

    Int_t areaCounter[nEtaBins][nPhiBins][nMaxJets];
    for(Int_t eI = 0; eI < nEtaBins; ++eI){
      for(Int_t pI = 0; pI < nPhiBins; ++pI){
	for(Int_t jI = 0; jI < njtSig_; ++jI){
	  areaCounter[eI][pI][jI] = 0;
	}
      }
    }

    for(Int_t cI = 0; cI < nconstSig_; ++cI){
      if(jtptSig_[constjtposSig_[cI]] < minJtPt) continue;
      if(TMath::Abs(jtetaSig_[constjtposSig_[cI]]) >= maxJtAbsEta) continue;

      Int_t etaPos = -1;
      Int_t phiPos = -1;
      for(Int_t eI = 0; eI < nEtaBins; ++eI){
	if(constetaSig_[cI] >= etaBins[eI] && constetaSig_[cI] < etaBins[eI+1]){
	  etaPos = eI;
	  break;
	}
      }

      for(Int_t pI = 0; pI < nPhiBins; ++pI){	  
	if(constphiSig_[cI] >= phiBins[pI] && constphiSig_[cI] < phiBins[pI+1]){
	  phiPos = pI;
	  break;
	}
      }


      ++(areaCounter[etaPos][phiPos][constjtposSig_[cI]]);
    }

    for(Int_t jI = 0; jI < njtSig_; ++jI){
      if(jtptSig_[jI] < minJtPt) continue;
      if(TMath::Abs(jtetaSig_[jI]) >= maxJtAbsEta) continue;

      Int_t areaTowers = 0;
      Int_t areaTowersFilled = 0;

      for(Int_t eI = 0; eI < nEtaBins; ++eI){
	Double_t etaCenter = (etaBins[eI] + etaBins[eI+1])/2.;
	
	for(Int_t pI = 0; pI < nPhiBins; ++pI){
	  Double_t phiCenter = (phiBins[pI] + phiBins[pI+1])/2.;

	  if(getDR(jtetaSig_[jI], jtphiSig_[jI], etaCenter, phiCenter) < deltaRJet){
	    ++areaTowers;
	    if(areaCounter[eI][pI][jI] >= 1) ++areaTowersFilled;
	  }
	}
      }

      Double_t area = ((Double_t)areaTowersFilled)/((Double_t)areaTowers);

      jetArea_p->Fill(area);
    }    
  }
  
  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  jetArea_p->Write("", TObject::kOverwrite);
  delete jetArea_p;
  
  outFile_p->Close();
  delete outFile_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/makeVacuumAreaHist.exe <inFileName>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += makeVacuumAreaHist(argv[1]);
  return retVal;
}
