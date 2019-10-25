//Author: Chris McGinn (2019.10.17)

//cpp
#include <iostream>
#include <string>

//ROOT
#include "TDatime.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TTree.h"

//JetMethods
#include "EvalPytHyd/include/getLogBins.h"
#include "EvalPytHyd/include/histDefUtility.h"
#include "EvalPytHyd/include/plotUtilities.h"

#include "MakePytHyd/include/checkMakeDir.h"

int jetStatTesting(const std::string inFileName)
{
  if(!checkFileExt(inFileName, "root")) return 1;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;
  
  std::string outFileName = inFileName.substr(0, inFileName.find(".root"));
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);

  outFileName = "output/" + dateStr + "/" + outFileName + "_JetStat_" + dateStr + ".root";

  const Int_t nToys = 90;
  const Double_t maxJtAbsEta = 2.0;
  const Double_t minJtPt = 60.;
  const Double_t minPartPt = 1.0;
  const Double_t maxPartPt = 60.;

  const Int_t nToyBins = 6;
  Double_t toyBins[nToyBins+1];
  getLogBins(minPartPt, maxPartPt, nToyBins, toyBins);
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1D* pullHist_p = new TH1D("pullHist_h", ";Pull;Counts", 51, -5.0, 5.0);
  TH1D* pullHist_PerPt_p[nToyBins];
  for(Int_t tI = 0; tI < nToyBins; ++tI){
    std::string ptStr = "Pt" + prettyString(toyBins[tI],1,true) + "to" + prettyString(toyBins[tI+1],1,true);
    pullHist_PerPt_p[tI] = new TH1D(("pullHist_" + ptStr + "_h").c_str(), ";Pull;Counts", 51, -5.0, 5.0);
  }  
  
  TH1D* toyHist1_p = new TH1D("toyHist1_h", ";Particle p_{T} [GeV/c];Counts", nToyBins, toyBins);
  TH1D* toyHist2_p = new TH1D("toyHist2_h", ";Particle p_{T} [GeV/c];Counts", nToyBins, toyBins);

  centerTitles({toyHist1_p, toyHist2_p});


  const Int_t exampleCap = 20;
  Int_t nExamples = 0;
  
  TH1D* globalHist1_p = new TH1D("globalHist1_h", ";Particle p_{T} [GeV/c];N_{particles}/Full Sample", nToyBins, toyBins);
  TH1D* globalHist2_p = new TH1D("globalHist2_h", ";Particle p_{T} [GeV/c];N_{particles}/FullSample", nToyBins, toyBins);
  TH1D* globalHist1_PerToy_p = new TH1D("globalHist1_PerToy_h", ";Particle p_{T} [GeV/c];N_{particles}/PYTHIA Ensemble", nToyBins, toyBins);
  TH1D* globalHist2_PerToy_p = new TH1D("globalHist2_PerToy_h", ";Particle p_{T} [GeV/c];N_{particles}/PYTHIA Ensemble", nToyBins, toyBins);
  TH1D* globalHist1_PerToyPerN_p = new TH1D("globalHist1_PerToyPerN_h", ";Particle p_{T} [GeV/c];N_{particles}/PYTHIA Jet", nToyBins, toyBins);
  TH1D* globalHist2_PerToyPerN_p = new TH1D("globalHist2_PerToyPerN_h", ";Particle p_{T} [GeV/c];N_{particles}/PYTHIA Jet", nToyBins, toyBins);

  std::vector<TH1*> histVect = {globalHist1_PerToy_p, globalHist2_PerToy_p, globalHist2_PerToy_p, globalHist2_PerToyPerN_p};
  setSumW2(histVect);
  histVect.push_back(globalHist1_p);
  histVect.push_back(globalHist2_p);
  centerTitles(histVect);

  
  Int_t toy1Counter = 0;
  Int_t toy2Counter = 0;
  Int_t global1Counter = 0;
  Int_t global2Counter = 0;
  
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

    std::vector<std::vector<double> > constPt;
    for(Int_t jI = 0; jI < njtSig_; ++jI){
      constPt.push_back({});
    }

    for(Int_t cI = 0; cI < nconstSig_; ++cI){
      constPt[constjtposSig_[cI]].push_back(constptSig_[cI]);
    }
    
    for(Int_t jI = 0; jI < njtSig_; ++jI){
      if(jtptSig_[jI] < minJtPt) continue;
      if(TMath::Abs(jtetaSig_[jI]) >= maxJtAbsEta) continue;

      if(toy1Counter >= nToys && toy2Counter >= nToys){
	for(Int_t tI = 0; tI < nToyBins; ++tI){
	  Double_t val = toyHist1_p->GetBinContent(tI+1) - toyHist2_p->GetBinContent(tI+1);
	  Double_t err = toyHist1_p->GetBinError(tI+1)*toyHist1_p->GetBinError(tI+1);
	  err += toyHist2_p->GetBinError(tI+1)*toyHist2_p->GetBinError(tI+1);
	  err = TMath::Sqrt(err);

	  pullHist_PerPt_p[tI]->Fill(val/err);
	  
	  pullHist_p->Fill(val/err);

	  for(Int_t bI = 0; bI < toyHist1_p->GetBinContent(tI+1); ++bI){
	    globalHist1_PerToy_p->Fill(toyHist1_p->GetBinCenter(tI+1));
	  }

	  for(Int_t bI = 0; bI < toyHist2_p->GetBinContent(tI+1); ++bI){
	    globalHist2_PerToy_p->Fill(toyHist2_p->GetBinCenter(tI+1));
	  }
	}

	if(nExamples < exampleCap){

	  outFile_p->cd();

	  toyHist1_p->Write(("exampleToy1_Ex" + std::to_string(nExamples) + "_h").c_str(), TObject::kOverwrite);
	  toyHist2_p->Write(("exampleToy2_Ex" + std::to_string(nExamples) + "_h").c_str(), TObject::kOverwrite);
	  
	  inFile_p->cd();	  
	  
	  ++nExamples;
	}
	
	delete toyHist1_p;
	delete toyHist2_p;

	outFile_p->cd();
	toyHist1_p = new TH1D("toyHist1_h", ";Particle p_{T} [GeV/c];Counts", nToyBins, toyBins);
	toyHist2_p = new TH1D("toyHist2_h", ";Particle p_{T} [GeV/c];Counts", nToyBins, toyBins);

	centerTitles({toyHist1_p, toyHist2_p});
	
	toy1Counter = 0;
	toy2Counter = 0;
      }
      
      if(toy1Counter < nToys){
	for(unsigned int cI = 0; cI < constPt[jI].size(); ++cI){
	  toyHist1_p->Fill(constPt[jI][cI]);
	  globalHist1_p->Fill(constPt[jI][cI]);
	  globalHist1_PerToyPerN_p->Fill(constPt[jI][cI]);
	}

	++global1Counter;
	++toy1Counter;
      }
      else if(toy2Counter < nToys){
	for(unsigned int cI = 0; cI < constPt[jI].size(); ++cI){
	  toyHist2_p->Fill(constPt[jI][cI]);
	  globalHist2_p->Fill(constPt[jI][cI]);
	  globalHist2_PerToyPerN_p->Fill(constPt[jI][cI]);
	}
	
	++global2Counter;
	++toy2Counter;
      }
    }    
  }
  
  
  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd(); 

  Int_t nTotalToys = pullHist_PerPt_p[0]->GetEntries();
  
  globalHist1_PerToy_p->Scale(1./(Double_t)nTotalToys);
  globalHist2_PerToy_p->Scale(1./(Double_t)nTotalToys);

  globalHist1_PerToyPerN_p->Scale(1./(Double_t)global1Counter);
  globalHist2_PerToyPerN_p->Scale(1./(Double_t)global2Counter);

  globalHist1_p->Write("", TObject::kOverwrite);
  globalHist2_p->Write("", TObject::kOverwrite);
  globalHist1_PerToy_p->Write("", TObject::kOverwrite);
  globalHist2_PerToy_p->Write("", TObject::kOverwrite);
  globalHist1_PerToyPerN_p->Write("", TObject::kOverwrite);
  globalHist2_PerToyPerN_p->Write("", TObject::kOverwrite);

  pullHist_p->Write("", TObject::kOverwrite);

  for(Int_t tI = 0; tI < nToyBins; ++tI){
    pullHist_PerPt_p[tI]->Write("", TObject::kOverwrite);
  }

  outFile_p->cd();
  TDirectoryFile* paramDir_p = (TDirectoryFile*)outFile_p->mkdir("paramDir");
  paramDir_p->cd();

  TNamed nToyName("nToys", std::to_string(nTotalToys));
  nToyName.Write("", TObject::kOverwrite);

  TNamed nPerToyName("nPerToy", std::to_string(nToys));
  nPerToyName.Write("", TObject::kOverwrite);
  
  paramDir_p->Close();
  delete paramDir_p;
  
  delete toyHist1_p;
  delete toyHist2_p;
  delete globalHist1_p;
  delete globalHist2_p;
  delete globalHist1_PerToy_p;
  delete globalHist2_PerToy_p;
  delete globalHist1_PerToyPerN_p;
  delete globalHist2_PerToyPerN_p;
  delete pullHist_p;

  for(Int_t tI = 0; tI < nToyBins; ++tI){
    delete pullHist_PerPt_p[tI];
  }
  
  outFile_p->Close();
  delete outFile_p;
  
  return 0; 
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/jetStatTesting <inFileName>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += jetStatTesting(argv[1]);
  return retVal;
}
