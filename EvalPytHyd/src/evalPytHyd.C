//Author: Chris McGinn (2019.10.08)

//cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TDatime.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TTree.h"

//JetMethods
#include "EvalPytHyd/include/histDefUtility.h"
#include "EvalPytHyd/include/plotUtilities.h"

#include "MakePytHyd/include/checkMakeDir.h"
#include "MakePytHyd/include/returnRootFileContentsList.h"

int evalPytHyd(const std::string inFileName)
{
  if(!checkFileExt(inFileName, ".root")) return 1;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;
  
  std::string outFileName = inFileName.substr(0, inFileName.find(".root"));
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  outFileName = "output/" + dateStr + "/" + outFileName + "_EvalHist_" + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  //Some binnings
  const Int_t nMaxJetAlgo = 5;

  const Double_t maxJtAbsEta = 2.0;
  
  const Int_t nJtPtBins = 8;
  const Double_t jtPtBinsLow[nJtPtBins] = {30, 40, 50, 60, 80, 100, 120, 160};
  const Double_t jtPtBinsHigh[nJtPtBins] = {40, 50, 60, 80, 100, 120, 160, 200};
  const Int_t jtPtInterval = jtPtBinsHigh[nJtPtBins-1] - jtPtBinsLow[0];
  
  const Int_t nCentBins = 4;
  const Int_t centBinsLow[nCentBins] = {50, 30, 10, 0};
  const Int_t centBinsHigh[nCentBins] = {100, 50, 30, 10};

  TH1D* signalSpectra_p[nCentBins];
  TH1D* algoSignalSpectra_p[nMaxJetAlgo][nCentBins];
  TH1D* algoEff_p[nMaxJetAlgo][nCentBins];
  TH1D* algoScale_p[nMaxJetAlgo][nJtPtBins][nCentBins];
  TH2D* algoScaleVCent_p[nMaxJetAlgo][nJtPtBins];
  
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
    std::cout << "Given inFileName \'" << inFileName << "\' contains more than 1 ttree. fix. return 1" << std::endl;
    return 1;
  }
  TTree* inTree_p = (TTree*)inFile_p->Get(treeList[0].c_str());
  TObjArray* branchList_p = (TObjArray*)inTree_p->GetListOfBranches();
  std::vector<std::string> branchList;
  for(Int_t bI = 0; bI < branchList_p->GetEntries(); ++bI){
    branchList.push_back(branchList_p->At(bI)->GetName());
  }

  bool hasCent = vectContainsStr("cent", &branchList);
  bool hasNJtSig = vectContainsStr("njtSig", &branchList);
  bool hasOneOtherAlgo = false;
  std::vector<std::string> algoStr = {"Sig"};
  
  std::cout << "Branches: " << std::endl;
  if(hasNJtSig){
    for(unsigned int bI = 0; bI < branchList.size(); ++bI){
      std::cout << " "<< bI << "/" << branchList.size() << ": " << branchList[bI] << std::endl;
      
      if(branchList[bI].size() >= 3){
	if(branchList[bI].find("njt") != std::string::npos && !isStrSame("njtSig", branchList[bI])){
	  std::string endStr = branchList[bI].substr(3, branchList[bI].size());
	  if(!vectContainsStr("jtpt" + endStr, &branchList)) continue;
	  if(!vectContainsStr("jteta" + endStr, &branchList)) continue;
	  if(!vectContainsStr("jtphi" + endStr, &branchList)) continue;
	  if(!vectContainsStr("sigpos" + endStr, &branchList)) continue;
	  
	  algoStr.push_back(endStr);
	  hasOneOtherAlgo = true;
	}
      }    
    }
  }

  if(hasOneOtherAlgo){
    if(!vectContainsStr("jtptSig", &branchList)) hasNJtSig = false;
    if(!vectContainsStr("jtphiSig", &branchList)) hasNJtSig = false;
    if(!vectContainsStr("jtetaSig", &branchList)) hasNJtSig = false;

    for(unsigned int aI = 0; aI < algoStr.size(); ++aI){
      if(isStrSame(algoStr[aI], "Sig")) continue;
      if(!vectContainsStr("has" + algoStr[aI] + "MatchSig", &branchList)) hasNJtSig = false;
    }
  }
  
  if(!hasCent) std::cout << "Given inFileName \'" << inFileName << "\' ttree \'" << treeList[0] << "\' does not contain centrality var \'cent\'. return 1" << std::endl;
  if(!hasNJtSig) std::cout << "Given inFileName \'" << inFileName << "\' ttree \'" << treeList[0] << "\' does not contain signal jet collection, as denoted by var \'njtSig\'. return 1" << std::endl;
  if(!hasOneOtherAlgo) std::cout << "Given inFileName \'" << inFileName << "\' ttree \'" << treeList[0] << "\' does not contain any other valid algo, denoted by var \'njt*\', not equal to \'njtSig\'. return 1" << std::endl;

  if(!hasCent || !hasNJtSig || !hasOneOtherAlgo) return 1;

  Int_t cent_;

  Double_t minCent_ = inTree_p->GetMinimum("cent");
  Double_t maxCent_ = inTree_p->GetMaximum("cent");

  Double_t centFactor = 100./(maxCent_ - minCent_ + 1);
  
  const Int_t nMaxJets = 500;
  Int_t njt_[nMaxJetAlgo];
  Float_t jtpt_[nMaxJetAlgo][nMaxJets];
  Float_t jteta_[nMaxJetAlgo][nMaxJets];
  Float_t jtphi_[nMaxJetAlgo][nMaxJets];
  Int_t sigpos_[nMaxJetAlgo][nMaxJets];
  Bool_t hasMatch_[nMaxJetAlgo][nMaxJets];

  inTree_p->SetBranchStatus("*", 0);
  inTree_p->SetBranchStatus("cent", 1);
  inTree_p->SetBranchAddress("cent", &cent_);

  while(algoStr.size() > (unsigned int)nMaxJetAlgo){
    std::cout << "algos size \'" << algoStr.size() << "\' exceeds allowed max \'" << nMaxJetAlgo << "\'. removing \'" << algoStr[algoStr.size()-1] << "\'." << std::endl;
    algoStr.erase(algoStr.begin() + algoStr.size()-1);
  }

  outFile_p->cd();
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHigh[cI]);

    signalSpectra_p[cI] = new TH1D(("signalSpectra_" + centStr + "_h").c_str(), ";Gen. Jet p_{T} [GeV/c];Counts", jtPtInterval, jtPtBinsLow[0], jtPtBinsHigh[nJtPtBins-1]);
  }
  
  for(unsigned int aI = 0; aI < algoStr.size(); ++aI){
    for(Int_t pI = 0; pI < nJtPtBins; ++pI){
      std::string ptStr = "GenJtPt" + prettyString(jtPtBinsLow[pI], 1, true) + "to" + prettyString(jtPtBinsHigh[pI], 1, true);
      std::string ptStr2 =  prettyString(jtPtBinsLow[pI], 1, false) + "<p_{T}^{Gen.}<" + prettyString(jtPtBinsHigh[pI], 1, false);
      
      algoScaleVCent_p[aI][pI] = new TH2D(("algo" + algoStr[aI] + "ScaleVCent_" + ptStr + "_h").c_str(), (";Centrality (%);Reco./Gen. (" + ptStr2 + ")").c_str(), 100, -0.5, 99.5, 150, 0.0, 3.0);
      centerTitles(algoScaleVCent_p[aI][pI]);

      for(Int_t cI = 0; cI < nCentBins; ++cI){
	std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHigh[cI]);

	algoScale_p[aI][pI][cI] = new TH1D(("algo" + algoStr[aI] + "Scale_" + centStr + "_" + ptStr + "_h").c_str(), (";Reco./Gen. (" + ptStr2 + ");Counts").c_str(), 150, 0.0, 3.0);
	centerTitles(algoScale_p[aI][pI][cI]);

      }
    }

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHigh[cI]);
      algoSignalSpectra_p[aI][cI] = new TH1D(("algo" + algoStr[aI] + "SignalSpectra_" + centStr + "_h").c_str(), ";Gen. Jet p_{T} [GeV/c];Counts", jtPtInterval, jtPtBinsLow[0], jtPtBinsHigh[nJtPtBins-1]);
      algoEff_p[aI][cI] = new TH1D(("algo" + algoStr[aI] + "Eff_" + centStr + "_h").c_str(), ";Gen. Jet p_{T} [GeV/c];Efficiency", jtPtInterval, jtPtBinsLow[0], jtPtBinsHigh[nJtPtBins-1]);
      centerTitles({algoSignalSpectra_p[aI][cI], algoEff_p[aI][cI]});
    }
  }
  
  inFile_p->cd();  
  for(unsigned int aI = 0; aI < algoStr.size(); ++aI){
    inTree_p->SetBranchStatus(("njt" + algoStr[aI]).c_str(), 1);
    inTree_p->SetBranchStatus(("jtpt" + algoStr[aI]).c_str(), 1);
    inTree_p->SetBranchStatus(("jtphi" + algoStr[aI]).c_str(), 1);
    inTree_p->SetBranchStatus(("jteta" + algoStr[aI]).c_str(), 1);

    inTree_p->SetBranchAddress(("njt" + algoStr[aI]).c_str(), &njt_[aI]);
    inTree_p->SetBranchAddress(("jtpt" + algoStr[aI]).c_str(), jtpt_[aI]);
    inTree_p->SetBranchAddress(("jtphi" + algoStr[aI]).c_str(), jtphi_[aI]);
    inTree_p->SetBranchAddress(("jteta" + algoStr[aI]).c_str(), jteta_[aI]);

    if(!isStrSame(algoStr[aI], "Sig")){
      inTree_p->SetBranchStatus(("sigpos" + algoStr[aI]).c_str(), 1);
      inTree_p->SetBranchStatus(("has" + algoStr[aI] + "MatchSig").c_str(), 1);

      inTree_p->SetBranchAddress(("sigpos" + algoStr[aI]).c_str(), sigpos_[aI]);
      inTree_p->SetBranchAddress(("has" + algoStr[aI] + "MatchSig").c_str(), hasMatch_[aI]);
    }
  }

  const Int_t nEntries = inTree_p->GetEntries();
  const Int_t nDiv = TMath::Max(1, nEntries/20);

  std::cout << "Processing " << nEntries << " events..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;
    inTree_p->GetEntry(entry);

    Double_t cent2 = cent_*centFactor;
    Int_t centPos = -1;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(cent2 >= centBinsLow[cI] && cent2 < centBinsHigh[cI]){
	centPos = cI;
	break;
      }
    }

    if(centPos >= 0){
      for(Int_t jI = 0; jI < njt_[0]; ++jI){
	signalSpectra_p[centPos]->Fill(jtpt_[0][jI]);
      }
    }
    
    for(unsigned int aI = 0; aI < algoStr.size(); ++aI){
      if(!isStrSame("Sig", algoStr[aI])){

	for(Int_t jI = 0; jI < njt_[aI]; ++jI){
	  if(sigpos_[aI][jI] >= 0){
	    Double_t genEta = jteta_[0][sigpos_[aI][jI]];	    
	    if(TMath::Abs(genEta) >= maxJtAbsEta) continue;
	    Int_t ptPos = -1;
	    Double_t genPt = jtpt_[0][sigpos_[aI][jI]];	    
	    for(Int_t pI = 0; pI < nJtPtBins; ++pI){
	      if(genPt >= jtPtBinsLow[pI] && genPt < jtPtBinsHigh[pI]){
		ptPos = pI;
		break;
	      }
	    }	    	    
	    if(ptPos >= 0){
	      algoScaleVCent_p[aI][ptPos]->Fill(cent2, jtpt_[aI][jI]/genPt);
	      if(centPos >= 0) algoScale_p[aI][ptPos][centPos]->Fill(jtpt_[aI][jI]/genPt);
	    }
	  }
	}
      }
    }
  }
  
  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    signalSpectra_p[cI]->Write("", TObject::kOverwrite);
  }
  
  TDirectoryFile* dirs_p[nMaxJetAlgo];
  for(unsigned int aI = 0; aI < algoStr.size(); ++aI){
    if(!isStrSame(algoStr[aI], "Sig")){
      outFile_p->cd();
      dirs_p[aI] = (TDirectoryFile*)outFile_p->mkdir(algoStr[aI].c_str());
      dirs_p[aI]->cd();

      for(Int_t pI = 0; pI < nJtPtBins; ++pI){
	algoScaleVCent_p[aI][pI]->Write("", TObject::kOverwrite);
	for(Int_t cI = 0; cI < nCentBins; ++cI){
	  algoScale_p[aI][pI][cI]->Write("", TObject::kOverwrite);
	}
      }

      for(Int_t cI = 0; cI < nCentBins; ++cI){
	algoSignalSpectra_p[aI][cI]->Write("", TObject::kOverwrite);
	algoEff_p[aI][cI]->Write("", TObject::kOverwrite);
      }
    }

    for(Int_t pI = 0; pI < nJtPtBins; ++pI){
      delete algoScaleVCent_p[aI][pI];
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	delete algoScale_p[aI][pI][cI];
      }
    }

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      delete algoSignalSpectra_p[aI][cI];
      delete algoEff_p[aI][cI];
    }
  }

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    delete signalSpectra_p[cI];
  }

  outFile_p->cd();

  //We need to append some parameters
  outNameList.push_back(TNamed("MAXJTABSETA", prettyString(maxJtAbsEta, 2, false).c_str()));
  outNameList.push_back(TNamed("NJTPTBINS", std::to_string(nJtPtBins).c_str()));
  std::string jtPtBinsLowStr = "";
  std::string jtPtBinsHighStr = "";
  std::string jtPtBinsStr = "";

  for(Int_t jI = 0; jI < nJtPtBins; ++jI){
    std::string low = prettyString(jtPtBinsLow[jI], 1, false);
    std::string high = prettyString(jtPtBinsHigh[jI], 1, false);
    jtPtBinsLowStr = jtPtBinsLowStr + low + ",";
    jtPtBinsHighStr = jtPtBinsHighStr + high + ",";
    jtPtBinsStr = jtPtBinsStr + "GenJtPt" + low + "to" + high + ",";
  }
  while(jtPtBinsStr.find(".") != std::string::npos){jtPtBinsStr.replace(jtPtBinsStr.find("."), 1, "p");}
  
  outNameList.push_back(TNamed("JTPTBINSLOW", jtPtBinsLowStr.c_str()));
  outNameList.push_back(TNamed("JTPTBINSHIGH", jtPtBinsHighStr.c_str()));
  outNameList.push_back(TNamed("JTPTBINSSTR", jtPtBinsStr.c_str()));

  outNameList.push_back(TNamed("NCENTBINS", std::to_string(nCentBins).c_str()));
  std::string centBinsLowStr = "";
  std::string centBinsHighStr = "";
  std::string centBinsStr = "";
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    int low = centBinsLow[cI];
    int high = centBinsHigh[cI];
    centBinsLowStr = centBinsLowStr + std::to_string(low) + ",";
    centBinsHighStr = centBinsHighStr + std::to_string(high) + ","; 
    centBinsStr = centBinsStr + "Cent" + std::to_string(low) + "to" + std::to_string(high) + ",";
  }
  
  outNameList.push_back(TNamed("CENTBINSLOW", centBinsLowStr.c_str()));
  outNameList.push_back(TNamed("CENTBINSHIGH", centBinsHighStr.c_str()));
  outNameList.push_back(TNamed("CENTBINSSTR", centBinsStr.c_str()));
  
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
    std::cout << "Usage: ./bin/evalPytHyd.exe <inFileName>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += evalPytHyd(argv[1]);
  return retVal;
}
