//Author: Chris McGinn (2019.10.05)

//cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

//JetMethods paths
#include "MakePytHyd/include/checkMakeDir.h"
#include "MakePytHyd/include/returnRootFileContentsList.h"
#include "MakePytHyd/include/stringUtil.h"

#include "MakePytHyd/include/configParser.h"
#include "MakePytHyd/include/generalTreeHandler.h"

int testConfigParser(std::string inConfigFileName)
{
  if(!checkFile(inConfigFileName)){
    std::cout << "Given inConfigFileName \'" << inConfigFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  configParser test(inConfigFileName);
  std::string pytFileName = test.GetParamValByName("INPYTFILENAME");
  std::string bkgdFileName = test.GetParamValByName("INBKGDFILENAME");
  std::string pytTreeName = test.GetParamValByName("INPYTTREENAME");
  std::string bkgdTreeName = test.GetParamValByName("INBKGDTREENAME");

  TFile* pytFile_p = new TFile(pytFileName.c_str(), "READ");
  std::vector<std::string> pytTreeNames = returnRootFileContentsList(pytFile_p, "TTree");
  if(!vectContainsStr(pytTreeName, &pytTreeNames)){
    std::cout << "Config file pytTreeName \'" << pytTreeName << "\' is not found. return 1" << std::endl;    
    pytFile_p->Close();
    delete pytFile_p;
    return 1;
  }

  TTree* pytTree_p = (TTree*)pytFile_p->Get(pytTreeName.c_str());
  generalTreeHandler pytHandler;
  if(!pytHandler.Init(test, true, pytTree_p)){
    std::cout << "INIT FAILED" << std::endl;
    pytFile_p->Close();
    delete pytFile_p;
    return 1;
  }

  const Int_t nEntries = pytTree_p->GetEntries();
  const Int_t nDiv = TMath::Max(1, nEntries/20);

  std::cout << "Processing " << nEntries << " events..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;

    pytTree_p->GetEntry(entry);
    pytHandler.Update();
    std::vector<float> ptVect_p = pytHandler.GetPtVect();
  }
  
  pytFile_p->Close();
  delete pytFile_p;  

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/testConfigParser.exe <inConfigFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += testConfigParser(argv[1]);
  return retVal;
}
