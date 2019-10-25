//Author: Chris McGinn (2019.10.05)

//cpp
#include <fstream>
#include <iostream>

//ROOT
#include "TNamed.h"

//JetMethods paths
#include "MakePytHyd/include/checkMakeDir.h"
#include "MakePytHyd/include/configParser.h"

bool configParser::Init(std::string inConfigName)
{
  Clean(); // just in case
  if(!checkFile(inConfigName)){
    std::cout << "CONFIGPARSER ERROR: Given inConfigName \'" << inConfigName << "\' is not valid file. configParser failed on init. return false" << std::endl;
    return false;
  }
  std::ifstream configFile(inConfigName.c_str());
  std::string tempStr;
  while(std::getline(configFile, tempStr)){
    while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, "");}
    if(tempStr.size() == 0) continue;
    if(tempStr.substr(0, 1).find("#") != std::string::npos) continue;
    while(tempStr.find("#") != std::string::npos){tempStr.replace(tempStr.rfind("#"), tempStr.size(), "");}
    
    std::cout << tempStr << std::endl;

    std::string frontStr = tempStr.substr(0, tempStr.find("="));
    std::string backStr = tempStr.substr(tempStr.find("=")+1, tempStr.size());

    paramNameToValMap[frontStr] = backStr;
  }
  configFile.close();
  isInitGood = true;
  return true;
}

std::string configParser::GetParamValByName(std::string name)
{
  bool mapHasStr = paramNameToValMap.count(name) > 0;
  if(!mapHasStr){
    std::cout << "CONFIGPARSER ERROR: name argument \'" << name << "\' is not found in map for GetParamValByName. returning \'\'" << std::endl;
    return "";
  }
  return paramNameToValMap[name];
}

bool configParser::WriteToFile(TFile* inFile_p, TDirectoryFile* inDir_p)
{
  if(inFile_p == nullptr){
    std::cout << "CONFIGPARSER ERROR: Given inFile to WriteToFile is null. return false" << std::endl;
    return false;
  }

  if(inDir_p == nullptr){
    inFile_p->cd();
    inDir_p = (TDirectoryFile*)inFile_p->mkdir("params");
  }

  inFile_p->cd();
  inDir_p->cd();

  for(auto const& iter : paramNameToValMap){
    TNamed tempName(iter.first.c_str(), iter.second.c_str());
    tempName.Write("", TObject::kOverwrite);
  }
  
  return true;
}

void configParser::Clean()
{
  isInitGood = false;
  paramNameToValMap.clear();
  
  return;
}

