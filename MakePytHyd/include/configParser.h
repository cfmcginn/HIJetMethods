//Author: Chris McGinn (2019.10.05)

#ifndef CONFIGPARSER_H
#define CONFIGPARSER_H

//cpp
#include <map>
#include <string>

//ROOT
#include "TDirectoryFile.h"
#include "TFile.h"

class configParser{
 public:
  configParser(){}
  configParser(std::string inConfigName){Init(inConfigName); return;}
  ~configParser(){}

  bool Init(std::string inConfigName);
  std::string GetParamValByName(std::string name);
  bool WriteToFile(TFile* inFile_p, TDirectoryFile* inDir_p);
  void Clean();
  
 private:
  bool isInitGood = false;
  std::map <std::string, std::string> paramNameToValMap;
};

#endif
