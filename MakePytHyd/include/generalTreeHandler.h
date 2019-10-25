//AUTHOR: Chris McGinn (2019.10.06)

#ifndef GENERALTREEHANDLER_H
#define GENERALTREEHANDLER_H

//cpp
#include <string>
#include <vector>

//ROOT
#include "TTree.h"

//JetMethods
#include "MakePytHyd/src/configParser.C"

class generalTreeHandler{
 public:
  generalTreeHandler(){return;}
  generalTreeHandler(configParser inConfig, bool inIsPYT, TTree* inTree_p){Init(inConfig, inIsPYT, inTree_p); return;}
  ~generalTreeHandler(){return;}
  
  bool Init(configParser inConfig, bool inIsPYT, TTree* inTree_p);  
  void Update();
  std::vector<float> GetPtVect(){return *ptVect_p;}
  std::vector<float> GetPhiVect(){return *phiVect_p;}
  std::vector<float> GetEtaVect(){return *etaVect_p;}
  std::vector<float> GetMassVect(){return *massVect_p;}
  std::vector<int> GetPdgVect(){return *pdgVect_p;}
  std::vector<short> GetChgVect(){return *chgVect_p;}

  std::vector<float> GetQGPtVect(){return *qgptVect_p;}
  std::vector<float> GetQGPhiVect(){return *qgphiVect_p;}
  std::vector<float> GetQGEtaVect(){return *qgetaVect_p;}
  std::vector<int> GetQGPdgVect(){return *qgpdgVect_p;}
  void Clean();
  
 private:
  std::vector<std::string> treeBranches;
  bool isGoodInit = false;
  bool isPYT;
  bool isARR;
  bool isQGARR;
  configParser config;
  bool hasMass = false;
  bool hasChg = false;
  
  bool IsBranchGood(std::string inBranchName);

  std::vector<float>* ptVect_p = nullptr;
  std::vector<float>* phiVect_p = nullptr;
  std::vector<float>* etaVect_p = nullptr;
  std::vector<float>* massVect_p = nullptr;
  std::vector<int>* pdgVect_p = nullptr;
  std::vector<short>* chgVect_p = nullptr;

  std::vector<float>* qgptVect_p = nullptr;
  std::vector<float>* qgphiVect_p = nullptr;
  std::vector<float>* qgetaVect_p = nullptr;
  std::vector<int>* qgpdgVect_p = nullptr;

  static const Int_t nMaxPart = 100000;
  Int_t nPart_ = 0;
  Float_t ptArr_[nMaxPart];
  Float_t phiArr_[nMaxPart];
  Float_t etaArr_[nMaxPart];
  Float_t massArr_[nMaxPart];
  Int_t pdgArr_[nMaxPart];
  Short_t chgArr_[nMaxPart];

  Int_t nQGPart_ = 0;
  Float_t qgptArr_[nMaxPart];
  Float_t qgphiArr_[nMaxPart];
  Float_t qgetaArr_[nMaxPart];
  Int_t qgpdgArr_[nMaxPart];
};

#endif
