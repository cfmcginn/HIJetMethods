//AUTHOR: Chris McGinn (2019.10.06)

//ROOT
#include "TDatabasePDG.h"
#include "TObjArray.h"
#include "TParticlePDG.h"

//JetMethods
#include "MakePytHyd/include/checkMakeDir.h"
#include "MakePytHyd/include/configParser.h"
#include "MakePytHyd/include/generalTreeHandler.h"
#include "MakePytHyd/include/stringUtil.h"

bool generalTreeHandler::Init(configParser inConfig, bool inIsPYT, TTree* inTree_p)
{
  Clean();
  if(inTree_p == nullptr){
    std::cout << "GENERALTREEHANDLER ERROR: given TTree on init is nullptr. init fails. return false" << std::endl;
    return false;
  }

  std::cout << "ARRAY: " << std::endl;
  TObjArray* tempObjArray = (TObjArray*)inTree_p->GetListOfBranches();
  for(Int_t oI = 0; oI < tempObjArray->GetEntries(); ++oI){
    treeBranches.push_back(tempObjArray->At(oI)->GetName());    
  }
  
  config = inConfig;
  isPYT = inIsPYT;

  std::string arrStr, nStr, ptStr, etaStr, phiStr, massStr, pdgStr, chgStr;
  std::string arrQGStr, nQGStr, qgptStr, qgetaStr, qgphiStr, qgpdgStr;
  
  if(isPYT){
    arrStr = config.GetParamValByName("PYTTYPE");
    if(isStrSame(arrStr, "ARR")) isARR = true;
    else if(isStrSame(arrStr, "VECT")) isARR = false;
    else{
      std::cout << "GENERALTREEHANDLER ERROR: PYTTYPE \'" << arrStr << "\' is not valid. Please pick \'VECT\' or \'ARR\'. return false" << std::endl;
      Clean();
      return false;
    }
    
    arrQGStr = config.GetParamValByName("PYTQGTYPE");
    if(isStrSame(arrQGStr, "ARR")) isQGARR = true;
    else if(isStrSame(arrQGStr, "VECT")) isQGARR = false;
    else{
      std::cout << "GENERALTREEHANDLER ERROR: PYTQGTYPE \'" << arrQGStr << "\' is not valid. Please pick \'VECT\' or \'QGARR\'. return false" << std::endl;
      Clean();
      return false;
    }

    isARR ? nStr = config.GetParamValByName("PYTN") : nStr = "";
    ptStr = config.GetParamValByName("PYTPT");
    etaStr = config.GetParamValByName("PYTETA");
    phiStr = config.GetParamValByName("PYTPHI");
    massStr = config.GetParamValByName("PYTMASS");
    pdgStr = config.GetParamValByName("PYTPDGID");
    chgStr = config.GetParamValByName("PYTCHG");    

    isQGARR ? nQGStr = config.GetParamValByName("PYTQGN") : nStr = "";
    qgptStr = config.GetParamValByName("PYTQGPT");
    qgetaStr = config.GetParamValByName("PYTQGETA");
    qgphiStr = config.GetParamValByName("PYTQGPHI");
    qgpdgStr = config.GetParamValByName("PYTQGPDGID");
  }
  else{
    arrStr = config.GetParamValByName("BKGDTYPE");
    if(isStrSame(arrStr, "ARR")) isARR = true;
    else if(isStrSame(arrStr, "VECT")) isARR = false;
    else{
      std::cout << "GENERALTREEHANDLER ERROR: BKGDTYPE \'" << arrStr << "\' is not valid. Please pick \'VECT\' or \'ARR\'. return false" << std::endl;
      Clean();
      return false;
    }
    
    isARR ? nStr = config.GetParamValByName("BKGDN") : nStr = "";
    ptStr = config.GetParamValByName("BKGDPT");
    etaStr = config.GetParamValByName("BKGDETA");
    phiStr = config.GetParamValByName("BKGDPHI");
    massStr = config.GetParamValByName("BKGDMASS");
    pdgStr = config.GetParamValByName("BKGDPDGID");
    chgStr = config.GetParamValByName("BKGDCHG");        
  }

  isGoodInit = true;
  bool essentials = IsBranchGood(ptStr) && IsBranchGood(etaStr) && IsBranchGood(phiStr) && IsBranchGood(pdgStr);
  if(isARR) essentials = essentials && IsBranchGood(nStr);
  
  if(!essentials){
    std::cout << "GENERALTREEHANDLER ERROR: Missing essentials branches. init failed. return false"  << std::endl;
    Clean();   
    return false;
  }

  hasChg = IsBranchGood(chgStr);
  hasMass = IsBranchGood(massStr);

  inTree_p->SetBranchStatus("*", 0);
  if(isARR) inTree_p->SetBranchStatus(nStr.c_str(), 1);
  inTree_p->SetBranchStatus(ptStr.c_str(), 1);
  inTree_p->SetBranchStatus(phiStr.c_str(), 1);
  inTree_p->SetBranchStatus(etaStr.c_str(), 1);
  inTree_p->SetBranchStatus(pdgStr.c_str(), 1);
  if(hasChg) inTree_p->SetBranchStatus(chgStr.c_str(), 1);
  if(hasMass) inTree_p->SetBranchStatus(massStr.c_str(), 1);

  if(isPYT){
    if(isQGARR) inTree_p->SetBranchStatus(nQGStr.c_str(), 1);
    inTree_p->SetBranchStatus(qgptStr.c_str(), 1);
    inTree_p->SetBranchStatus(qgphiStr.c_str(), 1);
    inTree_p->SetBranchStatus(qgetaStr.c_str(), 1);
    inTree_p->SetBranchStatus(qgpdgStr.c_str(), 1);

    if(isQGARR){
      qgptVect_p = new std::vector<float>;
      qgetaVect_p = new std::vector<float>;
      qgphiVect_p = new std::vector<float>;
      qgpdgVect_p = new std::vector<int>;
      
      inTree_p->SetBranchAddress(nQGStr.c_str(), &nQGPart_);
      inTree_p->SetBranchAddress(qgptStr.c_str(), qgptArr_);
      inTree_p->SetBranchAddress(qgphiStr.c_str(), qgphiArr_);
      inTree_p->SetBranchAddress(qgetaStr.c_str(), qgetaArr_);
      inTree_p->SetBranchAddress(qgpdgStr.c_str(), qgpdgArr_);
    }
    else{
      inTree_p->SetBranchAddress(qgptStr.c_str(), &qgptVect_p);
      inTree_p->SetBranchAddress(qgphiStr.c_str(), &qgphiVect_p);
      inTree_p->SetBranchAddress(qgetaStr.c_str(), &qgetaVect_p);
      inTree_p->SetBranchAddress(qgpdgStr.c_str(), &qgpdgVect_p);      
    }
  }
  
  if(isARR){
    ptVect_p = new std::vector<float>;
    etaVect_p = new std::vector<float>;
    phiVect_p = new std::vector<float>;
    massVect_p = new std::vector<float>;
    pdgVect_p = new std::vector<int>;
    chgVect_p = new std::vector<short>;

    inTree_p->SetBranchAddress(nStr.c_str(), &nPart_);
    inTree_p->SetBranchAddress(ptStr.c_str(), ptArr_);
    inTree_p->SetBranchAddress(phiStr.c_str(), phiArr_);
    inTree_p->SetBranchAddress(etaStr.c_str(), etaArr_);
    inTree_p->SetBranchAddress(pdgStr.c_str(), pdgArr_);
    if(hasChg) inTree_p->SetBranchAddress(chgStr.c_str(), chgArr_);
    if(hasMass) inTree_p->SetBranchAddress(massStr.c_str(), massArr_);
  }
  else{
    inTree_p->SetBranchAddress(ptStr.c_str(), &ptVect_p);
    inTree_p->SetBranchAddress(phiStr.c_str(), &phiVect_p);
    inTree_p->SetBranchAddress(etaStr.c_str(), &etaVect_p);
    inTree_p->SetBranchAddress(pdgStr.c_str(), &pdgVect_p);

    if(hasChg) inTree_p->SetBranchAddress(chgStr.c_str(), &chgVect_p);
    else chgVect_p = new std::vector<short>;

    if(hasMass) inTree_p->SetBranchAddress(massStr.c_str(), &massVect_p);
    else massVect_p = new std::vector<float>;
  }

  return true;
}

void generalTreeHandler::Update()
{
  if(!isGoodInit){
    std::cout << "GENERALTREEHANDLER ERROR: Update() called despite failed init. return" << std::endl;
    return;
  }

  TDatabasePDG* databasePDG_p = new TDatabasePDG();
  
  if(isARR){
    ptVect_p->clear();
    phiVect_p->clear();
    etaVect_p->clear();
    massVect_p->clear();
    pdgVect_p->clear();
    chgVect_p->clear();

    for(Int_t pI = 0; pI < nPart_; ++pI){
      ptVect_p->push_back(ptArr_[pI]);
      phiVect_p->push_back(phiArr_[pI]);
      etaVect_p->push_back(etaArr_[pI]);
      pdgVect_p->push_back(pdgArr_[pI]);      
      
      if(hasMass) massVect_p->push_back(massArr_[pI]);
      if(hasChg) chgVect_p->push_back(chgArr_[pI]);
    }
  }

  if(isPYT){
    if(isQGARR){
      ptVect_p->clear();
      phiVect_p->clear();
      etaVect_p->clear();
      pdgVect_p->clear();

      for(Int_t pI = 0; pI < nPart_; ++pI){
	ptVect_p->push_back(ptArr_[pI]);
	phiVect_p->push_back(phiArr_[pI]);
	etaVect_p->push_back(etaArr_[pI]);
	pdgVect_p->push_back(pdgArr_[pI]);
      }
    }
  }
  
  if(!hasMass){
    massVect_p->clear();
    for(unsigned int pI = 0; pI < pdgVect_p->size(); ++pI){
      TParticlePDG* tempPart_p = (TParticlePDG*)databasePDG_p->GetParticle(pdgVect_p->at(pI));
      massVect_p->push_back((float)tempPart_p->Mass());      
    }
  }

  if(!hasChg){
    chgVect_p->clear();
    for(unsigned int pI = 0; pI < pdgVect_p->size(); ++pI){
      TParticlePDG* tempPart_p = (TParticlePDG*)databasePDG_p->GetParticle(pdgVect_p->at(pI));
      chgVect_p->push_back((short)tempPart_p->Charge());      
    }
  }

  delete databasePDG_p;

  return;
}

void generalTreeHandler::Clean()
{
  if(nullptr != ptVect_p){
    ptVect_p->clear();
    delete ptVect_p;
    ptVect_p = nullptr;
  }
  if(nullptr != etaVect_p){
    etaVect_p->clear();
    delete etaVect_p;
    etaVect_p = nullptr;
  }
  if(nullptr != phiVect_p){
    phiVect_p->clear();
    delete phiVect_p;
    phiVect_p = nullptr;
  }
  if(nullptr != massVect_p){
    massVect_p->clear();
    delete massVect_p;
    massVect_p = nullptr;
  }
  if(nullptr != pdgVect_p){
    pdgVect_p->clear();
    delete pdgVect_p;
    pdgVect_p = nullptr;
  }
  if(nullptr != chgVect_p){
    chgVect_p->clear();
    delete chgVect_p;
    chgVect_p = nullptr;
  }
  
  if(nullptr != qgptVect_p){
    qgptVect_p->clear();
    delete qgptVect_p;
    qgptVect_p = nullptr;
  }
  if(nullptr != qgetaVect_p){
    qgetaVect_p->clear();
    delete qgetaVect_p;
    qgetaVect_p = nullptr;
  }
  if(nullptr != qgphiVect_p){
    qgphiVect_p->clear();
    delete qgphiVect_p;
    qgphiVect_p = nullptr;
  }
  if(nullptr != qgpdgVect_p){
    qgpdgVect_p->clear();
    delete qgpdgVect_p;
    qgpdgVect_p = nullptr;
  }

  hasMass = false;
  hasChg = false;
  nPart_ = 0;
  nQGPart_ = 0;
    
  isGoodInit = false;
  treeBranches.clear();
  config.Clean();  
  return;
}

bool generalTreeHandler::IsBranchGood(std::string inBranchName)
{
  if(!isGoodInit){
    std::cout << "GENERALTREEHANDLER ERROR: IsBranchGood is called despite failed init. return false" << std::endl;
    return false;
  } 
  if(inBranchName.size() == 0) return false;
  if(!vectContainsStr(inBranchName, &treeBranches)){
    std::cout << "GENERALTREEHANDLER ERROR: inBranchName \'" << inBranchName << "\' is not found. please check input file and tree. return false" << std::endl;
    return false;
  }
  
  return true;
}

