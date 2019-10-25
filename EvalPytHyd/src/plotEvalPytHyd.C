//Author: Chris McGinn (2019.10.22)

//cpp
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TCanvas.h"
#include "TDatime.h"
#include "TDirectoryFile.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TNamed.h"
#include "TPad.h"
#include "TStyle.h"

//Local
#include "EvalPytHyd/include/histDefUtility.h"
#include "EvalPytHyd/include/kirchnerPalette.h"
#include "EvalPytHyd/include/plotUtilities.h"

#include "MakePytHyd/include/checkMakeDir.h"
#include "MakePytHyd/include/returnRootFileContentsList.h"

std::vector<std::string> strVectValsFromStr(std::string inStr)
{
  std::vector<std::string> retVect;
  inStr = inStr + ",";
  while(inStr.find(",,") != std::string::npos){inStr.replace(inStr.find(",,"), 2, ",");}
  while(inStr.find(",") != std::string::npos){
    retVect.push_back(inStr.substr(0, inStr.find(",")));
    inStr.replace(0, inStr.find(",")+1, "");
  }
  return retVect;
}

std::vector<int> intVectValsFromStr(std::string inStr)
{
  std::vector<int> retVect;
  std::vector<std::string> tempVect = strVectValsFromStr(inStr);
  for(unsigned int vI = 0; vI < tempVect.size(); ++vI){
    retVect.push_back(std::stoi(tempVect[vI]));
  }
  return retVect;
}

std::vector<double> doubleVectValsFromStr(std::string inStr)
{
  std::vector<double> retVect;
  std::vector<std::string> tempVect = strVectValsFromStr(inStr);
  for(unsigned int vI = 0; vI < tempVect.size(); ++vI){
    retVect.push_back(std::stod(tempVect[vI]));
  }
  return retVect;
}

void getNXYPanel(Int_t nPanel, Int_t &nXPanel, Int_t &nYPanel)
{
  if(nPanel == 1){nXPanel = 1; nYPanel = 1;}
  else if(nPanel == 2){nXPanel = 2; nYPanel = 1;}
  else if(nPanel == 3){nXPanel = 3; nYPanel = 1;}
  else if(nPanel == 4){nXPanel = 4; nYPanel = 1;}
  else if(nPanel == 5){nXPanel = 3; nYPanel = 2;}
  else if(nPanel == 6){nXPanel = 3; nYPanel = 2;}
  else if(nPanel == 7){nXPanel = 4; nYPanel = 2;}
  else if(nPanel == 8){nXPanel = 4; nYPanel = 2;}
  else if(nPanel == 9){nXPanel = 3; nYPanel = 3;}
  else if(nPanel == 10){nXPanel = 5; nYPanel = 2;}
  else if(nPanel == 11){nXPanel = 4; nYPanel = 3;}
  else if(nPanel == 12){nXPanel = 4; nYPanel = 3;}
  else{
    std::cout << "GETNXYPANEL ERROR: Given nPanel \'" << nPanel << "\' does not have a valid x,y return vals. returning 0 0 for obvious failure in output. return" << std::endl;
    nXPanel = 0; nYPanel = 0;
  }
  
  return;
}

std::string labelHandler(std::string inStr)
{
  std::string cent = "Cent";
  std::string jt = "GenJtPt";
  if(inStr.find(cent) != std::string::npos){
    inStr.replace(inStr.find(cent), cent.size(), "");
    inStr.replace(inStr.find("to"), 2, "-");
    inStr = "#bf{" + inStr + "%}";
  }
  else if(inStr.find(jt) != std::string::npos){
    inStr.replace(inStr.find(jt), jt.size(), "");
    while(inStr.find("p") != std::string::npos){inStr.replace(inStr.find("p"), 1, ".");}
    inStr.replace(inStr.find("to"), 2, " < p_{T,Gen.} < ");
  }
  
  return inStr;
}

std::string algoToPrettyAlgo(std::string inStr)
{
  if(isStrSame(inStr, "SigBkgd")) inStr = "No Sub.";
  else if(isStrSame(inStr, "SigBkgdA")) inStr = "Area Sub.";
  else if(isStrSame(inStr, "SigBkgdSK")) inStr = "Softkiller";

  return inStr;
}

void plotScaleFromMap(std::map<std::string, std::map<std::string, std::vector<TH1D*> > > inMap, std::vector<std::string> tag1, std::vector<std::string> tag2, std::vector<std::string> tag3, std::string dateStr)
{
  const Int_t nStyles = 5;
  const Int_t styles[nStyles] = {24, 25, 28, 27, 46};

  const Int_t nColors = 4;
  const Int_t colors[nColors] = {0, 3, 6, 1};

  
  const Int_t nMaxPtBins = 100;  
  Int_t nPtBins = 0;
  Double_t ptBins[nMaxPtBins+1];
  for(unsigned pI = 0; pI < tag1.size(); ++pI){
    std::string frontStr = tag1[pI].substr(tag1[pI].find("Pt")+2, tag1[pI].size());
    frontStr.replace(frontStr.find("to"), frontStr.size(), "");
    std::string backStr = tag1[pI].substr(tag1[pI].find("to")+2, tag1[pI].size());
    while(frontStr.find("p") != std::string::npos){frontStr.replace(frontStr.find("p"), 1, ".");}
    while(backStr.find("p") != std::string::npos){backStr.replace(backStr.find("p"), 1, ".");}
    
    if(pI == 0) ptBins[nPtBins] = std::stod(frontStr);
    ++nPtBins;
    ptBins[nPtBins] = std::stod(backStr);
  }
  
  
  kirchnerPalette kPal;
  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(20);
  
  //Gonna do a little pre-proc
  const Int_t nMaxXPanel = 20;
  const Int_t nMaxYPanel = 20;
  const Double_t canvXWidth = 450.;
  const Double_t canvYWidth = canvXWidth;
  Int_t nPanel = tag1.size();
  Int_t nXPanel, nYPanel;
  getNXYPanel(nPanel, nXPanel, nYPanel);

  if(nXPanel > nMaxXPanel){
    std::cout << "PLOTSCALEFROMMAP ERROR: nXPanel requested \'" << nXPanel << "\' from total nPanel \'" << nPanel << "\' exceeds max allowed nMaxXPanel \'" << nMaxXPanel << "\'. return" << std::endl;
    return;
  }
  if(nYPanel > nMaxYPanel){
    std::cout << "PLOTSCALEFROMMAP ERROR: nYPanel requested \'" << nYPanel << "\' from total nPanel \'" << nPanel << "\' exceeds max allowed nMaxYPanel \'" << nMaxYPanel << "\'. return" << std::endl;
    return;
  }

  
  
  for(auto const & tag2Name : tag2){
    if(isStrSame(tag2Name, "SigBkgd")) continue;

    std::vector<std::vector<double> > mean, meanErr, stddev, stddevErr;
    
    TCanvas* canv_p = new TCanvas("canv_p", "", canvXWidth*nXPanel, canvYWidth*nYPanel);
    canv_p->SetTopMargin(0.01);
    canv_p->SetLeftMargin(0.01);
    canv_p->SetRightMargin(0.01);
    canv_p->SetBottomMargin(0.01);
    
    TLegend* leg_p = new TLegend(0.65, 0.6, 0.95, 0.85);
    leg_p->SetTextFont(43);
    leg_p->SetTextSize(20);
    leg_p->SetBorderSize(0);
    leg_p->SetFillColor(0);
    leg_p->SetFillStyle(0);

    TPad* pads_p[nMaxXPanel][nMaxYPanel];

    for(Int_t yI = 0; yI < nYPanel; ++yI){
      Double_t yLow = ((Double_t)(nYPanel - 1 - yI))/((Double_t)nYPanel);
      Double_t yHigh = ((Double_t)(nYPanel - yI))/((Double_t)nYPanel);

      for(Int_t xI = 0; xI < nXPanel; ++xI){
	Double_t xLow = ((Double_t)xI)/((Double_t)nXPanel);
	Double_t xHigh = ((Double_t)(xI+1))/((Double_t)nXPanel);

	canv_p->cd();
	pads_p[xI][yI] = new TPad(("pad_" + std::to_string(xI) + "_" + std::to_string(yI)).c_str(), "", xLow, yLow, xHigh, yHigh);
	pads_p[xI][yI]->SetTopMargin(0.01);
	pads_p[xI][yI]->SetRightMargin(0.01);
	pads_p[xI][yI]->SetBottomMargin(0.16);
	pads_p[xI][yI]->SetLeftMargin(0.16);
	
	pads_p[xI][yI]->Draw("SAME");
      }
    }

    for(unsigned int tI3 = 0; tI3 < tag3.size(); ++tI3){
      mean.push_back({});
      meanErr.push_back({});
      stddev.push_back({});
      stddevErr.push_back({});
    }
    
    for(unsigned int tI1 = 0; tI1 < tag1.size(); ++tI1){
      canv_p->cd();
      pads_p[tI1%nXPanel][tI1/nXPanel]->cd();
      
      Double_t padMax = -1;
      for(unsigned int tI3 = 0; tI3 < tag3.size(); ++tI3){
	if(isStrSame(tag3[tI3], "SigBkgd") && tag2Name.find("Cent50") == std::string::npos) continue;

	Double_t tempMax = getMax((inMap[tag1[tI1]])[tag2Name][tI3]);
	if(tempMax > padMax) padMax = tempMax;

	inMap[tag1[tI1]][tag2Name][tI3]->GetXaxis()->SetTitleFont(43);
	inMap[tag1[tI1]][tag2Name][tI3]->GetXaxis()->SetLabelFont(43);
	inMap[tag1[tI1]][tag2Name][tI3]->GetYaxis()->SetTitleFont(43);
	inMap[tag1[tI1]][tag2Name][tI3]->GetYaxis()->SetLabelFont(43);

	inMap[tag1[tI1]][tag2Name][tI3]->GetXaxis()->SetTitleSize(20);
	inMap[tag1[tI1]][tag2Name][tI3]->GetXaxis()->SetLabelSize(16);
	inMap[tag1[tI1]][tag2Name][tI3]->GetYaxis()->SetTitleSize(20);
	inMap[tag1[tI1]][tag2Name][tI3]->GetYaxis()->SetLabelSize(16);

	inMap[tag1[tI1]][tag2Name][tI3]->SetMarkerSize(1);
	inMap[tag1[tI1]][tag2Name][tI3]->SetMarkerStyle(styles[tI3%nStyles]);
	inMap[tag1[tI1]][tag2Name][tI3]->SetMarkerColor(kPal.getColor(colors[tI3%nColors]));
	inMap[tag1[tI1]][tag2Name][tI3]->SetLineColor(kPal.getColor(colors[tI3%nColors]));

	if(nYPanel == 2) inMap[tag1[tI1]][tag2Name][tI3]->GetXaxis()->SetTitleOffset(2.0);
      }

      bool isDrawn = false;
      for(unsigned int tI3 = 0; tI3 < tag3.size(); ++tI3){
	if(isStrSame(tag3[tI3], "SigBkgd") && tag2Name.find("Cent50") == std::string::npos) continue;
	
	(inMap[tag1[tI1]])[tag2Name][tI3]->SetMaximum(padMax*1.1);
	(inMap[tag1[tI1]])[tag2Name][tI3]->SetMinimum(0.0);
	
	if(isDrawn) (inMap[tag1[tI1]])[tag2Name][tI3]->DrawCopy("HIST E1 P SAME");
	else (inMap[tag1[tI1]])[tag2Name][tI3]->DrawCopy("HIST E1 P");

	TF1* fit_p = new TF1("fit_p", "gaus", 0.0, 3.0);
	fit_p->SetLineStyle(2);
	//	fit_p->SetLineColor(kPal.getColor(colors[tI3%nColors]));
	fit_p->SetLineColor(kPal.getColor(0));

	inMap[tag1[tI1]][tag2Name][tI3]->Fit("fit_p", "Q N M", "", 0.0, 3.0);
	fit_p->DrawCopy("SAME");

	mean[tI3].push_back(fit_p->GetParameter(1));
	meanErr[tI3].push_back(fit_p->GetParError(1));
	stddev[tI3].push_back(fit_p->GetParameter(2));
	stddevErr[tI3].push_back(fit_p->GetParError(2));
	
	delete fit_p;
	
	if(tI1 == 0){
	  leg_p->AddEntry((inMap[tag1[tI1]])[tag2Name][tI3], algoToPrettyAlgo(labelHandler(tag3[tI3])).c_str(), "P L");
	}
	
	isDrawn = true;
      }

      std::string labelName = labelHandler(tag1[tI1]);
      label_p->DrawLatex(0.2, 0.88, labelName.c_str());
      leg_p->Draw("SAME");
      gStyle->SetOptStat(0);
    }

    

    pads_p[0][0]->cd();
    std::string labelName = algoToPrettyAlgo(labelHandler(tag2Name));
    labelName = "#bf{PYTHIA+HYDJET, #color[" + std::to_string(kRed) + "]{" + labelName + "}}";
    label_p->DrawLatex(0.2, 0.95, labelName.c_str());      
    
    std::string saveName = "pdfDir/" + dateStr + "/scale_" + tag2Name + "_" + dateStr + ".pdf";
    quietSaveAs(canv_p, saveName);

    for(Int_t yI = 0; yI < nYPanel; ++yI){
      for(Int_t xI = 0; xI < nXPanel; ++xI){
	delete pads_p[xI][yI];
      }
    }

    delete canv_p;

    canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.01);
    canv_p->SetRightMargin(0.01);
    canv_p->SetBottomMargin(0.16);
    canv_p->SetLeftMargin(0.16);

    canv_p->cd();

    bool isDrawn = false;
    for(unsigned int tI3 = 0; tI3 < tag3.size(); ++tI3){
      if(isStrSame(tag3[tI3], "SigBkgd") && tag2Name.find("Cent50") == std::string::npos) continue;

      TH1D* histMean_p = new TH1D("histMean_p", ";Gen. Jet p_{T} [GeV/c];#LTReco./Gen.#GT", nPtBins, ptBins);
      centerTitles(histMean_p);

      histMean_p->GetXaxis()->SetTitleFont(43);
      histMean_p->GetXaxis()->SetLabelFont(43);
      histMean_p->GetYaxis()->SetTitleFont(43);
      histMean_p->GetYaxis()->SetLabelFont(43);
      
      histMean_p->GetXaxis()->SetTitleSize(20);
      histMean_p->GetXaxis()->SetLabelSize(16);
      histMean_p->GetYaxis()->SetTitleSize(20);
      histMean_p->GetYaxis()->SetLabelSize(16);
      
      histMean_p->SetMarkerSize(1);
      histMean_p->SetMarkerStyle(styles[tI3%nStyles]);
      histMean_p->SetMarkerColor(kPal.getColor(colors[tI3%nColors]));
      histMean_p->SetLineColor(kPal.getColor(colors[tI3%nColors]));
      
      for(unsigned int tI1 = 0; tI1 < tag1.size(); ++tI1){
	histMean_p->SetBinContent(tI1+1, mean[tI3][tI1]);
	histMean_p->SetBinError(tI1+1, meanErr[tI3][tI1]);
      }

      histMean_p->SetMaximum(1.35);
      histMean_p->SetMinimum(0.65);
      
      if(!isDrawn) histMean_p->DrawCopy("HIST E1 P");
      else histMean_p->DrawCopy("HIST E1 P SAME");

      isDrawn = true;

      delete histMean_p;
    }

    leg_p->SetY1NDC(0.2);
    leg_p->SetY2NDC(0.4);
    leg_p->Draw("SAME");
    gStyle->SetOptStat(0);

    label_p->DrawLatex(0.2, 0.95, labelName.c_str());

    saveName = "pdfDir/" + dateStr + "/scaleMean_" + tag2Name + "_" + dateStr + ".pdf";
    quietSaveAs(canv_p, saveName);
    delete canv_p;

    canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.01);
    canv_p->SetRightMargin(0.01);
    canv_p->SetBottomMargin(0.16);
    canv_p->SetLeftMargin(0.16);

    canv_p->cd();

    isDrawn = false;
    for(unsigned int tI3 = 0; tI3 < tag3.size(); ++tI3){
      if(isStrSame(tag3[tI3], "SigBkgd") && tag2Name.find("Cent50") == std::string::npos) continue;

      TH1D* histStdDev_p = new TH1D("histStdDev_p", ";Gen. Jet p_{T} [GeV/c];#sigma(Reco./Gen.)", nPtBins, ptBins);
      centerTitles(histStdDev_p);

      histStdDev_p->GetXaxis()->SetTitleFont(43);
      histStdDev_p->GetXaxis()->SetLabelFont(43);
      histStdDev_p->GetYaxis()->SetTitleFont(43);
      histStdDev_p->GetYaxis()->SetLabelFont(43);
      
      histStdDev_p->GetXaxis()->SetTitleSize(20);
      histStdDev_p->GetXaxis()->SetLabelSize(16);
      histStdDev_p->GetYaxis()->SetTitleSize(20);
      histStdDev_p->GetYaxis()->SetLabelSize(16);
      
      histStdDev_p->SetMarkerSize(1);
      histStdDev_p->SetMarkerStyle(styles[tI3%nStyles]);
      histStdDev_p->SetMarkerColor(kPal.getColor(colors[tI3%nColors]));
      histStdDev_p->SetLineColor(kPal.getColor(colors[tI3%nColors]));
      
      for(unsigned int tI1 = 0; tI1 < tag1.size(); ++tI1){
	histStdDev_p->SetBinContent(tI1+1, stddev[tI3][tI1]);
	histStdDev_p->SetBinError(tI1+1, stddevErr[tI3][tI1]);
      }

      histStdDev_p->SetMaximum(0.65);
      histStdDev_p->SetMinimum(0.0);
      
      if(!isDrawn) histStdDev_p->DrawCopy("HIST E1 P");
      else histStdDev_p->DrawCopy("HIST E1 P SAME");

      isDrawn = true;

      delete histStdDev_p;
    }

    leg_p->SetY1NDC(0.8);
    leg_p->SetY2NDC(0.94);
    leg_p->Draw("SAME");
    gStyle->SetOptStat(0);

    label_p->DrawLatex(0.2, 0.95, labelName.c_str());

    saveName = "pdfDir/" + dateStr + "/scaleStdDev_" + tag2Name + "_" + dateStr + ".pdf";
    quietSaveAs(canv_p, saveName);
    delete canv_p;

    canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.01);
    canv_p->SetRightMargin(0.01);
    canv_p->SetBottomMargin(0.16);
    canv_p->SetLeftMargin(0.16);

    canv_p->cd();

    isDrawn = false;
    for(unsigned int tI3 = 0; tI3 < tag3.size(); ++tI3){
      if(isStrSame(tag3[tI3], "SigBkgd") && tag2Name.find("Cent50") == std::string::npos) continue;

      TH1D* histStdDev_p = new TH1D("histStdDev_p", ";Gen. Jet p_{T} [GeV/c];#sigma(Reco./Gen.)/#mu", nPtBins, ptBins);
      centerTitles(histStdDev_p);

      histStdDev_p->GetXaxis()->SetTitleFont(43);
      histStdDev_p->GetXaxis()->SetLabelFont(43);
      histStdDev_p->GetYaxis()->SetTitleFont(43);
      histStdDev_p->GetYaxis()->SetLabelFont(43);
      
      histStdDev_p->GetXaxis()->SetTitleSize(20);
      histStdDev_p->GetXaxis()->SetLabelSize(16);
      histStdDev_p->GetYaxis()->SetTitleSize(20);
      histStdDev_p->GetYaxis()->SetLabelSize(16);
      
      histStdDev_p->SetMarkerSize(1);
      histStdDev_p->SetMarkerStyle(styles[tI3%nStyles]);
      histStdDev_p->SetMarkerColor(kPal.getColor(colors[tI3%nColors]));
      histStdDev_p->SetLineColor(kPal.getColor(colors[tI3%nColors]));
      
      for(unsigned int tI1 = 0; tI1 < tag1.size(); ++tI1){
	histStdDev_p->SetBinContent(tI1+1, stddev[tI3][tI1]/mean[tI3][tI1]);
	histStdDev_p->SetBinError(tI1+1, stddevErr[tI3][tI1]/mean[tI3][tI1]);
      }

      histStdDev_p->SetMaximum(0.65);
      histStdDev_p->SetMinimum(0.0);
      
      if(!isDrawn) histStdDev_p->DrawCopy("HIST E1 P");
      else histStdDev_p->DrawCopy("HIST E1 P SAME");

      isDrawn = true;

      delete histStdDev_p;
    }

    leg_p->SetY1NDC(0.8);
    leg_p->SetY2NDC(0.95);
    leg_p->Draw("SAME");
    gStyle->SetOptStat(0);

    label_p->DrawLatex(0.2, 0.95, labelName.c_str());

    saveName = "pdfDir/" + dateStr + "/scaleStdDev_" + tag2Name + "_" + dateStr + ".pdf";
    quietSaveAs(canv_p, saveName);
    delete canv_p;

    delete leg_p;
  }

  delete label_p;
  
  return;
}


int plotEvalPytHyd(std::string inFileName)
{
  if(!checkFileExt(inFileName, ".root")) return 1;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);
  
  //We will be plotting a ton of stuff, try to keep it modular
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> tdirList = returnRootFileContentsList(inFile_p, "TDirectoryFile");
  unsigned int pos = 0;
  while(pos < tdirList.size()){
    if(tdirList[pos].find("params") != std::string::npos) tdirList.erase(tdirList.begin()+pos);
    else ++pos;
  }

  std::vector<std::string> nameList = returnRootFileContentsList(inFile_p, "TNamed");
  std::vector<std::string> th1List = returnRootFileContentsList(inFile_p, "TH1D");
  std::vector<std::string> th2List = returnRootFileContentsList(inFile_p, "TH2D");

  //Declare the generic TH1* pointer vectors
  std::vector<TH1D*> th1s;
  std::vector<TH2D*> th2s;
  
  //param map
  std::map<std::string, std::string> paramMap;
  for(unsigned int nI = 0; nI < nameList.size(); ++nI){
    std::string name = nameList[nI].substr(nameList[nI].find("/")+1, nameList[nI].size());
    std::string title = ((TNamed*)inFile_p->Get(nameList[nI].c_str()))->GetTitle();
    paramMap[name] = title;
  }

  Int_t nCentBins = std::stoi(paramMap["NCENTBINS"]);
  std::vector<int> centBinsLow = intVectValsFromStr(paramMap["CENTBINSLOW"]);
  std::vector<int> centBinsHigh = intVectValsFromStr(paramMap["CENTBINSHIGH"]);
  std::vector<std::string> centBinsStr = strVectValsFromStr(paramMap["CENTBINSSTR"]);

  Int_t nJtPtBins = std::stoi(paramMap["NJTPTBINS"]);
  std::vector<int> jtPtBinsLow = intVectValsFromStr(paramMap["JTPTBINSLOW"]);
  std::vector<int> jtPtBinsHigh = intVectValsFromStr(paramMap["JTPTBINSHIGH"]);
  std::vector<std::string> jtPtBinsStr = strVectValsFromStr(paramMap["JTPTBINSSTR"]);

  std::cout << "NCENTBINS: " << nCentBins << std::endl;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::cout << " " << centBinsLow[cI] << "-" << centBinsHigh[cI] << std::endl;
  }

  std::cout << "NJTPTBINS: " << nJtPtBins << std::endl;
  for(Int_t jI = 0; jI < nJtPtBins; ++jI){
    std::cout << " " << jtPtBinsLow[jI] << "-" << jtPtBinsHigh[jI] << std::endl;
  }

  //We will use maps for simplified lookups
  std::map<std::string, std::map<std::string, std::vector<TH1D*> > > ptToAlgoToTH1Map, ptToCentToTH1Map;
  std::map<std::string, std::map<std::string, std::vector<TH2D*> > > ptToAlgoToTH2Map, ptToCentToTH2Map;
  //Build out the maps w/ empty vectors;
  for(unsigned int ptI = 0; ptI < jtPtBinsStr.size(); ++ptI){
    for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
      (ptToCentToTH1Map[jtPtBinsStr[ptI]])[centBinsStr[cI]] = {};
      (ptToCentToTH2Map[jtPtBinsStr[ptI]])[centBinsStr[cI]] = {};
    }

    for(unsigned int algoI = 0; algoI < tdirList.size(); ++algoI){
      (ptToAlgoToTH1Map[jtPtBinsStr[ptI]])[tdirList[algoI]] = {};
      (ptToAlgoToTH2Map[jtPtBinsStr[ptI]])[tdirList[algoI]] = {};
    }
  }

  
  //Lets do some pre-processing
  for(auto const & name : th1List){
    th1s.push_back((TH1D*)inFile_p->Get(name.c_str()));

    if(name.find("Scale_") != std::string::npos){
      th1s[th1s.size()-1]->Sumw2();
      th1s[th1s.size()-1]->Scale(1./th1s[th1s.size()-1]->Integral());
    }

    //Now we put it in the map -> first need ptStr, then cent or algo
    std::string histPtStr = "";
    for(unsigned int ptI = 0; ptI < jtPtBinsStr.size(); ++ptI){
      if(name.find("_" + jtPtBinsStr[ptI] + "_") != std::string::npos){
	histPtStr = jtPtBinsStr[ptI];
	break;
      }
    }

    std::string histCentStr = "";
    for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
      if(name.find("_" + centBinsStr[cI] + "_") != std::string::npos){
	histCentStr = centBinsStr[cI];
	break;
      }
    }

    std::string histAlgoStr = "";
    for(unsigned int algoI = 0; algoI < tdirList.size(); ++algoI){
      if(name.find(tdirList[algoI] + "/") != std::string::npos){
	histAlgoStr = tdirList[algoI];
	break;
      }
    }

    if(histPtStr.size() == 0) std::cout << "WARNING: COULD NOT FIND PT STRING IN HIST \'" << name << "\'" << std::endl;
    else{
      if(histCentStr.size() == 0) std::cout << "WARNING: COULD NOT FIND CENT STRING IN HIST \'" << name << "\'" << std::endl;
      else (ptToCentToTH1Map[histPtStr])[histCentStr].push_back(th1s[th1s.size()-1]);

      if(histAlgoStr.size() == 0) std::cout << "WARNING: COULD NOT FIND ALGO STRING IN HIST \'" << name << "\'" << std::endl;
      else (ptToAlgoToTH1Map[histPtStr])[histAlgoStr].push_back(th1s[th1s.size()-1]);
    }
  }
    
  //Validate our map building:

  std::cout << "Check maps: " << std::endl;
  for(auto const & iter : ptToCentToTH1Map){
    std::cout << " " << iter.first << std::endl;
    for(auto const & iter2 : iter.second){
      std::cout << "  " << iter2.first << std::endl;
      for(auto const & iter3 : iter2.second){
	std::cout << "   " << iter3->GetName() << std::endl;
      }
    }
  }
  
  plotScaleFromMap(ptToCentToTH1Map, jtPtBinsStr, centBinsStr, tdirList, dateStr);
  plotScaleFromMap(ptToAlgoToTH1Map, jtPtBinsStr, tdirList, centBinsStr, dateStr);

  
  inFile_p->Close();
  delete inFile_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/plotEvalPytHyd.exe <inFileName>. return 1" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += plotEvalPytHyd(argv[1]);
  return retVal;
}
