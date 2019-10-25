//Author: Chris McGinn (2019.10.17)

//cpp
#include <iostream>
#include <string>

//ROOT
#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TNamed.h"
#include "TPad.h"
#include "TStyle.h"

//JetMethods
#include "EvalPytHyd/include/getLogBins.h"
#include "EvalPytHyd/include/histDefUtility.h"
#include "EvalPytHyd/include/plotUtilities.h"

#include "MakePytHyd/include/checkMakeDir.h"
#include "MakePytHyd/include/returnRootFileContentsList.h"

int plotJetStatTesting(const std::string inFileName)
{
  if(!checkFileExt(inFileName, "root")) return 1;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(16);
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> histList = returnRootFileContentsList(inFile_p, "TH1D");

  TNamed* nToyName_p = (TNamed*)inFile_p->Get("paramDir/nToys");
  TNamed* nPerToyName_p = (TNamed*)inFile_p->Get("paramDir/nPerToy");
  
  Int_t nBinsForPullSum = 0;
  Float_t ptLow = 100000;
  Float_t ptHigh = -1;
  for(unsigned int hI = 0; hI < histList.size(); ++hI){
    if(histList[hI].find("pull") == std::string::npos) continue;
    if(histList[hI].find("pullHist_h") != std::string::npos) continue;

    std::string ptLowStr = histList[hI].substr(histList[hI].find("Pt")+2, histList[hI].size());
    ptLowStr.replace(ptLowStr.find("_"), ptLowStr.size(), "");

    std::string ptHighStr = ptLowStr.substr(ptLowStr.find("to")+2, ptLowStr.size());
    ptHighStr.replace(ptHighStr.find("p"), 1, ".");

    ptLowStr.replace(ptLowStr.find("to"), ptLowStr.size(), "");
    ptLowStr.replace(ptLowStr.find("p"), 1, ".");

    if(ptLow > std::stod(ptLowStr)) ptLow = std::stod(ptLowStr);
    if(ptHigh < std::stod(ptHighStr)) ptHigh = std::stod(ptHighStr);
    
    ++nBinsForPullSum;
  }

  std::cout << "nPtBins: " << nBinsForPullSum << ", " << ptLow << "-" << ptHigh << std::endl;

  const Int_t nBinsMax = 100;
  Double_t bins[nBinsMax+1];
  getLogBins(ptLow, ptHigh, nBinsForPullSum, bins);
  
  TH1D* pullMean_p = new TH1D("pullMean_h", ";Particle p_{T} [GeV/c];#LTPull#GT", nBinsForPullSum, bins);
  TH1D* pullSigma_p = new TH1D("pullSigma_h", ";Particle p_{T} [GeV/c];#sigma(Pull)", nBinsForPullSum, bins);
  centerTitles({pullMean_p, pullSigma_p});

  //exampleToy1_Ex0_h
  for(unsigned int hI = 0; hI < histList.size(); ++hI){
    if(histList[hI].find("exampleToy") == std::string::npos) continue;
    if(histList[hI].find("exampleToy2") != std::string::npos) continue;

    std::string histName2 = histList[hI];
    histName2.replace(histName2.find("exampleToy1"), std::string("exampleToy1").size(), "exampleToy2");

    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    TH1D* hist1_p = (TH1D*)inFile_p->Get(histList[hI].c_str());
    TH1D* hist2_p = (TH1D*)inFile_p->Get(histName2.c_str());

    hist1_p->SetMarkerColor(1);
    hist1_p->SetMarkerSize(1);
    hist1_p->SetMarkerStyle(24);
    hist1_p->SetLineColor(1);
    
    hist2_p->SetMarkerColor(kRed);
    hist2_p->SetMarkerSize(1);
    hist2_p->SetMarkerStyle(25);
    hist2_p->SetLineColor(kRed);

    hist1_p->SetMaximum(1.1*TMath::Max(hist1_p->GetMaximum(), hist2_p->GetMaximum()));

    hist1_p->GetXaxis()->SetTitleOffset(1.25);
    
    hist1_p->DrawCopy("HIST E1 P");
    hist2_p->DrawCopy("HIST E1 P SAME");

    gStyle->SetOptStat(0);
    gPad->SetLogx();
    
    histName2.replace(histName2.find("Toy2"), 4, "");
    std::string saveName = "pdfDir/" + dateStr + "/" + histName2 + "_" + dateStr;
    quietSaveAs(canv_p, saveName + ".pdf");
    quietSaveAs(canv_p, saveName + ".png");
    delete canv_p;
  }

  
  for(unsigned int hI = 0; hI < histList.size(); ++hI){
    if(histList[hI].find("pull") == std::string::npos) continue;

    std::string ptStr = "PtInclusive";
    std::string ptStr2 = "Inclusive p_{T}";

    Float_t ptLow = 100000;
    Float_t ptHigh = -1;
    
    if(histList[hI].find("Pt") != std::string::npos){
      ptStr = histList[hI].substr(histList[hI].find("Pt"), histList[hI].size());
      ptStr.replace(ptStr.find("_"), ptStr.size(), "");
      
      ptStr2 = ptStr.substr(2, ptStr.size());
      while(ptStr2.find("p") != std::string::npos){
	ptStr2.replace(ptStr2.find("p"), 1, ".");
      }
      ptStr2.replace(ptStr2.find("to"), 2, " < p_{T} < ");
      ptStr2 = ptStr2 + " GeV/c";

      if(ptStr.find("p") != 4){
	ptStr.replace(0, 2, "Pt0");
      }
      
      std::string ptLowStr = histList[hI].substr(histList[hI].find("Pt")+2, histList[hI].size());
      ptLowStr.replace(ptLowStr.find("_"), ptLowStr.size(), "");
      
      std::string ptHighStr = ptLowStr.substr(ptLowStr.find("to")+2, ptLowStr.size());
      ptHighStr.replace(ptHighStr.find("p"), 1, ".");
      
      ptLowStr.replace(ptLowStr.find("to"), ptLowStr.size(), "");
      ptLowStr.replace(ptLowStr.find("p"), 1, ".");

      ptLow = std::stod(ptLowStr);
      ptHigh = std::stod(ptHighStr);
    }

    
    std::cout << "Plotting \'" << histList[hI] << "\'" << std::endl;

    TH1D* hist_p = (TH1D*)inFile_p->Get(histList[hI].c_str());
    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.01);
    canv_p->SetRightMargin(0.01);
    canv_p->SetLeftMargin(0.12);
    canv_p->SetBottomMargin(0.12);
    
    centerTitles(hist_p);

    hist_p->SetLineColor(1);
    hist_p->SetMarkerColor(1);
    hist_p->SetMarkerSize(1);
    hist_p->SetMarkerStyle(21);
    hist_p->DrawCopy("HIST E1 P");

    gStyle->SetOptStat(0);
    
    label_p->DrawLatex(0.18, 0.95, ptStr2.c_str());
    label_p->DrawLatex(0.18, 0.90, ("#LTPull#GT=" + prettyString(hist_p->GetMean(), 2, false) + "#pm" + prettyString(hist_p->GetMeanError(), 2, false)).c_str());
    label_p->DrawLatex(0.18, 0.85, ("#sigma(Pull)=" + prettyString(hist_p->GetStdDev(), 2, false) + "#pm" + prettyString(hist_p->GetStdDevError(), 2, false)).c_str());

    label_p->DrawLatex(0.7, 0.95, ("N_{Toys}=2x" + std::string(nToyName_p->GetTitle())).c_str());
    label_p->DrawLatex(0.7, 0.90, ("N_{Jet}/Toy=" + std::string(nPerToyName_p->GetTitle())).c_str());
    
    std::string saveName = "pdfDir/" + dateStr + "/pull_" + ptStr + "_" + dateStr;

    quietSaveAs(canv_p, saveName + ".pdf");
    quietSaveAs(canv_p, saveName + ".png");
    delete canv_p;

    if(ptHigh > 0){
      Double_t ptCent = (ptHigh + ptLow)/2.;

      Int_t binPos = -1;
      for(Int_t pI = 0; pI < nBinsForPullSum; ++pI){
	if(ptCent >= bins[pI] && ptCent < bins[pI+1]){
	  binPos = pI;
	}
      }

      pullMean_p->SetBinContent(binPos+1, hist_p->GetMean());
      pullMean_p->SetBinError(binPos+1, hist_p->GetMeanError());

      pullSigma_p->SetBinContent(binPos+1, hist_p->GetStdDev());
      pullSigma_p->SetBinError(binPos+1, hist_p->GetStdDevError());
    }
  }

  
  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(0.12);
  canv_p->SetBottomMargin(0.12);

  for(Int_t pI = 0; pI < nBinsForPullSum; ++pI){
    pullMean_p->SetBinContent(pI+1, pullMean_p->GetBinContent(pI+1) + 1.0);
  }

  pullMean_p->GetYaxis()->SetTitle("1 + #LTPull#GT or #sigma(Pull)");

  TLegend* leg_p = new TLegend(0.7, 0.8, 0.95, 0.95);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(16);
  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);

  pullMean_p->SetMarkerStyle(24);
  pullMean_p->SetMarkerColor(1);
  pullMean_p->SetLineColor(1);
  pullMean_p->SetMarkerSize(1);

  pullSigma_p->SetMarkerStyle(25);
  pullSigma_p->SetMarkerColor(kRed);
  pullSigma_p->SetLineColor(kRed);
  pullSigma_p->SetMarkerSize(1);

  pullMean_p->SetMaximum(1.45);
  pullMean_p->SetMinimum(0.55);

  pullMean_p->GetXaxis()->SetTitleFont(43);
  pullMean_p->GetYaxis()->SetTitleFont(43);
  pullMean_p->GetXaxis()->SetLabelFont(43);
  pullMean_p->GetYaxis()->SetLabelFont(43);

  pullMean_p->GetXaxis()->SetTitleSize(16);
  pullMean_p->GetYaxis()->SetTitleSize(16);
  pullMean_p->GetXaxis()->SetLabelSize(16);
  pullMean_p->GetYaxis()->SetLabelSize(16);

  pullMean_p->GetXaxis()->SetTitleOffset(1.5);
  
  pullMean_p->DrawCopy("HIST E1 P");
  pullSigma_p->DrawCopy("HIST E1 P SAME");

  leg_p->AddEntry(pullMean_p, "1 + #LTPull#GT", "P L");
  leg_p->AddEntry(pullSigma_p, "#sigma(Pull)", "P L");

  leg_p->Draw("SAME");

  label_p->DrawLatex(0.7, 0.25, ("N_{Toys}=2x" + std::string(nToyName_p->GetTitle())).c_str());
  label_p->DrawLatex(0.7, 0.20, ("N_{Jet}/Toy=" + std::string(nPerToyName_p->GetTitle())).c_str());

  
  TLine* line_p = new TLine();
  line_p->SetLineStyle(2);
  line_p->DrawLine(bins[0], 1.0, bins[nBinsForPullSum], 1.0);
  delete line_p;
  
  gPad->SetLogx();

  std::string saveName = "pdfDir/" + dateStr + "/pullMeanAdSigma_" + dateStr;
  quietSaveAs(canv_p, saveName + ".pdf");
  quietSaveAs(canv_p, saveName + ".png");
  delete canv_p;

  delete leg_p;
  
  delete pullMean_p;
  delete pullSigma_p;
  
  inFile_p->Close();
  delete inFile_p;

  delete label_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/plotJetStatTesting.exe <inFileName>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += plotJetStatTesting(argv[1]);
  return retVal;
}
