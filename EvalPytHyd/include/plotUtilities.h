#ifndef PLOTUTILITIES_H
#define PLOTUTILITIES_H

//cpp dependencies
#include <iostream>
#include <string>

//ROOT dependencies
#include "TBox.h"
#include "TCanvas.h"
#include "TError.h"
#include "TH1.h"
#include "TMath.h"

std::string prettyString(double inVal, const int prec, const bool doDot)
{
  std::string minStr = "";
  if(inVal < 0) minStr = "-";
  inVal = TMath::Abs(inVal);
  std::string retStr = "";
  if(inVal < 1.){
    inVal *= TMath::Power(10, prec);
    retStr = std::to_string(inVal);
    //  retStr = retStr.substr(0, 8) + "." + retStr.substr(8, retStr.size());

    for(int i = 0; i < prec; ++i){
      int pos = retStr.find(".");

      if(pos == 0) retStr = ".0" + retStr.substr(1, retStr.size());
      else retStr = retStr.substr(0, pos-1) + "." + retStr.substr(pos-1, 1) + retStr.substr(pos+1, retStr.size());
    }
  }
  else retStr = std::to_string(inVal);

  if(retStr.find(".") == std::string::npos) retStr = retStr + ".";
  
  while(retStr.size() < (unsigned int)(prec+1)){retStr = retStr + "0";}
  while(retStr.find(".") < retStr.size()-1-prec){retStr.replace(retStr.size()-1, 1,"");}
  if(doDot) retStr.replace(retStr.find("."), 1, "p");

  if(minStr.size() != 0) retStr = minStr + retStr.substr(0, retStr.size()-1);
  if(retStr.substr(0,1).find("p") != std::string::npos) retStr = "0" + retStr;
  if(retStr.substr(0,1).find(".") != std::string::npos) retStr = "0" + retStr;
  if(retStr.substr(0,2).find("-p") != std::string::npos) retStr = "-0" + retStr.substr(1, retStr.size());
  if(retStr.substr(0,2).find("-.") != std::string::npos) retStr = "-0" + retStr.substr(1, retStr.size());
  
  return retStr;
}

std::string prettyStringE(const double inVal, const int prec, const bool doDot)
{
  std::string retStr = prettyString(inVal, prec, false);
  int tenScale = retStr.find(".") - 1;
  while(retStr.find(".") != std::string::npos){
    retStr.replace(retStr.find("."), 1, "");
  }
  if(!doDot) retStr = retStr.substr(0,1) + "." + retStr.substr(1,retStr.size());
  else retStr = retStr.substr(0,1) + "p" + retStr.substr(1,retStr.size());

  while(retStr.find(".")+prec+1 < retStr.size()){
    retStr = retStr.substr(0, retStr.size()-1);
  }
  
  retStr = retStr + "E";
  if(tenScale >= 0) retStr = retStr + "+" + std::to_string(std::abs(tenScale));
  else retStr = retStr + "-" + std::to_string(std::abs(tenScale));
  
  return retStr;
}


void prettyCanv(TCanvas* canv_p)
{
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(1.5*canv_p->GetLeftMargin());
  canv_p->SetBottomMargin(canv_p->GetLeftMargin());
  canv_p->SetTopMargin(canv_p->GetLeftMargin()/2.);

  return;
}


void prettyTH1(TH1* hist_p, const double size, const int style, const int col)
{
  hist_p->SetMarkerSize(size);
  hist_p->SetMarkerStyle(style);
  hist_p->SetMarkerColor(col);
  hist_p->GetXaxis()->CenterTitle();
  hist_p->GetYaxis()->CenterTitle();

  return;
}


void drawWhiteBox(Double_t x1, Double_t x2, Double_t y1, Double_t y2)
{
  TBox* tempBox_p = new TBox();
  tempBox_p->SetFillColor(0);
  tempBox_p->DrawBox(x1, y1, x2, y2);
  delete tempBox_p;
}


void quietSaveAs(TCanvas* canv_p, const std::string saveName)
{
  Int_t oldLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kWarning;
  canv_p->SaveAs(saveName.c_str());
  gErrorIgnoreLevel = oldLevel;
  
  return;
}


double getNearestFactor10Up(double inVal, UInt_t steps = 0)
{
  double val = 0.00000000000001;
 
  for(UInt_t i = 0; i < 28; ++i){
    //    std::cout << val << std::endl;
    if(inVal >= val && inVal < val*10){
      val *= 10;
      break;
    }
    else val *= 10;
  }

  for(UInt_t i = 0; i < steps; ++i){
    val *= 10;
  }
  
  return val;
}

double getNearestFactor10Down(double inVal, UInt_t steps = 0)
{
  double val = 0.00000000000001;
 
  for(UInt_t i = 0; i < 28; ++i){
    //    std::cout << val << std::endl;
    if(inVal >= val && inVal < val*10){
      break;
    }
    else val *= 10;
  }

  for(UInt_t i = 0; i < steps; ++i){
    val /= 10;
  }
  
  return val;
}

#endif
