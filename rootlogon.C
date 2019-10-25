{
  //For proper include paths w/ root interpreter code invocation
  gInterpreter->AddIncludePath(gSystem->ExpandPathName("$HIJETMETHODS"));
  gStyle->SetOptStat(0);
}
