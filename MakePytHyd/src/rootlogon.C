{
  //For proper include paths w/ root interpreter code invocation
  gInterpreter->AddIncludePath(gSystem->ExpandPathName("$JETMETHODS"));
  gStyle->SetOptStat(0);
}
