// --> >// $Id: rootlogon.C,v 1.1 2012/04/13 17:04:28 bendavid Exp $
{
 {
  TString libstr(Form("%s/lib/%s/%s",
                      gSystem->Getenv("CMSSW_BASE"),
                      gSystem->Getenv("SCRAM_ARCH"),
                      "libCondFormatsEgammaObjects.so"));

  gSystem->Load(libstr);
 }

  gSystem->AddIncludePath("-I$CMSSW_BASE/src/CondFormats/EgammaObjects/interface");

  gInterpreter->AddIncludePath((TString(":")+TString(gSystem->Getenv("CMSSW_BASE"))+
				TString("/src/CondFormats/EgammaObjects/interface")).Data());

}
