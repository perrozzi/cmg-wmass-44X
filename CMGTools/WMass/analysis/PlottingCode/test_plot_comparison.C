void test_plot_comparison(){
  
  bool useUpDown = true;

  // TFile *fdata = new TFile("../JobOutputs/test_RecoilCorr_RochCorr_RecoilCorr_EffSFCorr_PtSFCorr_PileupSFCorr/test_numbers_DATA/WanalysisOnDATA.root");
  // TString outputfolder="scale/";
  // TFile *fmc = new TFile("../JobOutputs/test_RecoilCorr_RochCorr_RecoilCorr_EffSFCorr_PtSFCorr_PileupSFCorr/test_numbers_MCDATALIKE/WanalysisOnDATA.root");
  // TFile *fmcup = new TFile("../JobOutputs/test_RecoilCorr_RochCorr_RecoilCorr_Scale1_EffSFCorr_PtSFCorr_PileupSFCorr/test_numbers_MCDATALIKE/WanalysisOnDATA.root");
  // TFile *fmcdown = new TFile("../JobOutputs/test_RecoilCorr_RochCorr_RecoilCorr_Scale-1_EffSFCorr_PtSFCorr_PileupSFCorr/test_numbers_MCDATALIKE/WanalysisOnDATA.root");
  
  TFile *fdata = new TFile("../JobOutputs/test_RecoilCorr_RochCorr_RecoilCorr_EffSFCorr_PileupSFCorr/test_numbers_DATA/WanalysisOnDATA.root");
  TString outputfolder="resolution/";
  TFile *fmc = new TFile("../JobOutputs/test_RecoilCorr_RochCorr_RecoilCorr_EffSFCorr_PileupSFCorr/test_numbers_MCDATALIKE/WanalysisOnDATA.root");
  TFile *fmcup = new TFile("../JobOutputs/test_RecoilCorr_RochCorr_RecoilCorr_Resol1_EffSFCorr_PileupSFCorr/test_numbers_MCDATALIKE/WanalysisOnDATA.root");
  TFile *fmcdown = new TFile("../JobOutputs/test_RecoilCorr_RochCorr_RecoilCorr_Resol-1_EffSFCorr_PileupSFCorr/test_numbers_MCDATALIKE/WanalysisOnDATA.root");
  useUpDown = true;
  

  gROOT->ProcessLine(Form(".! mkdir %s",outputfolder.Data()));
  TIter nextkey(fdata->GetListOfKeys());
  TKey *key;
  while (key = (TKey*)nextkey()) {
    if(!key->ReadObj()->InheritsFrom("TH1D")) continue;
    TString hname = key->ReadObj()->GetName();
    // if(!hname.Contains("NonScaled_8_")) continue;
    
    // if(hname.Contains("hPileUp_Fall11_noWeights")) continue;
    if(hname.Contains("QCD")) continue;
    if(hname.Contains("METScaled")) continue;
    if(hname.Contains("MtScaled")) continue;
    if(hname.Contains("PtScaled")) continue;
    cout << "processing histo " << hname << endl;
    
    TCanvas*c1=new TCanvas("c1","c1",1200,800);
    TPad *pad1 = new TPad("pad1", "",0,0.25,1,1);
    // pad1->SetLogx(logx);
    // pad1->SetLogy(logy);
    // pad1->SetLogz(logz);
    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad1->SetBottomMargin(0);
    pad1->SetTopMargin(0.075);
    pad1->Draw();
    TPad *pad2 = new TPad("pad2", "",0,0,1,0.25);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.35);
    pad2->Draw();

    pad1->cd();  

    // TH1D *hmc=(TH1D*)key->ReadObj();
    TH1D *hmc=(TH1D*)fmc->Get(key->ReadObj()->GetName());
    TH1D *hmc=(TH1D*)key->ReadObj();
    hmc->Scale(1/hmc->Integral());
    hmc->SetLineColor(2);
    hmc->SetMarkerColor(2);
    hmc->DrawNormalized("histo");
    TH1D *hmcup, *hmcdown;
    if(useUpDown){
      hmcup=(TH1D*)fmcup->Get(key->ReadObj()->GetName());
      hmcup->Scale(1/hmcup->Integral());
      hmcup->SetLineColor(4);
      hmcup->SetMarkerColor(4);
      hmcup->DrawNormalized("same histo");
      hmcdown=(TH1D*)fmcdown->Get(key->ReadObj()->GetName());
      hmcdown->Scale(1/hmcdown->Integral());
      hmcdown->SetLineColor(6);
      hmcdown->SetMarkerColor(6);
      hmcdown->DrawNormalized("same histo");
    }
    TH1D *hdata=(TH1D*)fdata->Get(key->ReadObj()->GetName());
    hdata->Scale(1/hdata->Integral());
    hdata->SetLineColor(1);
    hdata->SetMarkerColor(1);
    hdata->DrawNormalized("same");
    
    TLatex *t = new TLatex();
    t->SetNDC();
    t->SetTextFont(63);
    t->SetTextSizePixels(17);  
    if(useUpDown) t->DrawLatex(0.65,0.75,Form("MC #downarrow : Mean= %.2f RMS= %.2f",hmcdown->GetMean(),hmcdown->GetRMS()));
    t->DrawLatex(0.65,0.725,Form("MC    : Mean= %.2f RMS= %.2f",hmc->GetMean(),hmc->GetRMS()));
    t->DrawLatex(0.65,0.7, Form("DATA : Mean= %.2f RMS= %.2f",hdata->GetMean(),hdata->GetRMS()));
    if(useUpDown) t->DrawLatex(0.65,0.675,  Form("MC #uparrow : Mean= %.2f RMS= %.2f",hmcup->GetMean(),hmcup->GetRMS()));
    
    pad2->cd();
    TH1D hmcdataratio=(*hmc)/(*hdata);
    hmcdataratio.SetLineColor(2);
    hmcdataratio.SetMarkerColor(2);
    hmcdataratio.GetYaxis()->SetRangeUser(0.9,1.1);
    hmcdataratio.GetYaxis()->SetLabelSize(0.075);
    hmcdataratio.GetXaxis()->SetLabelSize(0.1);
    hmcdataratio.Draw("histo");
    TH1D hmcupdataratio, hmcdowndataratio;
    if(useUpDown){
      hmcupdataratio=(*hmcup)/(*hdata);
      hmcupdataratio.SetLineColor(4);
      hmcupdataratio.SetMarkerColor(4);
      hmcupdataratio.Draw("same histo");
      hmcdowndataratio=(*hmcdown)/(*hdata);
      hmcdowndataratio.SetLineColor(6);
      hmcdowndataratio.SetMarkerColor(6);
      hmcdowndataratio.Draw("same histo");
    }
    
    c1->SaveAs(Form("%s/%s.png",outputfolder.Data(),key->ReadObj()->GetName()));
   
  }

}