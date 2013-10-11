#include "../includes/common.h"
#include "../AnalysisCode/common_stuff.C"

void PlotWdistributionsMCvsDATA_stack(TString folderMCsig="",TString folderMCEWK="",TString folderMCTT="",TString folderMCQCD="",TString folderDATA=""){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  
  cout << "PlotWdistributionsMCvsDATA_stack.C working dir: " << gSystem->WorkingDirectory() << endl;

  TFile *fMCsig=new TFile(Form("%s/WanalysisOnDATA.root",folderMCsig.Data()));
  TFile *fMCEWK=new TFile(Form("%s/WanalysisOnDATA.root",folderMCEWK.Data()));
  TFile *fMCTT=new TFile(Form("%s/WanalysisOnDATA.root",folderMCTT.Data()));
  TFile *fMCQCD=new TFile(Form("%s/WanalysisOnDATA.root",folderMCQCD.Data()));
  TFile *fDATA=new TFile(Form("%s/WanalysisOnDATA.root",folderDATA.Data()));
  
  fMCsig->Print();
  fMCEWK->Print();
  fMCTT->Print();
  fDATA->Print();

  TString LegendEvTypeTeX="W#rightarrow#mu#nu";
  
  for(int i=0; i<WMass::etaMuonNSteps; i++){
    TString eta_str = Form("%.1f",WMass::etaMaxMuons[i]); eta_str.ReplaceAll(".","p");
    // for(int j=0; j<2*WMass::WMassNSteps+1; j++){
    for(int j=WMass::WMassNSteps; j<WMass::WMassNSteps+1; j++){
      int jWmass = WMass::WMassCentral_MeV-(WMass::WMassNSteps-j)*WMass::WMassStep_MeV;
      
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hWPos_METNonScaled_8_JetCut_pdf10042-0_eta%s_%d",eta_str.Data(),jWmass),0,0,0,1,";MET [GeV];Counts",-1,-1,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hnvtx_Sig_eta%s_%d",eta_str.Data(),jWmass),0,0,0,1,";number of vertices;Counts",0,50,1);
      common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hnvtx_noWeights",eta_str.Data(),jWmass),0,0,0,1,";number of vertices;Counts",-1,-1,1);
      common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hnvtx_1_TrgAndGoodVtx",eta_str.Data(),jWmass),0,0,0,1,";number of vertices;Counts",-1,-1,1);
      common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hnvtx_5_RecoCut",eta_str.Data(),jWmass),0,0,0,1,";number of vertices;Counts",-1,-1,1);
      common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hnvtx_6_METCut",eta_str.Data(),jWmass),0,0,0,1,";number of vertices;Counts",-1,-1,1);
      common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hnvtx_7_RecoilCut",eta_str.Data(),jWmass),0,0,0,1,";number of vertices;Counts",-1,-1,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hWPos_pt_Sig_eta%s_%d",eta_str.Data(),jWmass),0,0,0,1,";number of vertices;Counts",-1,-1,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hnoTrgMuonsLeadingPt_Sig_eta%s_%d",eta_str.Data(),jWmass),0,0,0,1,";Z mass [GeV];Counts",50,-1,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hpfMET_WPos_Sig_eta%s_%d",eta_str.Data(),jWmass),0,1,0,1,";MET [GeV];Counts",-1,-1,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hpfMETphi_WPos_Sig_eta%s_%d",eta_str.Data(),jWmass),0,1,0,1,";MET #phi [rad];Counts",-1,-1,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hWPos_pt_unzoomed_Sig_eta%s_%d",eta_str.Data(),jWmass),0,1,0,1,";W p_{T} [GeV];Counts",-1,-1,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hWPos_pt_unzoomed_Sig_eta%s_%d",eta_str.Data(),jWmass),0,1,0,1,";W p_{T} [GeV];Counts",-1,100,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hWPos_pt_Sig_eta%s_%d",eta_str.Data(),jWmass),0,1,0,1,";W p_{T} [GeV];Counts",-1,-1,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hWPos_phi_Sig_eta%s_%d",eta_str.Data(),jWmass),0,1,0,1,";W #phi [rad];Counts",-1,-1,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hWPos_mt_Sig_eta%s_%d",eta_str.Data(),jWmass),0,1,0,1,";W m_{T} [GeV];Counts",-1,-1,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hWPos_mt_QCD_eta%s_%d",eta_str.Data(),jWmass),0,1,0,1,";W m_{T} QCD [GeV];Counts",-1,-1,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hMupt_WPos_Sig_eta%s_%d",eta_str.Data(),jWmass),0,1,0,1,";#mu p_{T} [GeV];Counts",-1,-1,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hMueta_WPos_Sig_eta%s_%d",eta_str.Data(),jWmass),0,1,0,1,";#mu #eta;Counts",-1,-1,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hMuphi_WPos_Sig_eta%s_%d",eta_str.Data(),jWmass),0,0,0,1,";#mu #phi [rad];Counts",-1,-1,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hMulogiso_WPos_Sig_eta%s_%d",eta_str.Data(),jWmass),0,0,0,1,";log rel iso;Counts",-1,-1,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hJetpt_WPos_Sig_eta%s_%d",eta_str.Data(),jWmass),0,1,0,1,";u1 (GeV);Counts",-1,-1,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hJeteta_WPos_Sig_eta%s_%d",eta_str.Data(),jWmass),0,1,0,1,";u2 (GeV);Counts",-1,-1,1);
      // common_stuff::plotAndSaveHisto1D_stack(LegendEvTypeTeX,fMCsig,fMCEWK,fMCTT,fMCQCD,fDATA,Form("hJetphi_WPos_Sig_eta%s_%d",eta_str.Data(),jWmass),0,1,0,1,";u2 (GeV);Counts",-1,-1,1);

    }

  }
  
  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

void plotAndSaveHisto1D(TFile*f1, TString str1, TFile*f2, TString str2, int logx, int logy, int logz, int scaleMCtoDATA){

  TH1D*h1=(TH1D*)f1->Get(str1.Data());
  h1->SetLineColor(2);
  TH1D*h2=(TH1D*)f2->Get(str2.Data());
  
  TCanvas *c1 = new TCanvas("c"+str1);
  c1->SetLogx(logx);
  c1->SetLogy(logy);
  c1->SetLogz(logz);
  
  if(scaleMCtoDATA){
    h1->DrawNormalized();
    h2->DrawNormalized("same");
  }else{
    h1->Draw();
    h2->Draw("same");    
  }
  
  c1->SaveAs(str1+".png");
  
}

