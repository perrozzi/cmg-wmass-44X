// UNCOMMENT TO USE PDF REWEIGHTING
//#define LHAPDF_ON

#ifdef LHAPDF_ON
  #include "LHAPDF/LHAPDF.h"
#endif 

#define Wanalysis_cxx
#include "Wanalysis.h"
#include "../includes/common.h"
#include "common_stuff.h"
#include "RecoilCorrector.h"
// #include "rochcor_42X.h"
#include "rochcor_44X_v3.h"
#include "MuScleFitCorrector.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include <vector>
#include <TGraphAsymmErrors.h>
#include <ctime>
#include <time.h>

void Wanalysis::Loop(int chunk, int Entry_ini, int Entry_fin, int IS_MC_CLOSURE_TEST, int isMCorDATA, TString outputdir, int useMomentumCorr, int smearRochCorrByNsigma, int useEffSF, int usePtSF, int useVtxSF, int controlplots, TString sampleName, int generated_PDF_set, int generated_PDF_member, int contains_PDF_reweight, int usePhiMETCorr, int useRecoilCorr)
{
  if (fChain == 0) return;

  std::map<std::string, TH1D*> h_1d;
  std::map<std::string, TH2D*> h_2d;

  cout << "generated_PDF_set= "<<generated_PDF_set
       << " generated_PDF_member= " << generated_PDF_member
       << " contains_PDF_reweight= " << contains_PDF_reweight
       << " WMass::NtoysMomCorr= " << WMass::NtoysMomCorr
       << endl;  

  TRandom3 *r = new TRandom3(0);
      
  #ifdef LHAPDF_ON
    // LHAPDF::initPDFSet(1,"CT10nnlo.LHgrid");
    if(!sampleName.Contains("DATA")){
      cout << "inizializing LHAPDF::initPDFSet(1)" << endl;
      LHAPDF::initPDFSet(1,generated_PDF_set,generated_PDF_member); // CMSSW DEFAULT
      cout << "finished inizializing LHAPDF" << endl;

      cout << "inizializing LHAPDF::initPDFSet(0)" << endl;
      // LHAPDF::initPDFSet();
      if(WMass::PDF_sets==11200)
        // LHAPDF::initPDFSet(0,"CT10nnlo.LHgrid");
        LHAPDF::initPDFSet(0,11200,0);
      else if(WMass::PDF_sets==232000)
        // LHAPDF::initPDFSet(0,"NNPDF23_nnlo_as_0118.LHgrid");
        LHAPDF::initPDFSet(0,232000,0);
      else if(WMass::PDF_sets==21200)
        // LHAPDF::initPDFSet(0,"MSTW2008nnlo68cl.LHgrid");
        LHAPDF::initPDFSet(0,21200,0);
      else if(WMass::PDF_sets==21241)
        LHAPDF::initPDFSet(0,21241,0);  // else if(WMass::PDF_sets<0)
        // LHAPDF::initPDFSet(0,generated_PDF_set,generated_PDF_member);
    }
    
  #endif 
  
  TFile*finEffSF, *finPileupSF,*finPtSF;
  TGraphAsymmErrors*hEffSF_MuId_eta_2011[2],*hEffSF_Iso_eta_2011[2],*hEffSF_HLT_eta_2011/* ,*hEffSF_Iso_vtx_2011A,*hEffSF_Iso_vtx_2011B*/;
  TH1D*hPileupSF,*hWPtSFPos,*hWPtSFNeg;

  // retrieve efficiencies SF
  if(useEffSF && (IS_MC_CLOSURE_TEST || isMCorDATA==0)){
    finEffSF = new TFile("../utils/MuonEfficiencies_SF_2011_44X_DataMC.root"); // used only to build templates
    hEffSF_MuId_eta_2011[0]=(TGraphAsymmErrors*)finEffSF->Get("SF_TIGHT_nL8_2011A_eta__pt>20");
    hEffSF_MuId_eta_2011[1]=(TGraphAsymmErrors*)finEffSF->Get("SF_TIGHT_nL8_2011B_eta__pt>20");
    hEffSF_Iso_eta_2011[0]=(TGraphAsymmErrors*)finEffSF->Get("combRelPFISO12_2011A_eta__pt>20");
    hEffSF_Iso_eta_2011[1]=(TGraphAsymmErrors*)finEffSF->Get("combRelPFISO12_2011B_eta__pt>20");
    hEffSF_HLT_eta_2011=(TGraphAsymmErrors*)finEffSF->Get("SF_HLT_MuIso24_2011_eta__pt>30");
    // hEffSF_Iso_vtx_2011A=(TH1D*)finEffSF->Get("combRelPFISO12_2011A_vtx__pt>20_abseta<2.4");
    // hEffSF_Iso_vtx_2011B=(TH1D*)finEffSF->Get("combRelPFISO12_2011B_vtx__pt>20_abseta<2.4");
    if(!finEffSF){
      cout << "file MuonEfficiencies_SF_2011_44X_DataMC.root is missing, impossible to retrieve efficiency scale factors" << endl;
      return;
    }
  }
  // retrieve pileup SF
  bool pileup_reweighting_npu = true;
  if(useVtxSF && (IS_MC_CLOSURE_TEST || isMCorDATA==0)){
    TString vtx_str = sampleName; vtx_str.ReplaceAll("Sig",""); vtx_str.ReplaceAll("Fake","");
    // finPileupSF = new TFile(Form("../utils/pileup_reweighting_%s.root",vtx_str.Data())); // used only to build templates
    finPileupSF = new TFile(Form("../utils/pileup/pileup_reweighting_Fall11.root")); // used only to build templates
    // finPileupSF = new TFile(Form("../utils/pileup/pileup_reweighting_Fall11_radoW.root")); pileup_reweighting_npu=false; // used only to build templates
    // finPileupSF = new TFile(Form("../utils/pileup/pileup_reweighting_Fall11_afterRecoilCut.root")); // used only to build templates
    if(!finPileupSF){
      cout << "file " << Form("../utils/pileup/pileup_reweighting_Fall11_afterRecoilCut.root") << " is missing, impossible to retrieve pileup reweighting factors" << endl;
      return;
    }else{
      hPileupSF=(TH1D*)finPileupSF->Get("hpileup_reweighting_Fall11");
    }
  }
  // retrieve pileup SF
  if(usePtSF && (IS_MC_CLOSURE_TEST || isMCorDATA==0)){
    // cout << "APPLYING W PT RESCALING" << endl;
    finPtSF = new TFile(Form("../utils/Wpt_reweighting.root")); // used only to build templates
    if(!finPtSF){
      cout << "file " << Form("../utils/Wpt_reweighting.root") << " is missing, impossible to retrieve W pt reweighting factors" << endl;
      // return;
    }else{
      hWPtSFPos=(TH1D*)finPtSF->Get("hWPos_pt_Sig_eta0p6");
      // hWPtSFNeg=(TH1D*)finPtSF->Get("hWNeg_pt_Sig_eta0p6");
    }
  }
  

  static const int nbins=75;
  double bins_scaled[3][nbins+1]={{0.}};
  double bins_Notscaled[3][nbins+1]={{0.}};
  double binsize1=0.01,binsize2=0.04;
  double binsize;
  double xmin=0.6,xmax=1.8, x;
  for(int k=0;k<3;k++){
    x=xmin;
    binsize=binsize1;
    for(int i=0;i<nbins;i++){
      bins_scaled[k][i]=x;
      bins_Notscaled[k][i]=x*80/(k==1 ? 1 : 2); // mT has double range wrt pt, met
      if(x>1.2-binsize) binsize=binsize2;
      x+=binsize;
      // cout << "bins_scaled["<<k<<"]["<<i<<"]= "<< bins_scaled[k][i] << endl;
      // cout << "bins_Notscaled["<<k<<"]["<<i<<"]= "<< bins_Notscaled[k][i] << endl;
    }
    bins_scaled[k][nbins]=xmax;
    bins_Notscaled[k][nbins]=xmax*80/(k==1 ? 1 : 2);
    // cout << "bins_scaled["<<k<<"]["<<nbins<<"]= " <<bins_scaled[k][nbins] << endl;
    // cout << "bins_Notscaled["<<k<<"]["<<nbins<<"]= " <<bins_Notscaled[k][nbins] << endl;
    // cout << endl;
  }
  
  // Long64_t first_entry = 0;
  // Long64_t nentries = fChain->GetEntriesFast();
  Long64_t first_entry = Entry_ini;
  Long64_t nentries = Entry_fin;
  // int chunk = nentries != fChain->GetEntries() && first_entry>0 ? (int)nentries/1e6:0;
  TString chunk_str = chunk>0? Form("_chunk%d",chunk) : "";
  ofstream outTXTfile;
  outTXTfile.open(Form("%s/Wanalysis_EVlog%s.log",outputdir.Data(),chunk_str.Data()));
  if(!outputdir.Contains("../")) outputdir = "../"+outputdir;
  cout << "output filename= " << Form("%s/Wanalysis%s.root",outputdir.Data(),chunk_str.Data()) << endl;

  // cout << "chunk " << chunk << endl;return;

  if(IS_MC_CLOSURE_TEST==1 && isMCorDATA==1) first_entry=nentries/2; // in case of closure test, DATA runs from N/2 to N
  if(IS_MC_CLOSURE_TEST==1 && isMCorDATA==0) nentries=nentries/2; // in case of closure test, MC runs from 0 to N/2
  if(IS_MC_CLOSURE_TEST==1) lumi_scaling=lumi_scaling*2; // in case of closure test, scaling must be multiplied by 2

  //To get the central value of the momentum correction
  int random_seed_start=67525;
  rochcor_44X_v3 *rmcor44X = WMass::NtoysMomCorr>1? new rochcor_44X_v3(random_seed_start) : new rochcor_44X_v3();  // make the pointer of rochcor class
  TString MuscleCard = (IS_MC_CLOSURE_TEST || isMCorDATA==0) ? "MuScleFit_2011_MC_44X" : "MuScleFit_2011_DATA_44X";
  TString fitParametersFile = MuscleCard+".txt";
  MuScleFitCorrector *corrector;
  if(useMomentumCorr==2){
    cout << "using MuscleFit card " << fitParametersFile << endl;
    corrector = new MuScleFitCorrector(fitParametersFile);
  }
  
  //////////////
  ////////   Initialize recoil corrections
  ////////

  std::string fileCorrectToPos = "../RecoilCode/recoilfit_genWpos_inc.root";
  std::string fileZmmData = "../RecoilCode/recoilfit_DATA_inc.root";
  std::string fileZmmMC = "../RecoilCode/recoilfit_genZ_inc.root";

  RecoilCorrector::RecoilCorrector*  correctorRecoilWPos_;

  correctorRecoilWPos_ = new RecoilCorrector(fileCorrectToPos.c_str(),123456); // this file is used to read the jet mutliplicity, will be a dummy file                      
  correctorRecoilWPos_->addDataFile(fileZmmData.c_str());
  correctorRecoilWPos_->addMCFile(fileZmmMC.c_str());

  // the following variables are dummy, but necessary to call the RecoilCorrector.
  double u1_dummy = 0;
  double u2_dummy = 0;
  double fluc_dummy = 0;
  double zero_dummy = 0;
  int jetMult = 0; // set to zero;
  //for the lepPt, lepPhi, 2: lepton is on leg2;
  
  
  // first_entry=469;
  // nentries=470;
  
  // start the actual event loop
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=first_entry; jentry<nentries;jentry++) {
  // for (Long64_t jentry=0; jentry<1e5;jentry++) { // for testing purposes
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    if(jentry%250000==0) cout <<"Analyzed entry "<<jentry<<"/"<<nentries<<endl;
    // if(jentry%500==0) cout <<"Analyzed entry "<<jentry<<"/"<<nentries<<endl;
    if(jentry%50000==0){
      time_t now = time(0);
      TString dt = ctime(&now); dt.ReplaceAll("\n"," ");
      outTXTfile << dt << "\t - \t Analyzed entry "<<jentry<<"/"<<nentries<<endl;
    }

    
    double evt_weight_original = lumi_scaling;
    if(useVtxSF && (IS_MC_CLOSURE_TEST || isMCorDATA==0) && npu>0) evt_weight_original=lumi_scaling*hPileupSF->GetBinContent(hPileupSF->GetXaxis()->FindBin(pileup_reweighting_npu?npu:nvtx));
    if(usePtSF && (IS_MC_CLOSURE_TEST || isMCorDATA==0) && sampleName.Contains("WJetsSig")){
      // cout << "TRYING TO APPLY W PT RESCALING" << endl;
      if(MuGen_charge>0||Mu_charge>0) evt_weight_original*=hWPtSFPos->GetBinContent(hWPtSFPos->GetXaxis()->FindBin(W_pt>0?W_pt:WGen_pt));
      // else if(MuGen_charge<0||Mu_charge<0) evt_weight_original*=hWPtSFNeg->GetBinContent(hWPtSFNeg->GetXaxis()->FindBin(W_pt>0?W_pt:WGen_pt));
    }
    
    int runopt = r->Rndm()<0.457451 ? 0 : 1;
    double Mu_tight_muon_SF = 1;
    // THIS MUST BE CHECKED !!!!!
    if(useEffSF && (IS_MC_CLOSURE_TEST || isMCorDATA==0)){
      Mu_tight_muon_SF = hEffSF_MuId_eta_2011[runopt]->Eval(Mu_eta)*hEffSF_Iso_eta_2011[runopt]->Eval(Mu_eta)*hEffSF_HLT_eta_2011->Eval(Mu_eta);
    }
    // cout << "Mu_tight_muon_SF= " << Mu_tight_muon_SF << endl;
    // if(!(IS_MC_CLOSURE_TEST || isMCorDATA==0) && run<175832) continue; // TEMPORARY TO TEST ROCHESTER CORRECTIONS
    
    // if((IS_MC_CLOSURE_TEST || isMCorDATA==0) && controlplots) hPileUp_Fall11->Fill(npu);
    if(controlplots){
      common_stuff::plot1D("hnvtx_noWeights", nvtx, 1, h_1d, 50,0,50 );
      if(IS_MC_CLOSURE_TEST || isMCorDATA==0)
          common_stuff::plot1D("hPileUp_Fall11_noWeights", npu, 1, h_1d, 50,0,50 );
    }
    
    for(int i=0; i<WMass::etaMuonNSteps; i++){
      TString eta_str = Form("%.1f",WMass::etaMaxMuons[i]); eta_str.ReplaceAll(".","p");
      
      for(int j=0; j<2*WMass::WMassNSteps+1; j++){
        int jWmass = (WMass::WMassCentral_MeV-(WMass::WMassNSteps-j)*WMass::WMassStep_MeV);
        
        if(!sampleName.Contains("WJetsSig") && WMass::WMassNSteps!=j) continue;
        double iWmass = (WMass::WMassCentral_MeV-(WMass::WMassNSteps-j)*WMass::WMassStep_MeV)/1e3;
        
        // cout << "events i= " << i << " j= " << j << " iWmass= " << iWmass<< " evt_weight_original= " << evt_weight_original << endl;
        
        // BW REWEIGHTING
        double shat=0,gamma=2.141 /*HARD CODED TO PDG VALUE*/,mw0=0,mw_i=0,weight_i=1;
        if(useGenVar){
          shat=WGen_m*WGen_m;
          mw0=WMass::WMassCentral_MeV/1e3;
          mw_i=iWmass;
          // ((shat - mw0^2)^2 + gamma^2 mw0^2) / ((shat - mw_i^2)^2 + gamma^2 mw_i^2)
          weight_i=(TMath::Power(shat - mw0*mw0,2) + TMath::Power(gamma*mw0,2)) / (TMath::Power(shat - mw_i*mw_i,2) + TMath::Power(gamma*mw_i,2));
          // cout << "WGen_m = " << WGen_m << " mw0= " << mw0 << " mw_i= " << mw_i << " weight_i= " << weight_i << endl;
        }
        double evt_weight=evt_weight_original*weight_i;
          
        // GEN STUFF IF REQUESTED
        if(useGenVar){
          
          if(WGen_m>0){ // only for signal event
            if(MuGen_charge>0){ // only for positive muons
              TString MuGenCharge_str = MuGen_charge>0? "Pos" : "Neg"; 
              
              double MuGen_var_jacobian[3] = {2*MuGen_pt/iWmass,WGen_mt/iWmass,2*NuGen_pt/iWmass,};
              // AVOID OVERFLOW BIN TO BE FILLED
              for(int k=0;k<3;k++)
                if(MuGen_var_jacobian[k]>=xmax) MuGen_var_jacobian[k]=xmax-binsize2/2;
              
              for(int k=0;k<3;k++){
                // hWPos_VarScaled_1_Gen[k][i][j]->Fill(MuGen_var_jacobian[k],evt_weight);
                common_stuff::plot1D(Form("hW%s_%sScaled_1_Gen_eta%s_%d",MuGenCharge_str.Data(),WMass::FitVar_str[k].Data(),eta_str.Data(),jWmass),
                                    MuGen_var_jacobian[k], evt_weight, h_1d, 
                                    // nbins, bins_scaled[k] );
                                    50, WMass::fit_xmin[k]/(WMass::WMassCentral_MeV/1e3),WMass::fit_xmax[k]/(WMass::WMassCentral_MeV/1e3) );
                // mass cut not meaningful for W case
                // hWPos_VarScaled_2_ZGenMassCut[k][i][j]->Fill(MuGen_var_jacobian[k],evt_weight);
                common_stuff::plot1D(Form("hW%s_%sScaled_2_ZGenMassCut_eta%s_%d",MuGenCharge_str.Data(),WMass::FitVar_str[k].Data(),eta_str.Data(),jWmass),
                                    MuGen_var_jacobian[k], evt_weight, h_1d, 
                                    // nbins, bins_scaled[k] );
                                    50, WMass::fit_xmin[k]/(WMass::WMassCentral_MeV/1e3),WMass::fit_xmax[k]/(WMass::WMassCentral_MeV/1e3) );
              }
                
              if(TMath::Abs(MuGen_eta)<WMass::etaMaxMuons[i]){
                for(int k=0;k<3;k++){
                  // hWPos_VarScaled_3_Mu1GenCut[k][i][j]->Fill(MuGen_var_jacobian[k],evt_weight);
                  common_stuff::plot1D(Form("hW%s_%sScaled_3_Mu1GenCut_eta%s_%d",MuGenCharge_str.Data(),WMass::FitVar_str[k].Data(),eta_str.Data(),jWmass),
                                  MuGen_var_jacobian[k], evt_weight, h_1d, 
                                  // nbins, bins_scaled[k] );
                                  50, WMass::fit_xmin[k]/(WMass::WMassCentral_MeV/1e3),WMass::fit_xmax[k]/(WMass::WMassCentral_MeV/1e3) );
                  // second lepton (i.e. neutrino) not detected
                  // if(TMath::Abs(MuNegGen_eta)<2.4){
                    // hWPos_VarScaled_4_Mu2GenCut[k][i][j]->Fill(MuGen_var_jacobian[k],evt_weight);
                    common_stuff::plot1D(Form("hW%s_%sScaled_4_Mu2GenCut_eta%s_%d",MuGenCharge_str.Data(),WMass::FitVar_str[k].Data(),eta_str.Data(),jWmass),
                                  MuGen_var_jacobian[k], evt_weight, h_1d, 
                                  // nbins, bins_scaled[k] );
                                  50, WMass::fit_xmin[k]/(WMass::WMassCentral_MeV/1e3),WMass::fit_xmax[k]/(WMass::WMassCentral_MeV/1e3) );
                  // }
                }
              }
            }
          }
        }
          
        if(!useGenVar || W_mt>0){ // dummy thing to separate between sig and bkg in W+Jets (useless)
          TString MuCharge_str = Mu_charge>0? "Pos" : "Neg"; 
          
          for(int m=0; m<WMass::NtoysMomCorr; m++){
            TString toys_str = "";
            if(WMass::NtoysMomCorr>1) toys_str = Form("_MomCorrToy%d",m);
            
            // good reco event selection
            if(evtHasGoodVtx && evtHasTrg /* && Mu_charge>0 */){
            
              if(controlplots && m==0 && i==0 && j==0 ){
                common_stuff::plot1D("hnvtx_1_TrgAndGoodVtx", nvtx, evt_weight, h_1d, 50,0,50 );
                if(IS_MC_CLOSURE_TEST || isMCorDATA==0)
                  common_stuff::plot1D("hPileUp_Fall11_Sig_1_TrgAndGoodVtx", npu, evt_weight, h_1d, 50,0,50 );
              }
              
              TLorentzVector mu; //TLorentzVector of the reconstructed muon
              //Set TLorentzVector of mu
              mu.SetPtEtaPhiM(Mu_pt,Mu_eta,Mu_phi,Mu_mass);
              // mu.Print();
              // cout << "before roch= " << m << " " << Mu_pt << " " << Mu_eta << " " << Mu_phi << " "<< Mu_mass << endl;
              //use rochester correction if required
              if(useMomentumCorr==1){ // use Rochester Momentum scale corrections if required
                if(IS_MC_CLOSURE_TEST || isMCorDATA==0){
                  // IN THE LAST VERSION NO RUN DEPENDENCE IN MC (TO be CHECKED)
                  // int runopt = r->Rndm()<0.457451 ? 0 : 1; // smear MC according to Run2011A and Run2011B statistics (if cut on pileup the 0.457... must be changed accordingly!!!
                  rmcor44X->momcor_mc(mu, Mu_charge, smearRochCorrByNsigma/* , runopt */);
                }
                else{
                  rmcor44X->momcor_data(mu, Mu_charge, smearRochCorrByNsigma, run<175832 ? 0 : 1 );
                }
              }else if(useMomentumCorr==2){ // use MuscleFit Momentum scale corrections if required
                corrector->applyPtCorrection(mu,Mu_charge);
              }
              
              //------------------------------------------------------------------------------------------------
              // Apply recoil corrections
              //------------------------------------------------------------------------------------------------

              // double pfmet_corr=nu.Pt();
              // double pfmetphi_corr=nu.Phi();

              // double pfmet_corr_minus = nu.Pt();
              // double pfmetphi_corr_minus = nu.Phi();

              // double pfmet_corr_plus = nu.Pt();
              // double pfmetphi_corr_plus = nu.Phi();

              // double u1_corr= 0;
              // double u2_corr = 0;

              // double u1_corr_plus = 0;
              // double u2_corr_plus = 0;

              // double u1_corr_minus = 0;
              // double u2_corr_minus = 0;

              int vtxBin=nvtx;
              if(nvtx==0) vtxBin=1;
              if(nvtx>20) vtxBin=20;

              // correctorRecoilWPos_->CorrectType1( pfmet_corr, pfmetphi_corr,
                  // WGen_pt, WGen_phi,
                  // mu.Pt(), mu.Phi(), 
                  // u1_corr, u2_corr, 0., 0.,
                  // vtxBin);

              // correctorRecoilWPos_->CorrectType1( pfmet_corr_Z_plus, pfmetphi_corr_Z_plus,
                  // WGen_pt, WGen_phi,
                  // mu.Pt(), mu.Phi(), 
                  // u1_corr_plus, u2_corr_plus, 1., 1.,
                  // vtxBin);

              // correctorRecoilWPos_->CorrectType1( pfmet_corr_Z_minus, pfmetphi_corr_Z_minus,
                  // ZGen_pt, ZGen_phi,
                  // mu.Pt(), mu.Phi(), 
                  // u1_corr_plus, u2_corr_plus, 1., 1.,
                  // vtxBin);

              // TLorentzVector nu_corr,W_corr; //TLorentzVector of the reconstructed W
              // nu_corr.SetPtEtaPhiM(pfmet_corr,0,pfmet_phi_corr,0);
              // W_corr = mu + nu_corr;
              
              
              
              if(useRecoilCorr==1 && sampleName.Contains("WJetsSig")){ // use Rochester Momentum scale corrections if required
                if(Mu_charge>0)
                  correctorRecoilWPos_->CorrectType1( pfmet, pfmet_phi,
                                                      WGen_pt, WGen_phi, 
                                                      mu.Pt(), mu.Phi(), 
                                                      u1_dummy, u2_dummy, 0, 0,
                                                      vtxBin);
              }
              // else if(useRecoilCorr==2){
                // correctorRecoilWPos_->CorrectType2( pfmet, pfmet_phi,
                                                      // ZGen_pt, ZGen_phi, 
                                                      // Z_pt, Z_phi, 
                                                      // u1_dummy, u2_dummy, fluc_dummy, zero_dummy,
                                                      // jetMult);
              // }else if(useRecoilCorr==3){
                // correctorRecoilWPos_->CorrectAll( pfmet, pfmet_phi,
                                                      // ZGen_pt, ZGen_phi, 
                                                      // Z_pt, Z_phi, 
                                                      // u1_dummy, u2_dummy, fluc_dummy, zero_dummy,
                                                      // jetMult);
              // }
              if(usePhiMETCorr==1){ // use Rochester Momentum scale corrections if required
                pair<double, double> pfmet_phicorr = common_stuff::getPhiCorrMET( pfmet, pfmet_phi, nvtx, !sampleName.Contains("DATA"));
                // cout << "0 pfmet= " << pfmet << " pfmet_phi= " << pfmet_phi << endl;
                pfmet = pfmet_phicorr.first;
                pfmet_phi = pfmet_phicorr.second;
                // cout << "1 pfmet= " << pfmet << " pfmet_phi= " << pfmet_phi << endl;
              }
              // cout << "after roch= " << m << " " << mu.Pt() << " " << mu.Eta() << " " << mu.Phi() << " "<< Mu_mass << endl;
              TLorentzVector nu,W; //TLorentzVector of the reconstructed W
              nu.SetPtEtaPhiM(pfmet,0,pfmet_phi,0);
              W = mu + nu;

              double Mu_var_jacobian[3] = {2*mu.Pt()/iWmass,W.Mt()/iWmass,2*nu.Pt()/iWmass}; // SCALED VARIABLE
              double Mu_var_NotScaled[3] = {mu.Pt(),W.Mt(),nu.Pt()}; // SCALED VARIABLE 
              // double Mu_var_NotScaled[3] = {mu.Pt(),W.Mt()*W.Mt(),nu.Pt()}; // SCALED VARIABLE LUCA TEMP!!!
              double dphiGen = MuGen_phi-NuGen_phi;
              if(dphiGen > TMath::Pi()) dphiGen -= 2*TMath::Pi();
              if(dphiGen < -TMath::Pi()) dphiGen += 2*TMath::Pi();
              double Mu_var_NotScaledGen[3] = {MuGen_pt,WGen_mt,dphiGen}; // SCALED VARIABLE
              // LUCA ADD TO AVOID OVERFLOW
              for(int k=0;k<3;k++)
                if(Mu_var_jacobian[k]>=xmax) Mu_var_jacobian[k]=xmax-binsize2/2;
              
              int wmass1 = iWmass*1e3;

              // good event with mu from W candidate within acceptance
              if( /* Z_mass>50  not possible to apply something similar to Z mass cut*/
                  TMath::Abs(mu.Eta())<WMass::etaMaxMuons[i] && mu.Pt()>30
                ){
                // muon candidate is passing tight, iso, dxy requirements
                if(MuIsTightAndIso && MuRelIso<0.12 && Mu_dxy<0.02 && noTrgMuonsLeadingPt<10 ){
                  
                  if(controlplots && m==0 && i==0 && j==0 ){
                    common_stuff::plot1D("hnvtx_5_RecoCut", nvtx, evt_weight, h_1d, 50,0,50 );
                    if(IS_MC_CLOSURE_TEST || isMCorDATA==0) 
                      common_stuff::plot1D("hPileUp_Fall11_5_RecoCut", npu, evt_weight, h_1d, 50,0,50 );
                  }

                  for(int k=0;k<3;k++)
                    // if(m==0) hWPos_VarScaled_5_RecoCut[k][i][j]->Fill(Mu_var_jacobian[k],evt_weight*Mu_tight_muon_SF);
                    if(m==0) common_stuff::plot1D(Form("hW%s_%sScaled_5_RecoCut_eta%s_%d",MuCharge_str.Data(),WMass::FitVar_str[k].Data(),eta_str.Data(),jWmass),
                                    Mu_var_jacobian[k],evt_weight*Mu_tight_muon_SF, h_1d,
                                    // nbins, bins_scaled[k] );
                                    50, WMass::fit_xmin[k]/(WMass::WMassCentral_MeV/1e3),WMass::fit_xmax[k]/(WMass::WMassCentral_MeV/1e3) );

                  if(pfmet>25){
                  
                    if(controlplots && m==0 && i==0 && j==0 ){
                      common_stuff::plot1D("hnvtx_6_METCut", nvtx, evt_weight, h_1d, 50,0,50 );
                      if(IS_MC_CLOSURE_TEST || isMCorDATA==0)
                        common_stuff::plot1D("hPileUp_Fall11_6_METCut", npu, evt_weight, h_1d, 50,0,50 );
                    }

                    for(int k=0;k<3;k++)
                      // if(m==0) hWPos_VarScaled_6_METCut[k][i][j]->Fill(Mu_var_jacobian[k],evt_weight*Mu_tight_muon_SF);
                      if(m==0) common_stuff::plot1D(Form("hW%s_%sScaled_6_METCut_eta%s_%d",MuCharge_str.Data(),WMass::FitVar_str[k].Data(),eta_str.Data(),jWmass),
                                    Mu_var_jacobian[k],evt_weight*Mu_tight_muon_SF, h_1d, 
                                    // nbins, bins_scaled[k] );
                                    50, WMass::fit_xmin[k]/(WMass::WMassCentral_MeV/1e3),WMass::fit_xmax[k]/(WMass::WMassCentral_MeV/1e3) );
                    
                    // if(W_pt<20){
                    if(W_pt<1e6){
                      for(int k=0;k<3;k++)
                        // if(m==0) hWPos_VarScaled_7_RecoilCut[k][i][j]->Fill(Mu_var_jacobian[k],evt_weight*Mu_tight_muon_SF);
                        if(m==0) common_stuff::plot1D(Form("hW%s_%sScaled_7_RecoilCut_eta%s_%d",MuCharge_str.Data(),WMass::FitVar_str[k].Data(),eta_str.Data(),jWmass),
                                    Mu_var_jacobian[k],evt_weight*Mu_tight_muon_SF, h_1d, nbins,bins_scaled[k] );
                      
                      // if(true){ // no jet pt cut at the moment
                      if(
                          Mu_var_NotScaled[0] > WMass::fit_xmin[0] && Mu_var_NotScaled[0] < WMass::fit_xmax[0]
                          && Mu_var_NotScaled[1] > WMass::fit_xmin[1] && Mu_var_NotScaled[1] < WMass::fit_xmax[1]
                          && Mu_var_NotScaled[2] > WMass::fit_xmin[2] && Mu_var_NotScaled[2] < WMass::fit_xmax[2]
                         ){
                      // if(Jet_leading_pt<30){
                        // cout << (mu.Pt()<xmax*80/2 ? mu.Pt() : (xmax-binsize2/2)*80/2 )<< endl;
                        
                        if(controlplots && m==0 && i==0 && j==0 ){
                          common_stuff::plot1D("hnvtx_7_RecoilCut", nvtx, evt_weight, h_1d, 50,0,50 );
                          if(IS_MC_CLOSURE_TEST || isMCorDATA==0)
                            common_stuff::plot1D("hPileUp_Fall11_7_RecoilCut", npu, evt_weight, h_1d, 50,0,50 );
                        }
                        
                                    // VERY DUMMY REWEIGHTING
                        // hWPos_VarNonScaled_8_JetCut[i][j]->Fill(mu.Pt()*iWmass/(WMass::WMassCentral_MeV/1e3)<xmax*80/2 ? mu.Pt()*iWmass/(WMass::WMassCentral_MeV/1e3) : (xmax-binsize2/2)*80/2 ,evt_weight*Mu_tight_muon_SF);
                        
                        // std::cout << "event= " << jentry << " mw0= " << mw0 << " iWmass= " << iWmass << " WGen_m= " << WGen_m << " weight_i= " << weight_i << std::endl;
                        // cout << "filling pt= " << (mu.Pt()*iWmass/(WMass::WMassCentral_MeV/1e3)<xmax*80/2 ? mu.Pt() : (xmax-binsize2/2)*80/2) <<" evt_weight= " << evt_weight << " Mu_tight_muon_SF= " << Mu_tight_muon_SF << endl;
                        double lha_weight = 1;
                        // double lha_weight = LHAPDF::xfx(0,x1,Q,fl1)*LHAPDF::xfx(0,x2,Q,fl2) / (LHAPDF::xfx(1,x1,Q,fl1)*LHAPDF::xfx(1,x2,Q,fl2));
                        double weight_old = 1;
                        #ifdef LHAPDF_ON
                          weight_old = !sampleName.Contains("DATA") ? (LHAPDF::xfx(1,parton1_x,scalePDF,parton1_pdgId)*LHAPDF::xfx(1,parton2_x,scalePDF,parton2_pdgId)) : 1;
                        #endif
                        if(m==0){
                          // hPDF_x1->Fill(TMath::Log10(parton1_x));
                          common_stuff::plot1D("hPDF_x1",
                                    TMath::Log10(parton1_x), 1, h_1d, 1000,-4,0 );
                          // hPDF_x1unweighted->Fill(TMath::Log10(parton1_x),1/weight_old);
                          common_stuff::plot1D("hPDF_x1unweighted",
                                    TMath::Log10(parton1_x),1/weight_old, h_1d, 1000,-4,0 );
                          // hPDF_x2->Fill(TMath::Log10(parton2_x));
                          common_stuff::plot1D("hPDF_x2",
                                    TMath::Log10(parton2_x), 1, h_1d, 1000,-4,0 );
                          // hPDF_x2unweighted->Fill(TMath::Log10(parton2_x),1/weight_old);
                          common_stuff::plot1D("hPDF_x2unweighted",
                                    TMath::Log10(parton2_x),1/weight_old, h_1d, 1000,-4,0 );
                        }
                        // cout << "scalePDF= " << scalePDF << " parton1_x= " << parton1_x << " parton1_pdgId= " << parton1_pdgId 
                        // << " parton2_x= " << parton2_x << " parton2_pdgId= " << parton2_pdgId << endl;
                        // cout << " LHAPDF::xfx(0,parton1_x,scalePDF,parton1_pdgId)= LHAPDF::xfx(0,"<<parton1_x<<","<<scalePDF<<","<<parton1_pdgId<<")= " << LHAPDF::xfx(0,parton1_x,scalePDF,parton1_pdgId) << endl;
                        // cout << " LHAPDF::xfx(0,parton2_x,scalePDF,parton2_pdgId)= LHAPDF::xfx(0,"<<parton2_x<<","<<scalePDF<<","<<parton2_pdgId<<")= " << LHAPDF::xfx(0,parton2_x,scalePDF,parton2_pdgId) << endl;
                        // cout << " LHAPDF::xfx(1,parton1_x,scalePDF,parton1_pdgId)= LHAPDF::xfx(1,"<<parton1_x<<","<<scalePDF<<","<<parton1_pdgId<<")= " << LHAPDF::xfx(1,parton1_x,scalePDF,parton1_pdgId) << endl;
                        // cout << " LHAPDF::xfx(1,parton2_x,scalePDF,parton2_pdgId)= LHAPDF::xfx(1,"<<parton2_x<<","<<scalePDF<<","<<parton2_pdgId<<")= " << LHAPDF::xfx(1,parton2_x,scalePDF,parton2_pdgId) << endl;

                        for(int h=0; h<WMass::PDF_members; h++){
                          if(!sampleName.Contains("DATA") && WMass::PDF_sets>0 && WMass::PDF_sets!=generated_PDF_set && WMass::PDF_members!=generated_PDF_member){
                            double weight_new = 1;
                            #ifdef LHAPDF_ON
                              LHAPDF::usePDFMember(0,h);
                              weight_new = (LHAPDF::xfx(0,parton1_x,scalePDF,parton1_pdgId)*LHAPDF::xfx(0,parton2_x,scalePDF,parton2_pdgId));
                            #endif
                            lha_weight = weight_new/weight_old;
                            // if(m==0) hPDF_weights[h]->Fill(lha_weight);
                            if(m==0) common_stuff::plot1D(Form("hPDF_weights_%d_%d",WMass::PDF_sets<0?generated_PDF_set:WMass::PDF_sets,h),
                                    lha_weight, 1, h_1d, 1000, 0, 2 );
                            cout << " lha_weight= " << lha_weight << endl;
                          }
                          for(int k=0;k<3;k++){
                            // hWPos_VarScaled_8_JetCut[m][h][k][i][j]->Fill(Mu_var_jacobian[k],evt_weight*Mu_tight_muon_SF*lha_weight);
                            common_stuff::plot1D(Form("hW%s_%sScaled_8_JetCut_pdf%d-%d%s_eta%s_%d",MuCharge_str.Data(),WMass::FitVar_str[k].Data(),WMass::PDF_sets<0?generated_PDF_set:WMass::PDF_sets,h,toys_str.Data(),eta_str.Data(),jWmass),
                                    Mu_var_jacobian[k],evt_weight*Mu_tight_muon_SF*lha_weight, h_1d, 
                                    // nbins, bins_scaled[k]  );
                                    50, WMass::fit_xmin[k]/(WMass::WMassCentral_MeV/1e3),WMass::fit_xmax[k]/(WMass::WMassCentral_MeV/1e3) );
                            // hWPos_VarNonScaled_8_JetCut[m][h][k][i][j]->Fill(Mu_var_NotScaled[k]*iWmass/(WMass::WMassCentral_MeV/1e3)<xmax*80/(k==1 ? 1 : 2) ? Mu_var_NotScaled[k] : (xmax-binsize2/2)*80/(k==1 ? 1 : 2) ,evt_weight*Mu_tight_muon_SF*lha_weight);
                            // TEST WITH GEN VARIABLE
                            // hWPos_VarNonScaled_8_JetCut[m][h][k][i][j]->Fill(Mu_var_NotScaled[k]*iWmass/(WMass::WMassCentral_MeV/1e3)<xmax*80/(k==1 ? 1 : 2) ? Mu_var_NotScaledGen[k] : (xmax-binsize2/2)*80/(k==1 ? 1 : 2) ,evt_weight*Mu_tight_muon_SF*lha_weight);
                            
                            common_stuff::plot1D(Form("hW%s_%sNonScaled_8_JetCut_pdf%d-%d%s_eta%s_%d",MuCharge_str.Data(),WMass::FitVar_str[k].Data(),WMass::PDF_sets<0?generated_PDF_set:WMass::PDF_sets,h,toys_str.Data(),eta_str.Data(),jWmass),
                                    // Mu_var_NotScaled[k]*iWmass/(WMass::WMassCentral_MeV/1e3)<xmax*80/(k==1 ? 1 : 2) ? Mu_var_NotScaled[k] : (xmax-binsize2/2)*80/(k==1 ? 1 : 2) ,evt_weight*Mu_tight_muon_SF*lha_weight, h_1d, 
                                    Mu_var_NotScaled[k] ,evt_weight*Mu_tight_muon_SF*lha_weight, h_1d, 
                                    // nbins, bins_Notscaled[k] );
                                    50, WMass::fit_xmin[k],WMass::fit_xmax[k] );
                          }
                          if(controlplots && m==0 && i==0 && j==0 ){
                            common_stuff::plot1D("hnvtx_8_JetCut", nvtx, evt_weight, h_1d, 50,0,50 );
                            if(IS_MC_CLOSURE_TEST || isMCorDATA==0)
                              common_stuff::plot1D("hPileUp_Fall11_8_JetCut", npu, evt_weight, h_1d, 50,0,50 );
                          }

                          common_stuff::plot2D(Form("hW%s_%sVs%sScaled_8_JetCut_pdf%d-%d%s_eta%s_%d",MuCharge_str.Data(),WMass::FitVar_str[1].Data(),WMass::FitVar_str[0].Data(),WMass::PDF_sets<0?generated_PDF_set:WMass::PDF_sets,h,toys_str.Data(),eta_str.Data(),jWmass),
                                  Mu_var_NotScaled[0],Mu_var_NotScaled[1],evt_weight*Mu_tight_muon_SF*lha_weight, h_2d, 
                                  // nbins, bins_Notscaled[0], nbins, bins_Notscaled[1]  );
                                  50, WMass::fit_xmin[0], WMass::fit_xmax[0], 50, WMass::fit_xmin[1], WMass::fit_xmax[1] );
                          common_stuff::plot2D(Form("hW%s_%sVs%sScaled_8_JetCut_pdf%d-%d%s_eta%s_%d",MuCharge_str.Data(),WMass::FitVar_str[2].Data(),WMass::FitVar_str[0].Data(),WMass::PDF_sets<0?generated_PDF_set:WMass::PDF_sets,h,toys_str.Data(),eta_str.Data(),jWmass),
                                  Mu_var_NotScaled[0], Mu_var_NotScaled[2],evt_weight*Mu_tight_muon_SF*lha_weight, h_2d, 
                                  // nbins, bins_Notscaled[0], nbins, bins_Notscaled[2]  );
                                  50, WMass::fit_xmin[0], WMass::fit_xmax[0], 50, WMass::fit_xmin[2], WMass::fit_xmax[2] );
                          common_stuff::plot2D(Form("hW%s_%sVs%sScaled_8_JetCut_pdf%d-%d%s_eta%s_%d",MuCharge_str.Data(),WMass::FitVar_str[2].Data(),WMass::FitVar_str[1].Data(),WMass::PDF_sets<0?generated_PDF_set:WMass::PDF_sets,h,toys_str.Data(),eta_str.Data(),jWmass),
                                  Mu_var_NotScaled[1],Mu_var_NotScaled[2],evt_weight*Mu_tight_muon_SF*lha_weight, h_2d, 
                                  // nbins, bins_Notscaled[1], nbins, bins_Notscaled[2]  );
                                  50, WMass::fit_xmin[1], WMass::fit_xmax[1], 50, WMass::fit_xmin[2], WMass::fit_xmax[2] );
                        }
                        
                        // hWPos_VarNonScaled_8_JetCut[i][j]->Fill(mu.Pt()*iWmass/(WMass::WMassCentral_MeV/1e3)<xmax*80/2 ? mu.Pt() : (xmax-binsize2/2)*80/2 ,evt_weight*Mu_tight_muon_SF);

                        // cout << wmass1 << " " << WMass::WMassCentral_MeV << " " << (wmass1 - WMass::WMassCentral_MeV) << endl;
                        // cout << WMass::etaMaxMuons[i]  << " " << 2.1 << " " << ((WMass::etaMaxMuons[i] - 2.1)) << endl;
                        // if( (TMath::Abs(wmass1 - WMass::WMassCentral_MeV) > 1)
                          // // || (TMath::Abs(WMass::etaMaxMuons[i] - 2.1) > 1e-3 ) 
                          // ) 
                          // continue;
                          
                        // cout << "wmass1 " << wmass1 << " WMass::etaMaxMuons[i]= " << WMass::etaMaxMuons[i] << endl;

                        if(controlplots){
                          // if(useGenVar) hWPos_MuDRgen[0][i][j]->Fill(TMath::Log10(MuDRGenP),evt_weight*Mu_tight_muon_SF);
                          if(useGenVar) common_stuff::plot1D(Form("hW%s_MuDRgen_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[0].Data(),eta_str.Data(),jWmass),
                                    TMath::Log10(MuDRGenP),evt_weight*Mu_tight_muon_SF,
                                    h_1d, 1000,-6,1 );
                          // hnvtx[0][i][j]->Fill(nvtx,evt_weight); // TO FIT MET IN THE WHOLE RANGE!!!!
                          common_stuff::plot1D(Form("hnvtx_%s_eta%s_%d",WMass::nSigOrQCD_str[0].Data(),eta_str.Data(),jWmass),
                                    nvtx, evt_weight,
                                    h_1d, 50, 0, 50 );
                          // hnoTrgMuonsLeadingPt[0][i][j]->Fill(noTrgMuonsLeadingPt,evt_weight*Mu_tight_muon_SF); // TO FIT MET IN THE WHOLE RANGE!!!!
                          common_stuff::plot1D(Form("hnoTrgMuonsLeadingPt_%s_eta%s_%d",WMass::nSigOrQCD_str[0].Data(),eta_str.Data(),jWmass),
                                    noTrgMuonsLeadingPt, evt_weight*Mu_tight_muon_SF,
                                    h_1d, 300, -100, 200 );
                          // hpfMETphi_WPos[0][i][j]->Fill(pfmet_phi,evt_weight*Mu_tight_muon_SF);
                          common_stuff::plot1D(Form("hpfMETphi_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[0].Data(),eta_str.Data(),jWmass),
                                    pfmet_phi,evt_weight*Mu_tight_muon_SF,
                                    h_1d, 100, -TMath::Pi(), TMath::Pi() );
                          // hWPos_pt[0][i][j]->Fill(W_pt,evt_weight*Mu_tight_muon_SF);
                          common_stuff::plot1D(Form("hW%s_pt_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[0].Data(),eta_str.Data(),jWmass),
                                    W_pt, evt_weight*Mu_tight_muon_SF,
                                    h_1d, 100, 0, 25 );
                          // hWPos_phi[0][i][j]->Fill(W_phi,evt_weight*Mu_tight_muon_SF);
                          common_stuff::plot1D(Form("hW%s_phi_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[0].Data(),eta_str.Data(),jWmass),
                                    W_phi, evt_weight*Mu_tight_muon_SF,
                                    h_1d, 100,-TMath::Pi(),TMath::Pi() );
                          // hWPos_mt[0][i][j]->Fill(W_mt,evt_weight*Mu_tight_muon_SF);
                          common_stuff::plot1D(Form("hW%s_mt_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[0].Data(),eta_str.Data(),jWmass),
                                    W_mt, evt_weight*Mu_tight_muon_SF,
                                    h_1d, 200,0,200 );
                          // hMupt_WPos[0][i][j]->Fill(mu.Pt(),evt_weight*Mu_tight_muon_SF);
                          common_stuff::plot1D(Form("hMupt_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[0].Data(),eta_str.Data(),jWmass),
                                    mu.Pt(),evt_weight*Mu_tight_muon_SF,
                                    h_1d, 200,0,200 );
                          // hMueta_WPos[0][i][j]->Fill(mu.Eta(),evt_weight*Mu_tight_muon_SF);
                          common_stuff::plot1D(Form("hMueta_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[0].Data(),eta_str.Data(),jWmass),
                                    mu.Eta(),evt_weight*Mu_tight_muon_SF,
                                    h_1d, 100,-2.5,2.5 );
                          // hMuphi_WPos[0][i][j]->Fill(mu.Phi(),evt_weight*Mu_tight_muon_SF);
                          common_stuff::plot1D(Form("hMuphi_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[0].Data(),eta_str.Data(),jWmass),
                                    mu.Phi(),evt_weight*Mu_tight_muon_SF,
                                    h_1d, 100,-TMath::Pi(),TMath::Pi() );
                          // hMulogiso_WPos[0][i][j]->Fill(TMath::Log10(MuRelIso),evt_weight*Mu_tight_muon_SF);
                          common_stuff::plot1D(Form("hMulogiso_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[0].Data(),eta_str.Data(),jWmass),
                                    TMath::Log10(MuRelIso),evt_weight*Mu_tight_muon_SF,
                                    h_1d, 1000,-5,3 );
                          // hJetpt_WPos[0][i][j]->Fill(Jet_leading_pt,evt_weight*Mu_tight_muon_SF);
                          common_stuff::plot1D(Form("hJetpt_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[0].Data(),eta_str.Data(),jWmass),
                                    Jet_leading_pt,evt_weight*Mu_tight_muon_SF,
                                    h_1d, 100,0,50 );
                          // hJeteta_WPos[0][i][j]->Fill(Jet_leading_eta,evt_weight*Mu_tight_muon_SF);
                          common_stuff::plot1D(Form("hJeteta_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[0].Data(),eta_str.Data(),jWmass),
                                    Jet_leading_eta,evt_weight*Mu_tight_muon_SF,
                                    h_1d, 100,-2.5,2.5 );
                          // hJetphi_WPos[0][i][j]->Fill(Jet_leading_phi,evt_weight*Mu_tight_muon_SF);
                          common_stuff::plot1D(Form("hJetphi_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[0].Data(),eta_str.Data(),jWmass),
                                    Jet_leading_phi,evt_weight*Mu_tight_muon_SF,
                                    h_1d, 100,-TMath::Pi(),TMath::Pi() );
                        }
                      }
                    }
                  }
                  if(controlplots && W_pt<20 && m==0){ // do not cut on MET to let the fit have an handle
                          
                    // if(TMath::Abs(wmass1 - WMass::WMassCentral_MeV) > 1) continue;
                      
                    // hpfMET_WPos[0][i][j]->Fill(pfmet,evt_weight*Mu_tight_muon_SF); // TRICK TO FIT MET IN THE WHOLE RANGE!!!!
                    common_stuff::plot1D(Form("hpfMET_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[0].Data(),eta_str.Data(),jWmass),
                                    pfmet,evt_weight*Mu_tight_muon_SF, 
                                    h_1d, 200,0,200 );
                      
                  }
                }else if(controlplots && m==0){ // muon candidate is failing either tight ID, iso or dxy: QCD enriched region
                  
                  if(pfmet>25 && W_pt<20){
                    // if( (TMath::Abs(wmass1 - WMass::WMassCentral_MeV) > 1)) continue;
                    // hWPos_logiso_vs_logdxy[i][j]->Fill(TMath::Log10(Mu_dxy),TMath::Log10(MuRelIso),evt_weight*Mu_tight_muon_SF);
                    // hWPos_iso_vs_dxy[i][j]->Fill(Mu_dxy,MuRelIso,evt_weight*Mu_tight_muon_SF);
                  }
                    
                  if(
                    MuRelIso>0.12 // single muon cuts (inverted iso (is <0.5 for signal) , no tight requirement)
                    && Mu_dxy>0.02 // single muon cuts (inverted iso (is <0.5 for signal) , no tight requirement)
                    && W_pt<20
                  ){
                  
                    // hpfMET_WPos[1][i][j]->Fill(pfmet,evt_weight*Mu_tight_muon_SF);  // TO FIT MET IN THE WHOLE RANGE!!!!  
                    common_stuff::plot1D(Form("hpfMET_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[1].Data(),eta_str.Data(),jWmass),
                                    pfmet,evt_weight*Mu_tight_muon_SF, 
                                    h_1d, 200,0,200 );
                    
                    if(pfmet>25){
                      for(int k=0;k<3;k++)
                        // hWPos_VarScaled_QCD[k][i][j]->Fill(Mu_var_jacobian[k],evt_weight*Mu_tight_muon_SF);
                        common_stuff::plot1D(Form("hW%s_%sScaled_QCD_eta%s_%d",MuCharge_str.Data(),WMass::FitVar_str[k].Data(),eta_str.Data(),jWmass),
                                    Mu_var_jacobian[k],evt_weight*Mu_tight_muon_SF, h_1d, 
                                    // nbins, bins_scaled[k] );
                                    50, WMass::fit_xmin[k]/(WMass::WMassCentral_MeV/1e3),WMass::fit_xmax[k]/(WMass::WMassCentral_MeV/1e3) );
                      
                      // if( (TMath::Abs(wmass1 - WMass::WMassCentral_MeV) > 1)
                        // )
                        // continue;
                      
                      // if(useGenVar) hWPos_MuDRgen[1][i][j]->Fill(TMath::Log10(MuDRGenP),evt_weight*Mu_tight_muon_SF);
                      common_stuff::plot1D(Form("hW%s_MuDRgen_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[1].Data(),eta_str.Data(),jWmass),
                                    TMath::Log10(MuDRGenP),evt_weight*Mu_tight_muon_SF,
                                    h_1d, 1000,-6,1 );
                      // hnvtx[1][i][j]->Fill(nvtx,evt_weight);  // TO FIT MET IN THE WHOLE RANGE!!!!  
                      common_stuff::plot1D(Form("hnvtx_%s_eta%s_%d",WMass::nSigOrQCD_str[1].Data(),eta_str.Data(),jWmass),
                                    nvtx, evt_weight,
                                    h_1d, 50, 0, 50 );
                      // hnoTrgMuonsLeadingPt[1][i][j]->Fill(noTrgMuonsLeadingPt,evt_weight);  // TO FIT MET IN THE WHOLE RANGE!!!!  
                      common_stuff::plot1D(Form("hnoTrgMuonsLeadingPt_%s_eta%s_%d",WMass::nSigOrQCD_str[1].Data(),eta_str.Data(),jWmass),
                                    noTrgMuonsLeadingPt, evt_weight*Mu_tight_muon_SF,
                                    h_1d, 300, -100, 200 );
                      // hpfMETphi_WPos[1][i][j]->Fill(pfmet_phi,evt_weight*Mu_tight_muon_SF);
                      common_stuff::plot1D(Form("hpfMETphi_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[1].Data(),eta_str.Data(),jWmass),
                                    pfmet_phi,evt_weight*Mu_tight_muon_SF,
                                    h_1d, 100, -TMath::Pi(), TMath::Pi() );
                      // hWPos_pt[1][i][j]->Fill(W_pt,evt_weight*Mu_tight_muon_SF);
                      common_stuff::plot1D(Form("hW%s_pt_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[1].Data(),eta_str.Data(),jWmass),
                                    W_pt, evt_weight*Mu_tight_muon_SF,
                                    h_1d, 100, 0, 25 );
                      // hWPos_phi[1][i][j]->Fill(W_phi,evt_weight*Mu_tight_muon_SF);
                      common_stuff::plot1D(Form("hW%s_phi_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[1].Data(),eta_str.Data(),jWmass),
                                    W_phi, evt_weight*Mu_tight_muon_SF,
                                    h_1d, 100,-TMath::Pi(),TMath::Pi() );
                      // hWPos_mt[1][i][j]->Fill(W_mt,evt_weight*Mu_tight_muon_SF);
                      common_stuff::plot1D(Form("hW%s_mt_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[1].Data(),eta_str.Data(),jWmass),
                                    W_mt, evt_weight*Mu_tight_muon_SF,
                                    h_1d, 200,0,200 );
                      // hMupt_WPos[1][i][j]->Fill(mu.Pt(),evt_weight*Mu_tight_muon_SF);
                      common_stuff::plot1D(Form("hMupt_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[1].Data(),eta_str.Data(),jWmass),
                                    mu.Pt(),evt_weight*Mu_tight_muon_SF,
                                    h_1d, 200,0,200 );
                      // hMueta_WPos[1][i][j]->Fill(mu.Eta(),evt_weight*Mu_tight_muon_SF);
                      common_stuff::plot1D(Form("hMueta_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[1].Data(),eta_str.Data(),jWmass),
                                    mu.Eta(),evt_weight*Mu_tight_muon_SF,
                                    h_1d, 100,-2.5,2.5 );
                      // hMuphi_WPos[1][i][j]->Fill(mu.Phi(),evt_weight*Mu_tight_muon_SF);
                      common_stuff::plot1D(Form("hMuphi_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[1].Data(),eta_str.Data(),jWmass),
                                    mu.Phi(),evt_weight*Mu_tight_muon_SF,
                                    h_1d, 100,-TMath::Pi(),TMath::Pi() );
                      // hMulogiso_WPos[1][i][j]->Fill(TMath::Log10(MuRelIso),evt_weight*Mu_tight_muon_SF);
                      common_stuff::plot1D(Form("hMulogiso_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[1].Data(),eta_str.Data(),jWmass),
                                    TMath::Log10(MuRelIso),evt_weight*Mu_tight_muon_SF,
                                    h_1d, 1000,-5,3 );
                      // hJetpt_WPos[1][i][j]->Fill(Jet_leading_pt,evt_weight*Mu_tight_muon_SF);
                      common_stuff::plot1D(Form("hJetpt_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[1].Data(),eta_str.Data(),jWmass),
                                    Jet_leading_pt,evt_weight*Mu_tight_muon_SF,
                                    h_1d, 100,0,50 );
                      // hJeteta_WPos[1][i][j]->Fill(Jet_leading_eta,evt_weight*Mu_tight_muon_SF);
                      common_stuff::plot1D(Form("hJeteta_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[1].Data(),eta_str.Data(),jWmass),
                                    Jet_leading_eta,evt_weight*Mu_tight_muon_SF,
                                    h_1d, 100,-2.5,2.5 );
                      // hJetphi_WPos[1][i][j]->Fill(Jet_leading_phi,evt_weight*Mu_tight_muon_SF);
                      common_stuff::plot1D(Form("hJetphi_W%s_%s_eta%s_%d",MuCharge_str.Data(),WMass::nSigOrQCD_str[1].Data(),eta_str.Data(),jWmass),
                                    Jet_leading_phi,evt_weight*Mu_tight_muon_SF,
                                    h_1d, 100,-TMath::Pi(),TMath::Pi() );
                    }
                  }
                }
              } // end if for good muon in acceptance
            } // end good reco event
          } // end momentum correction toys
        } // end dummy separation sig/bkg
      } // end W mass loop
    } // end muon eta loop
  } // end event loop
  // }
  
  outTXTfile.close();
  
  TFile*fout = new TFile(Form("%s/Wanalysis%s.root",outputdir.Data(),chunk_str.Data()),"RECREATE");
  
    if(!sampleName.Contains("WJetsSig")){
      for(int i=0; i<WMass::etaMuonNSteps; i++){
        TString eta_str = Form("%.1f",WMass::etaMaxMuons[i]); eta_str.ReplaceAll(".","p");
        for(int j=0; j<2*WMass::WMassNSteps+1; j++){
          // if(WMass::WMassNSteps!=j){
            int jWmass = (WMass::WMassCentral_MeV-(WMass::WMassNSteps-j)*WMass::WMassStep_MeV);
            for(int k=0;k<3;k++){
              for(int h=0; h<WMass::PDF_members; h++){
                for(int m=0; m<WMass::NtoysMomCorr; m++){
                  TString toys_str = "";
                  if(WMass::NtoysMomCorr>1) toys_str = Form("_MomCorrToy%d",m);
                  TString MuCharge_str[]={"Pos","Neg"};
                  for(int MuCharge = 0; MuCharge < 2; MuCharge++){
                    // cout << "creating " << Form("hW%s_%sScaled_8_JetCut_pdf%d-%d%s_eta%s_%d",MuCharge_str[MuCharge].Data(),WMass::FitVar_str[k].Data(),WMass::PDF_sets<0?generated_PDF_set:WMass::PDF_sets,h,toys_str.Data(),eta_str.Data(),jWmass) << endl;
                    common_stuff::cloneHisto1D(Form("hW%s_%sScaled_8_JetCut_pdf%d-%d%s_eta%s_%d",MuCharge_str[MuCharge].Data(),WMass::FitVar_str[k].Data(),WMass::PDF_sets<0?generated_PDF_set:WMass::PDF_sets,h,toys_str.Data(),eta_str.Data(),WMass::WMassCentral_MeV), 
                                              Form("hW%s_%sScaled_8_JetCut_pdf%d-%d%s_eta%s_%d",MuCharge_str[MuCharge].Data(),WMass::FitVar_str[k].Data(),WMass::PDF_sets<0?generated_PDF_set:WMass::PDF_sets,h,toys_str.Data(),eta_str.Data(),jWmass), 
                                              h_1d);
                    // hWPos_VarScaled_8_JetCut[m][h][k][i][j]=(TH1D*)hWPos_VarScaled_8_JetCut[m][h][k][i][WMass::WMassNSteps]->Clone(Form("hWPos_%sScaled_8_JetCut_pdf%d-%d%s_eta%s_%d",WMass::FitVar_str[k].Data(),WMass::PDF_sets<0?generated_PDF_set:WMass::PDF_sets,h,toys_str.Data(),eta_str.Data(),jWmass));
                    // hWPos_VarScaled_8_JetCut[m][h][k][i][j]->SetName(Form("hWPos_%sScaled_8_JetCut_pdf%d-%d%s_eta%s_%d",WMass::FitVar_str[k].Data(),WMass::PDF_sets<0?generated_PDF_set:WMass::PDF_sets,h,toys_str.Data(),eta_str.Data(),jWmass));
                    // hWPos_VarScaled_8_JetCut[m][h][k][i][j]->SetTitle(Form("hWPos_%sScaled_8_JetCut_pdf%d-%d%s_eta%s_%d",WMass::FitVar_str[k].Data(),WMass::PDF_sets<0?generated_PDF_set:WMass::PDF_sets,h,toys_str.Data(),eta_str.Data(),jWmass));
                    // hWPos_VarScaled_8_JetCut[m][h][k][i][j]->Write();
                    
                    common_stuff::cloneHisto1D(Form("hW%s_%sNonScaled_8_JetCut_pdf%d-%d%s_eta%s_%d",MuCharge_str[MuCharge].Data(),WMass::FitVar_str[k].Data(),WMass::PDF_sets<0?generated_PDF_set:WMass::PDF_sets,h,toys_str.Data(),eta_str.Data(),WMass::WMassCentral_MeV), 
                                              Form("hW%s_%sNonScaled_8_JetCut_pdf%d-%d%s_eta%s_%d",MuCharge_str[MuCharge].Data(),WMass::FitVar_str[k].Data(),WMass::PDF_sets<0?generated_PDF_set:WMass::PDF_sets,h,toys_str.Data(),eta_str.Data(),jWmass), 
                                              h_1d);
                    // hWPos_VarNonScaled_8_JetCut[m][h][k][i][j]=(TH1D*)hWPos_VarNonScaled_8_JetCut[m][h][k][i][WMass::WMassNSteps]->Clone(Form("hWPos_%sNonScaled_8_JetCut_pdf%d-%d%s_eta%s_%d",WMass::FitVar_str[k].Data(),WMass::PDF_sets<0?generated_PDF_set:WMass::PDF_sets,h,toys_str.Data(),eta_str.Data(),jWmass));
                    // hWPos_VarNonScaled_8_JetCut[m][h][k][i][j]->SetName(Form("hWPos_%sNonScaled_8_JetCut_pdf%d-%d%s_eta%s_%d",WMass::FitVar_str[k].Data(),WMass::PDF_sets<0?generated_PDF_set:WMass::PDF_sets,h,toys_str.Data(),eta_str.Data(),jWmass));
                    // hWPos_VarNonScaled_8_JetCut[m][h][k][i][j]->SetTitle(Form("hWPos_%sNonScaled_8_JetCut_pdf%d-%d%s_eta%s_%d",WMass::FitVar_str[k].Data(),WMass::PDF_sets<0?generated_PDF_set:WMass::PDF_sets,h,toys_str.Data(),eta_str.Data(),jWmass));
                    // hWPos_VarNonScaled_8_JetCut[m][h][k][i][j]->Write();
                  }
                }
              }
            }
          // }
        }
      }
    }

    common_stuff::writeOutHistos( fout, h_1d, h_2d );
  
  fout->cd();
  
  fout->Write();
  fout->Close();
  
}
