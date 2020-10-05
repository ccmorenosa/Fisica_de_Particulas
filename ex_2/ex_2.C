#include <iostream>
#include <string>
#include <stdio.h>

#pragma link C++ enum  myns::myenum;

void ex_2 () {
  // Open de Root File
  TFile *file = TFile::Open("../ROOT_Files/Tracks_Clusters.root");

  // Get the Tree from the Root File
  TTree *tree = (TTree*) file->Get("JetRecoTree");
  tree->Print();

  // Pointers and variables to stores the branches data
  UInt_t npv = -1;
  Float_t mu_avg = -1;
  float evtw = -1;
  vector<float> *reco_R4_pt = {};
  vector<float> *truth_R4_pt = {};
  vector<float> *track_R4_pt = {};
  vector<float> *reco_R4_jvf = {};

  vector<float> *reco_R4_eta = {};
  vector<float> *reco_R4_phi = {};
  vector<float> *reco_R4_m = {};

  vector<float> *truth_R4_eta = {};
  vector<float> *truth_R4_phi = {};
  vector<float> *truth_R4_m = {};

  vector<float> *track_R4_eta = {};
  vector<float> *track_R4_phi = {};
  vector<float> *track_R4_m = {};

  vector<float> *reco_R10_pt = {};
  vector<float> *reco_R10_trimmed_pt = {};
  vector<float> *reco_R10_m = {};
  vector<float> *reco_R10_trimmed_m = {};

  vector<float> *truth_R10_pt = {};
  vector<float> *truth_R10_trimmed_pt = {};
  vector<float> *truth_R10_m = {};
  vector<float> *truth_R10_trimmed_m = {};


  // Get the branches address into the variables
  tree->SetBranchAddress("NPV", &npv);

  tree->SetBranchAddress("mu_average", &mu_avg);

  tree->SetBranchAddress("EventWeight", &evtw);
  tree->SetBranchAddress("RecoJets_R4_pt", &reco_R4_pt);
  tree->SetBranchAddress("TruthJets_R4_pt", &truth_R4_pt);
  tree->SetBranchAddress("TrackJets_R4_pt", &track_R4_pt);

  tree->SetBranchAddress("RecoJets_R4_jvf", &reco_R4_jvf);

  tree->SetBranchAddress("RecoJets_R4_eta", &reco_R4_eta);
  tree->SetBranchAddress("RecoJets_R4_phi", &reco_R4_phi);
  tree->SetBranchAddress("RecoJets_R4_m", &reco_R4_m);

  tree->SetBranchAddress("TruthJets_R4_eta", &truth_R4_eta);
  tree->SetBranchAddress("TruthJets_R4_phi", &truth_R4_phi);
  tree->SetBranchAddress("TruthJets_R4_m", &truth_R4_m);

  tree->SetBranchAddress("TrackJets_R4_eta", &track_R4_eta);
  tree->SetBranchAddress("TrackJets_R4_phi", &track_R4_phi);
  tree->SetBranchAddress("TrackJets_R4_m", &track_R4_m);

  tree->SetBranchAddress("RecoJets_R10_pt", &reco_R10_pt);
  tree->SetBranchAddress("RecoJets_R10_Trimmed_pt", &reco_R10_trimmed_pt);
  tree->SetBranchAddress("RecoJets_R10_m", &reco_R10_m);
  tree->SetBranchAddress("RecoJets_R10_Trimmed_m", &reco_R10_trimmed_m);

  tree->SetBranchAddress("TruthJets_R10_pt", &truth_R10_pt);
  tree->SetBranchAddress("TruthJets_R10_Trimmed_pt", &truth_R10_trimmed_pt);
  tree->SetBranchAddress("TruthJets_R10_m", &truth_R10_m);
  tree->SetBranchAddress("TruthJets_R10_Trimmed_m", &truth_R10_trimmed_m);

  // Create the canvas to plot the histograms
  TCanvas *canvas = new TCanvas(
    "Canvas",
    "a first way to plot a variable",
    800,
    600
  );


  //--------------------------------------------------------------------------//
  //                               Leading Jets                               //
  //--------------------------------------------------------------------------//

  // With weigth
  TH1F *hist_reco_pt_w = new TH1F(
    "Reco-jet",
    "jet pT [Weigth]; pT(GeV);Events",
    50,
    10,
    200
  );

  TH1F *hist_truth_pt_w = new TH1F(
    "Truth-jet",
    "jet pT [Weigth]; pT(GeV);Events",
    50,
    10,
    200
  );

  TH1F *hist_leadreco_pt_w = new TH1F(
    "Lead Reco-jet [Weight]",
    "Leading jet pT [Weigth]; pT(GeV);Events",
    50,
    10,
    200
  );

  TH1F *hist_leadtruth_pt_w = new TH1F(
    "Lead Truth-jet [Weight]",
    "Leading jet pT [Weigth]; pT(GeV);Events",
    50,
    10,
    200
  );


  // Without weigth
  TH1F *hist_leadreco_pt = new TH1F(
    "Lead Reco-jet",
    "Leading jet pT; pT(GeV);Events",
    50,
    10,
    200
  );

  TH1F *hist_leadtruth_pt = new TH1F(
    "Lead Truth-jet",
    "Leading jet pT; pT(GeV);Events",
    50,
    10,
    200
  );


  int nentries, nbytes, i;
  nentries = (Int_t)tree->GetEntries();

  for (i = 0; i < nentries; i++) {
    nbytes = tree->GetEntry(i);

    if(reco_R4_pt->size()>0){
      hist_leadreco_pt_w->Fill(reco_R4_pt->at(0)/1000.,evtw);
      hist_leadtruth_pt_w->Fill(truth_R4_pt->at(0)/1000.,evtw);
      hist_leadreco_pt->Fill(reco_R4_pt->at(0)/1000.);
      hist_leadtruth_pt->Fill(truth_R4_pt->at(0)/1000.);

      for(int j=0; j<reco_R4_pt->size(); j++){
        hist_reco_pt_w->Fill(reco_R4_pt->at(j)/1000.,evtw);
      }

      for(int j=0; j<truth_R4_pt->size(); j++){
        hist_truth_pt_w->Fill(truth_R4_pt->at(j)/1000.,evtw);
      }
    }
  }

  std::cout << "Done Leading Raco/Truth" << std::endl;

  hist_reco_pt_w->SetMarkerStyle(20);
  hist_reco_pt_w->SetMarkerColor(kRed);
  hist_truth_pt_w->SetMarkerStyle(21);
  hist_truth_pt_w->SetMarkerColor(kBlue);

  hist_leadreco_pt_w->SetMarkerStyle(20);
  hist_leadreco_pt_w->SetMarkerColor(kRed);
  hist_leadtruth_pt_w->SetMarkerStyle(21);
  hist_leadtruth_pt_w->SetMarkerColor(kBlue);

  hist_leadreco_pt->SetFillColorAlpha(kRed, 0.5);
  hist_leadtruth_pt->SetFillColorAlpha(kBlue, 0.5);

  hist_reco_pt_w->Draw();
  hist_truth_pt_w->Draw("same");
  canvas->SetLogy();
  canvas->Print("truth_reco_w.svg");
  canvas->Clear();

  hist_leadreco_pt_w->Draw();
  hist_leadtruth_pt_w->Draw("same");
  canvas->SetLogy();
  canvas->Print("Lead_truth_reco_w.svg");
  canvas->Clear();

  hist_leadreco_pt->Draw();
  hist_leadtruth_pt->Draw("same");
  canvas->SetLogy();
  canvas->Print("Lead_truth_reco.svg");
  canvas->Clear();


  //--------------------------------------------------------------------------//
  //                              Jets vs Pileup                              //
  //--------------------------------------------------------------------------//

  // Histograms
  TH2F *hist_r_jetpt_npv = new TH2F(
    "Reco-jet pT vs. NPV",
    ";NPV; jet pT",
    50,
    10,
    100,
    20,
    0,
    60
  );

  TH2F *hist_r_jetpt_mu = new TH2F(
    "Reco-jet pT vs. mu average",
    ";mu_avg; jet pT",
    50,
    10,
    100,
    20,
    0,
    60
  );

  TH2F *hist_t_jetpt_npv = new TH2F(
    "Truth-jet pT vs. NPV",
    ";NPV; jet pT",
    50,
    10,
    100,
    20,
    0,
    60
  );

  TH2F *hist_t_jetpt_mu = new TH2F(
    "Truth-jet pT vs. mu average",
    ";mu_avg; jet pT",
    50,
    10,
    100,
    20,
    0,
    60
  );


  // Profiles
  TProfile *prof_r_jetpt_npv = new TProfile(
    "Profile Reco-jet pT vs. NPV",
    ";NPV; jet pT",
    50,
    1,
    50,
    20,
    200
  );

  TProfile *prof_r_jetpt_mu = new TProfile(
    "Profile Reco-jet pT vs. mu average",
    ";mu_avg; jet pT",
    50,
    1,
    50,
    30,
    200
  );

  TProfile *prof_t_jetpt_npv = new TProfile(
    "Profile Truth-jet pT vs. NPV",
    ";NPV; jet pT",
    50,
    1,
    50,
    18,
    200
  );

  TProfile *prof_t_jetpt_mu = new TProfile(
    "Profile Truth-jet pT vs. mu average",
    ";mu_avg; jet pT",
    50,
    1,
    50,
    30,
    200
  );

  for (i = 0; i < nentries; i++) {
    nbytes = tree->GetEntry(i);

    if(reco_R4_pt->size()!=0 && reco_R4_pt->at(0)>20000.){
      for(int j=0; j<reco_R4_pt->size(); j++){
        hist_r_jetpt_npv->Fill(reco_R4_pt->at(j)/1000.,npv,evtw);
        prof_r_jetpt_npv->Fill(reco_R4_pt->at(j)/1000.,npv,evtw);
        hist_r_jetpt_mu->Fill(reco_R4_pt->at(j)/1000.,mu_avg,evtw);
        prof_r_jetpt_mu->Fill(reco_R4_pt->at(j)/1000.,mu_avg,evtw);
      }
    }

    if(truth_R4_pt->size()!=0 && truth_R4_pt->at(0)>20000.){
      for(int j=0; j<truth_R4_pt->size(); j++){
        hist_t_jetpt_npv->Fill(truth_R4_pt->at(j)/1000.,npv,evtw);
        prof_t_jetpt_npv->Fill(truth_R4_pt->at(j)/1000.,npv,evtw);
        hist_t_jetpt_mu->Fill(truth_R4_pt->at(j)/1000.,mu_avg,evtw);
        prof_t_jetpt_mu->Fill(truth_R4_pt->at(j)/1000.,mu_avg,evtw);
      }
    }
  }

  std::cout << "Done Jets vs. Pileup" << std::endl;

  // Draw Histograms
  hist_r_jetpt_npv->Draw("colz");
  canvas->SetLogy(0);
  canvas->Print("reco_npv.svg");
  canvas->Clear();

  hist_r_jetpt_mu->Draw("colz");
  canvas->Print("reco_mu.svg");
  canvas->Clear();

  hist_t_jetpt_npv->Draw("colz");
  canvas->Print("truth_npv.svg");
  canvas->Clear();

  hist_t_jetpt_mu->Draw("colz");
  canvas->Print("truth_mu.svg");
  canvas->Clear();


  // Draw Profiles
  prof_r_jetpt_npv->Draw("colz");
  canvas->Print("prof_reco_npv.svg");
  canvas->Clear();

  prof_r_jetpt_mu->Draw("colz");
  canvas->Print("prof_reco_mu.svg");
  canvas->Clear();

  prof_t_jetpt_npv->Draw("colz");
  canvas->Print("prof_truth_npv.svg");
  canvas->Clear();

  prof_t_jetpt_mu->Draw("colz");
  canvas->Print("prof_truth_mu.svg");
  canvas->Clear();

  //--------------------------------------------------------------------------//
  //                                 JVF Cuts                                 //
  //--------------------------------------------------------------------------//

  // With cuts
  TH1F *hist_leadreco_pt_compare_c = new TH1F(
    "Lead Reco-jet [Cuts]",
    "Leading jet pT; pT (GeV);Events",
    20,
    0,
    200
  );

  TH1F *hist_leadtruth_pt_compare_c = new TH1F(
    "Lead Truth-jet [Cuts]",
    "Leading jet pT; pT (GeV);Events",
    20,
    0,
    200
  );

  TH1F *hist_leadtrack_pt_compare_c = new TH1F(
    "Lead Track-jet [Cuts]",
    "Leading jet pT; pT (GeV);Events",
    20,
    0,
    200
  );


  // Without cuts
  TH1F *hist_leadreco_pt_compare = new TH1F(
    "Lead Reco-jet Compare",
    "Leading jet pT; pT (GeV);Events",
    20,
    0,
    200
  );

  TH1F *hist_leadtruth_pt_compare = new TH1F(
    "Lead Truth-jet Compare",
    "Leading jet pT; pT (GeV);Events",
    20,
    0,
    200
  );

  TH1F *hist_leadtrack_pt_compare = new TH1F(
    "Lead Track-jet Compare",
    "Leading jet pT; pT (GeV);Events",
    20,
    0,
    200
  );

  for (i = 0; i < nentries; i++) {
    nbytes = tree->GetEntry(i);
    if(reco_R4_pt->size()!=0 &&
    reco_R4_pt->at(0)>20000.){
      hist_leadreco_pt_compare->Fill(reco_R4_pt->at(0)/1000., evtw);
      hist_leadtruth_pt_compare->Fill(truth_R4_pt->at(0)/1000., evtw);

      if(std::abs(reco_R4_jvf->at(0)) < 0.5){
        hist_leadreco_pt_compare_c->Fill(reco_R4_pt->at(0)/1000., evtw);
        hist_leadtruth_pt_compare_c->Fill(truth_R4_pt->at(0)/1000., evtw);
      }
    }

    if(track_R4_pt->size()!=0){
      hist_leadtrack_pt_compare->Fill(track_R4_pt->at(0)/1000., evtw);

      if(reco_R4_pt->size()!=0 &&
      std::abs(reco_R4_jvf->at(0)) < 0.5){
        hist_leadtrack_pt_compare_c->Fill(track_R4_pt->at(0)/1000., evtw);
      }
    }

  }

  std::cout << "Done Cuts" << std::endl;
  canvas->SetLogy();

  // Compare Truth and Reco

  // Draw histograms with cuts
  hist_leadreco_pt_compare_c->SetMarkerStyle(20);
  hist_leadreco_pt_compare_c->SetMarkerColor(kRed);
  hist_leadreco_pt_compare_c->Draw("");

  hist_leadtruth_pt_compare_c->SetMarkerStyle(20);
  hist_leadtruth_pt_compare_c->SetMarkerColor(kBlue);
  hist_leadtruth_pt_compare_c->Draw("Same");

  // Draw histograms without cuts
  hist_leadreco_pt_compare->SetMarkerStyle(20);
  hist_leadreco_pt_compare->SetMarkerColor(kOrange);
  hist_leadreco_pt_compare->Draw("Same");

  hist_leadtruth_pt_compare->SetMarkerStyle(20);
  hist_leadtruth_pt_compare->SetMarkerColor(kCyan);
  hist_leadtruth_pt_compare->Draw("Same");

  canvas->Print("lead_truth_reco_cut.svg");
  canvas->Clear();

  // Compare Track and Reco

  // Draw histograms with cuts
  hist_leadreco_pt_compare_c->SetMarkerStyle(20);
  hist_leadreco_pt_compare_c->SetMarkerColor(kRed);
  hist_leadreco_pt_compare_c->Draw("");

  hist_leadtrack_pt_compare_c->SetMarkerStyle(20);
  hist_leadtrack_pt_compare_c->SetMarkerColor(kBlue);
  hist_leadtrack_pt_compare_c->Draw("Same");

  // Draw histograms without cuts
  hist_leadreco_pt_compare->SetMarkerStyle(20);
  hist_leadreco_pt_compare->SetMarkerColor(kOrange);
  hist_leadreco_pt_compare->Draw("Same");

  hist_leadtrack_pt_compare->SetMarkerStyle(20);
  hist_leadtrack_pt_compare->SetMarkerColor(kCyan);
  hist_leadtrack_pt_compare->Draw("Same");

  canvas->Print("lead_track_reco_cut.svg");
  canvas->Clear();

  //--------------------------------------------------------------------------//
  //                                  DeltaR                                  //
  //--------------------------------------------------------------------------//

  TH1F *hist_DR_reco_truth = new TH1F(
    "Delta R reco",
    "Delta R; #Delta R; Events",
    20,
    0,
    2
  );

  TH1F *hist_DR_reco_cut_truth = new TH1F(
    "Delta R reco [Cuts]",
    "Delta R; #Delta R; Events",
    20,
    0,
    2
  );

  TH1F *hist_DR_track_truth = new TH1F(
    "Delta R tracks",
    "Delta R; #Delta R; Events",
    20,
    0,
    2
  );

  for (i = 0; i < nentries; i++) {
    nbytes = tree->GetEntry(i);

    if(truth_R4_pt->size()!=0 && truth_R4_pt->at(0)>20000.){
      TLorentzVector truthJet;
      truthJet.SetPtEtaPhiM(
        truth_R4_pt->at(0),
        truth_R4_eta->at(0),
        truth_R4_phi->at(0),
        truth_R4_m->at(0)
      );

      if(reco_R4_pt->size()!=0 && fabs(reco_R4_jvf->at(0))>0.5){
        TLorentzVector recoJet;
        recoJet.SetPtEtaPhiM(
          reco_R4_pt->at(0),
          reco_R4_eta->at(0),
          reco_R4_phi->at(0),
          reco_R4_m->at(0)
        );

        //Plot the Delta R
        hist_DR_reco_cut_truth->Fill(truthJet.DeltaR(recoJet),evtw);
      }

      if(reco_R4_pt->size()!=0){
        TLorentzVector recoJet;
        recoJet.SetPtEtaPhiM(
          reco_R4_pt->at(0),
          reco_R4_eta->at(0),
          reco_R4_phi->at(0),
          reco_R4_m->at(0)
        );

        //Plot the Delta R
        hist_DR_reco_truth->Fill(truthJet.DeltaR(recoJet),evtw);
      }

      if(track_R4_pt->size()!=0){
        TLorentzVector trackJet;
        trackJet.SetPtEtaPhiM(
          track_R4_pt->at(0),
          track_R4_eta->at(0),
          track_R4_phi->at(0),
          track_R4_m->at(0)
        );

        //Plot the Delta R
        hist_DR_track_truth->Fill(truthJet.DeltaR(trackJet),evtw);
      }
    }
  }

  std::cout << "Done DeltaR" << std::endl;

  hist_DR_reco_truth->SetMarkerStyle(20);
  hist_DR_reco_truth->SetMarkerColor(kRed);
  hist_DR_reco_truth->Scale(1/hist_DR_reco_truth->Integral());
  hist_DR_reco_truth->DrawNormalized("");

  hist_DR_reco_cut_truth->SetMarkerStyle(21);
  hist_DR_reco_cut_truth->SetMarkerColor(kBlue);
  hist_DR_reco_cut_truth->Scale(1/hist_DR_reco_cut_truth->Integral());
  hist_DR_reco_cut_truth->DrawNormalized("Same");

  hist_DR_track_truth->SetMarkerStyle(22);
  hist_DR_track_truth->SetMarkerColor(kOrange);
  hist_DR_track_truth->Scale(1/hist_DR_track_truth->Integral());
  hist_DR_track_truth->DrawNormalized("Same");
  canvas->Print("DeltaR.svg");
  canvas->Clear();

  //--------------------------------------------------------------------------//
  //                               Matched Jets                               //
  //--------------------------------------------------------------------------//

  TH1F *hist_matched_reco_truth_20 = new TH1F(
    "Matched Delta R [reco] > 20 GeV",
    "Matched Delta R; pTjet/pTtruth; Events",
    20,
    0,
    5
  );

  TH1F *hist_matched_reco_truth_100 = new TH1F(
    "Matched Delta R [reco] > 100 GeV",
    "Matched Delta R; pTjet/pTtruth; Events",
    20,
    0,
    5
  );

  TH1F *hist_matched_reco_truth_500 = new TH1F(
    "Matched Delta R [reco] > 500 GeV",
    "Matched Delta R; pTjet/pTtruth; Events",
    20,
    0,
    5
  );

  TH1F *hist_matched_track_truth_20 = new TH1F(
    "Matched Delta R [track] > 20 GeV",
    "Matched Delta R; pTjet/pTtruth; Events",
    20,
    0,
    5
  );

  TH1F *hist_matched_track_truth_100 = new TH1F(
    "Matched Delta R [track] > 100 GeV",
    "Matched Delta R; pTjet/pTtruth; Events",
    20,
    0,
    5
  );

  TH1F *hist_matched_track_truth_500 = new TH1F(
    "Matched Delta R [track] > 500 GeV",
    "Matched Delta R; pTjet/pTtruth; Events",
    20,
    0,
    5
  );

  for (i = 0; i < nentries; i++) {
    nbytes = tree->GetEntry(i);

    if(truth_R4_pt->size()!=0){
      TLorentzVector truthJet;
      truthJet.SetPtEtaPhiM(
        truth_R4_pt->at(0),
        truth_R4_eta->at(0),
        truth_R4_phi->at(0),
        truth_R4_m->at(0)
      );

      if(reco_R4_pt->size() != 0){
        TLorentzVector recoJet;
        recoJet.SetPtEtaPhiM(
          reco_R4_pt->at(0),
          reco_R4_eta->at(0),
          reco_R4_phi->at(0),
          reco_R4_m->at(0)
        );

        if (truthJet.DeltaR(recoJet) < 0.3) {
          if (truth_R4_pt->at(0) > 20000.) {
            hist_matched_reco_truth_20->Fill(
              truth_R4_pt->at(0) / reco_R4_pt->at(0),
              evtw
            );
          }

          if (truth_R4_pt->at(0) > 100000.) {
            hist_matched_reco_truth_100->Fill(
              truth_R4_pt->at(0) / reco_R4_pt->at(0),
              evtw
            );
          }

          if (truth_R4_pt->at(0) > 500000.) {
            hist_matched_reco_truth_500->Fill(
              truth_R4_pt->at(0) / reco_R4_pt->at(0),
              evtw
            );
          }
        }
      }


      if(track_R4_pt->size() != 0){
        TLorentzVector trackJet;
        trackJet.SetPtEtaPhiM(
          track_R4_pt->at(0),
          track_R4_eta->at(0),
          track_R4_phi->at(0),
          track_R4_m->at(0)
        );

        //Plot the Delta R
        if (truthJet.DeltaR(trackJet) < 0.3) {
          if (truth_R4_pt->at(0)>20000.) {
            hist_matched_track_truth_20->Fill(
              truth_R4_pt->at(0) / track_R4_pt->at(0),
              evtw
            );
          }

          if (truth_R4_pt->at(0)>100000.) {
            hist_matched_track_truth_100->Fill(
              truth_R4_pt->at(0) / track_R4_pt->at(0),
              evtw
            );
          }

          if (truth_R4_pt->at(0)>500000.) {
            hist_matched_track_truth_500->Fill(
              truth_R4_pt->at(0) / track_R4_pt->at(0),
              evtw
            );
          }
        }
      }

    }
  }

  std::cout << "Done Matched Jets" << std::endl;

  hist_matched_reco_truth_20->SetMarkerStyle(20);
  hist_matched_reco_truth_20->SetMarkerColor(kRed);
  hist_matched_reco_truth_20->Draw("");

  hist_matched_reco_truth_100->SetMarkerStyle(21);
  hist_matched_reco_truth_100->SetMarkerColor(kBlue);
  hist_matched_reco_truth_100->Draw("Same");

  hist_matched_reco_truth_500->SetMarkerStyle(22);
  hist_matched_reco_truth_500->SetMarkerColor(kOrange);
  hist_matched_reco_truth_500->Draw("Same");
  canvas->Print("Matched_reco.svg");
  canvas->Clear();

  hist_matched_track_truth_20->SetMarkerStyle(20);
  hist_matched_track_truth_20->SetMarkerColor(kRed);
  hist_matched_track_truth_20->Draw("");

  hist_matched_track_truth_100->SetMarkerStyle(21);
  hist_matched_track_truth_100->SetMarkerColor(kBlue);
  hist_matched_track_truth_100->Draw("Same");

  hist_matched_track_truth_500->SetMarkerStyle(22);
  hist_matched_track_truth_500->SetMarkerColor(kOrange);
  hist_matched_track_truth_500->Draw("Same");
  canvas->Print("Matched_track.svg");
  canvas->Clear();

  //--------------------------------------------------------------------------//
  //                                Large Jets                                //
  //--------------------------------------------------------------------------//

  TH1F *hist_reco_R10_pt = new TH1F(
    "Leading jet reco R10 pT",
    "Leading jet reco R10 pT R10; pT(GeV); Events",
    20,
    0,
    100
  );

  TH1F *hist_reco_R10_trimmed_pt = new TH1F(
    "Leading jet reco R10 pT [Trimmed]",
    "Leading jet reco R10 pT R10; pT(GeV); Events",
    20,
    0,
    100
  );

  TH1F *hist_truth_R10_pt = new TH1F(
    "Leading jet truth R10 pT",
    "Leading jet truth R10 pT; pT(GeV); Events",
    20,
    0,
    100
  );

  TH1F *hist_truth_R10_trimmed_pt = new TH1F(
    "Leading jet truth R10 pT [Trimmed]",
    "Leading jet truth R10 pT; pT(GeV); Events",
    20,
    0,
    100
  );

  TH1F *hist_reco_R10_m = new TH1F(
    "Leading jet reco R10 mass",
    "Leading jet reco R10 mass; Mass; Events",
    20,
    0,
    100
  );

  TH1F *hist_reco_R10_trimmed_m = new TH1F(
    "Leading jet reco R10 mass [Trimmed]",
    "Leading jet reco R10 mass; Mass; Events",
    20,
    0,
    100
  );

  TH1F *hist_truth_R10_m = new TH1F(
    "Leading jet truth R10 mass",
    "Leading jet truth R10 mass; Mass; Events",
    20,
    0,
    100
  );

  TH1F *hist_truth_R10_trimmed_m = new TH1F(
    "Leading jet truth R10 mass [Trimmed]",
    "Leading jet truth R10 mass; Mass; Events",
    20,
    0,
    100
  );

  for (i = 0; i < nentries; i++) {
    nbytes = tree->GetEntry(i);

    if(reco_R10_pt->size() > 0){
      hist_reco_R10_pt->Fill(reco_R10_pt->at(0)/1000.,evtw);
      hist_reco_R10_trimmed_pt->Fill(reco_R10_trimmed_pt->at(0)/1000., evtw);

      hist_reco_R10_m->Fill(reco_R10_m->at(0)/1000., evtw);
      hist_reco_R10_trimmed_m->Fill(reco_R10_trimmed_m->at(0)/1000., evtw);
    }

    if (truth_R10_m->size() > 0) {
      hist_truth_R10_pt->Fill(truth_R10_pt->at(0)/1000.,evtw);
      hist_truth_R10_trimmed_pt->Fill(truth_R10_trimmed_pt->at(0)/1000., evtw);

      hist_truth_R10_m->Fill(truth_R10_m->at(0)/1000.);
      hist_truth_R10_trimmed_m->Fill(truth_R10_trimmed_m->at(0)/1000., evtw);
    }
  }

  std::cout << "Done Large Jets" << std::endl;


  // pT plots
  hist_reco_R10_pt->SetMarkerStyle(20);
  hist_reco_R10_pt->SetMarkerColor(kRed);
  hist_reco_R10_pt->Draw("");

  hist_reco_R10_trimmed_pt->SetMarkerStyle(21);
  hist_reco_R10_trimmed_pt->SetMarkerColor(kBlue);
  hist_reco_R10_trimmed_pt->Draw("Same");
  canvas->Print("reco_R10_pt.svg");
  canvas->Clear();

  hist_truth_R10_pt->SetMarkerStyle(20);
  hist_truth_R10_pt->SetMarkerColor(kRed);
  hist_truth_R10_pt->Draw("");

  hist_truth_R10_trimmed_pt->SetMarkerStyle(21);
  hist_truth_R10_trimmed_pt->SetMarkerColor(kBlue);
  hist_truth_R10_trimmed_pt->Draw("Same");
  canvas->Print("truth_R10_pt.svg");
  canvas->Clear();

  // Mass plots
  hist_reco_R10_m->SetMarkerStyle(20);
  hist_reco_R10_m->SetMarkerColor(kRed);
  hist_reco_R10_m->Draw("");

  hist_reco_R10_trimmed_m->SetMarkerStyle(21);
  hist_reco_R10_trimmed_m->SetMarkerColor(kBlue);
  hist_reco_R10_trimmed_m->Draw("Same");
  canvas->Print("reco_R10_m.svg");
  canvas->Clear();

  hist_truth_R10_m->SetMarkerStyle(20);
  hist_truth_R10_m->SetMarkerColor(kRed);
  hist_truth_R10_m->Draw("");

  hist_truth_R10_trimmed_m->SetMarkerStyle(21);
  hist_truth_R10_trimmed_m->SetMarkerColor(kBlue);
  hist_truth_R10_trimmed_m->Draw("Same");
  canvas->Print("truth_R10_m.svg");
}
