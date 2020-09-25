#include <iostream>
#include <string>
#include <stdio.h>

#pragma link C++ enum  myns::myenum;

void ex_1 () {
  // Open de Root File
  TFile *file = TFile::Open("../ROOT_Files/Tracks_Clusters.root");

  // Get the Tree from the Root File
  TTree *tree = (TTree*) file->Get("JetRecoTree");
  tree->Print();

  // Pointers and variables to stores the branches data
  UInt_t npv = -1;
  Float_t mu_avg = -1;
  vector<float> *tracks_pt = {};
  vector<float> *clusters_pt = {};
  vector<float> *tracks_eta = {};
  vector<float> *tracks_phi = {};
  vector<float> *tracks_m = {};
  vector<float> *tracks_vtx = {};
  vector<float> *clusters_eta = {};
  vector<float> *clusters_phi = {};
  vector<float> *clusters_m = {};
  vector<float> *particles_pt = {};
  vector<float> *particles_eta = {};
  vector<float> *particles_phi = {};
  vector<float> *particles_m = {};
  vector<float> *particles_pdgID = {};

  // Get the branches address into the variables
  tree->SetBranchAddress("NPV", &npv);
  tree->SetBranchAddress("mu_average", &mu_avg);

  tree->SetBranchAddress("Tracks_pt", &tracks_pt);

  tree->SetBranchAddress("Clusters_pt", &clusters_pt);

  tree->SetBranchAddress("Tracks_eta", &tracks_eta);
  tree->SetBranchAddress("Tracks_phi", &tracks_phi);
  tree->SetBranchAddress("Tracks_m", &tracks_m);
  tree->SetBranchAddress("Tracks_vtx", &tracks_vtx);
  tree->SetBranchAddress("Clusters_eta", &clusters_eta);
  tree->SetBranchAddress("Clusters_phi", &clusters_phi);
  tree->SetBranchAddress("Clusters_m", &clusters_m);
  tree->SetBranchAddress("Particles_pt", &particles_pt);
  tree->SetBranchAddress("Particles_eta", &particles_eta);
  tree->SetBranchAddress("Particles_phi", &particles_phi);
  tree->SetBranchAddress("Particles_m", &particles_m);
  tree->SetBranchAddress("Particles_pdgID", &particles_pdgID);


  // Create the canvas to plot the histograms
  TCanvas *canvas = new TCanvas(
    "Canvas",
    "a first way to plot a variable",
    800,
    600
  );

  TH1F *hist_mu_avg = new TH1F(
    "mu_average",
    "Average events; mu_average ; Events ",
    50,
    1,
    90
  );

  int nentries, nbytes;
  nentries = (Int_t)tree->GetEntries();

  for (int i = 0; i < nentries; i++) {
    nbytes = tree->GetEntry(i);
    hist_mu_avg->Fill(mu_avg);
  }
  std::cout << "Done mu" << '\n';

  hist_mu_avg->SetFillColor(kRed);
  hist_mu_avg->Draw();
  canvas->Print("H_mu.svg");
  canvas->Clear();

  TH2F *hist_npv_mu = new TH2F(
    "NPV_mu",
    "Number of primary vertices vs Average events; NPV; mu_average ; Events ",
    50,
    1,
    50,
    50,
    1,
    90
  );

  for (int i = 0; i < nentries; i++) {
    nbytes = tree->GetEntry(i);
    hist_npv_mu->Fill(npv, mu_avg);
  }
  std::cout << "Done NPV - mu" << '\n';

  hist_npv_mu->SetFillColor(kMagenta);
  hist_npv_mu->Draw("LEGO1");
  canvas->Print("H_npv_mu.svg");
  canvas->Clear();

  TH2F *hist_npv_tracks = new TH2F(
    "NPV_tracks",
    "Number of primary vertices vs Number of tracks; NPV; NTracks ; Events ",
    50,
    1,
    50,
    50,
    0,
    1000
  );

  TH2F *hist_npv_clusters = new TH2F(
    "NPV_clustres",
    "Number of primary vertices of clusters; NPV; NClusters ; Events ",
    50,
    1,
    50,
    50,
    0,
    1000
  );

  TH2F *hist_mu_tracks = new TH2F(
    "mu_tracks",
    "Average events vs Number of tracks; NPV; NTracks ; Events ",
    50,
    1,
    50,
    50,
    0,
    1000
  );

  TH2F *hist_mu_clusters = new TH2F(
    "mu_clustres",
    "Average events vs Number of clusters; NPV; NClusters ; Events ",
    50,
    1,
    50,
    50,
    0,
    1000
  );

  for (int i = 0; i < nentries; i++) {
    nbytes = tree->GetEntry(i);
    hist_npv_tracks->Fill(npv, tracks_pt->size());
    hist_npv_clusters->Fill(npv, clusters_pt->size());
    hist_mu_tracks->Fill(mu_avg, tracks_pt->size());
    hist_mu_clusters->Fill(mu_avg, clusters_pt->size());
  }

  std::cout << "Done tracks/clusters" << '\n';
  hist_npv_tracks->SetFillColor(kRed);
  hist_npv_tracks->Draw("LEGO1");
  canvas->Print("H_npv_tracks.svg");
  canvas->Clear();

  std::cout << "Done tracks/clusters" << '\n';
  hist_npv_clusters->SetFillColor(kCyan);
  hist_npv_clusters->Draw("LEGO1");
  canvas->Print("H_npv_clusters.svg");
  canvas->Clear();

  std::cout << "Done tracks/clusters" << '\n';
  hist_mu_tracks->SetFillColor(kYellow);
  hist_mu_tracks->Draw("LEGO1");
  canvas->Print("H_mu_tracks.svg");
  canvas->Clear();

  std::cout << "Done tracks/clusters" << '\n';
  hist_mu_clusters->SetFillColor(kOrange);
  hist_mu_clusters->Draw("LEGO1");
  canvas->Print("H_mu_clusters.svg");


  TH1F *hist_track_pt = new TH1F( "Track_pt", "Track pt; pt ; Events ", 50, 400, 1500);
  TH1F *hist_track_eta = new TH1F( "Track_eta", "Track eta; eta ; Events ", 50, -3, 3);
  TH1F *hist_track_phi = new TH1F( "Track_phi", "Track phi; phi ; Events ", 50, -4, 4);
  TH1F *hist_track_m = new TH1F( "Track_m", "Track m; m ; Events ", 50, 0, 200);
  TH1F *hist_track_vtx = new TH1F( "Track_vtx", "Track vtx; vtx ; Events ", 50, -5, 40);

  TH1F *hist_clusters_pt = new TH1F( "Clusters_pt", "Clusters pt; pt ; Events ", 50, 100, 1500);
  TH1F *hist_clusters_eta = new TH1F( "Clusters_eta", "Clusters eta; eta ; Events ", 50, -5, 5);
  TH1F *hist_clusters_phi = new TH1F( "Clusters_phi", "Clusters phi; phi ; Events ", 50, -.1, .1);
  TH1F *hist_clusters_m = new TH1F( "Clusters_m", "Clusters m; m ; Events ", 50, -0.5, 0.5);

  TH1F *hist_particles_pt = new TH1F( "Particles_pt", "Particles pt; pt ; Events ", 50, -5, 1000);
  TH1F *hist_particles_eta = new TH1F( "Particles_eta", "Particles eta; eta ; Events ", 50, 6, -6);
  TH1F *hist_particles_phi = new TH1F( "Particles_phi", "Particles phi; phi ; Events ", 50, -4, 4);
  TH1F *hist_particles_m = new TH1F( "Particles_m", "Particles m; m ; Events ", 50, -10, 1600);
  TH1F *hist_particles_pdgID = new TH1F( "Particles_pdgID", "Particles pdgID; pdgID ; Events ", 50, -4000, 4000);

  for (int i = 0; i < nentries; i++) {
    nbytes = tree->GetEntry(i);
    for(int tr=0; tr<tracks_pt->size(); tr++) {
      hist_track_pt->Fill(tracks_pt->at(tr));
      hist_track_eta->Fill(tracks_eta->at(tr));
      hist_track_phi->Fill(tracks_phi->at(tr));
      hist_track_m->Fill(tracks_m->at(tr));
      hist_track_vtx->Fill(tracks_vtx->at(tr));
    }


    for(int tr=0; tr<clusters_pt->size(); tr++) {
      hist_clusters_pt->Fill(clusters_pt->at(tr));
      hist_clusters_eta->Fill(clusters_eta->at(tr));
      hist_clusters_phi->Fill(clusters_phi->at(tr));
      hist_clusters_m->Fill(clusters_m->at(tr));
    }


    for(int tr=0; tr<particles_pt->size(); tr++) {
      hist_particles_pt->Fill(particles_pt->at(tr));
      hist_particles_eta->Fill(particles_eta->at(tr));
      hist_particles_phi->Fill(particles_phi->at(tr));
      hist_particles_m->Fill(particles_m->at(tr));
      hist_particles_pdgID->Fill(particles_pdgID->at(tr));
    }
  }

  std::cout << "Begin Tracks" << '\n';
  hist_track_pt->SetFillColor(kCyan);
  hist_track_pt->Draw();
  canvas->Print("H_track_pt.svg");
  canvas->Clear();

  hist_track_eta->SetFillColor(kCyan);
  hist_track_eta->Draw();
  canvas->Print("H_track_eta.svg");
  canvas->Clear();

  hist_track_phi->SetFillColor(kCyan);
  hist_track_phi->Draw();
  canvas->Print("H_track_phi.svg");
  canvas->Clear();

  hist_track_m->SetFillColor(kCyan);
  hist_track_m->Draw();
  canvas->Print("H_track_m.svg");
  canvas->Clear();

  hist_track_vtx->SetFillColor(kCyan);
  hist_track_vtx->Draw();
  canvas->Print("H_track_vtx.svg");
  canvas->Clear();
  std::cout << "Done Tracks" << '\n';



  std::cout << "Begin Clusters" << '\n';
  hist_clusters_pt->SetFillColor(kCyan);
  hist_clusters_pt->Draw();
  canvas->Print("H_clusters_pt.svg");
  canvas->Clear();

  hist_clusters_eta->SetFillColor(kCyan);
  hist_clusters_eta->Draw();
  canvas->Print("H_clusters_eta.svg");
  canvas->Clear();

  hist_clusters_phi->SetFillColor(kCyan);
  hist_clusters_phi->Draw();
  canvas->Print("H_clusters_phi.svg");
  canvas->Clear();

  hist_clusters_m->SetFillColor(kCyan);
  hist_clusters_m->Draw();
  canvas->Print("H_clusters_m.svg");
  canvas->Clear();
  std::cout << "Done Clusters" << '\n';



  std::cout << "Begin Particles" << '\n';
  hist_particles_pt->SetFillColor(kCyan);
  hist_particles_pt->Draw();
  canvas->Print("H_particles_pt.svg");
  canvas->Clear();

  hist_particles_eta->SetFillColor(kCyan);
  hist_particles_eta->Draw();
  canvas->Print("H_particles_eta.svg");
  canvas->Clear();

  hist_particles_phi->SetFillColor(kCyan);
  hist_particles_phi->Draw();
  canvas->Print("H_particles_phi.svg");
  canvas->Clear();

  hist_particles_m->SetFillColor(kCyan);
  hist_particles_m->Draw();
  canvas->Print("H_particles_m.svg");
  canvas->Clear();

  hist_particles_pdgID->SetFillColor(kCyan);
  hist_particles_pdgID->Draw();
  canvas->Print("H_pdgID.svg");
  canvas->Clear();
  canvas->Close();
  std::cout << "Done Particles" << '\n';

}
