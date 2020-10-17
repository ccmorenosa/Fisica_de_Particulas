#include <iostream>
#include <string>
#include <stdio.h>

void ex_3 () {
  // Open de Root File
  TFile *file = TFile::Open("../ROOT_Files/Data_8TeV.root");

  // Get the Tree from the Root File
  TTree *tree = (TTree*) file->Get("mini");
  tree->Print();

  Bool_t e_trig;
  Bool_t mu_trig;
  Bool_t good_vtx;

  UInt_t lep_n;
  UInt_t jet_n;

  Float_t MET;
  Float_t MET_phi;

  Float_t lep_pt[10];
  Float_t lep_eta[10];
  Float_t lep_phi[10];
  Float_t lep_E[10];

  Int_t lep_type[10];

  Float_t lep_ptcone30[10];
  Float_t lep_etcone20[10];

  Float_t jet_pt[10];
  Float_t jet_eta[10];
  Float_t jet_jvf[10];
  Float_t jet_mv1[10];

  tree->SetBranchAddress("trigE", &e_trig);
  tree->SetBranchAddress("trigM", &mu_trig);
  tree->SetBranchAddress("hasGoodVertex", &good_vtx);
  tree->SetBranchAddress("lep_n", &lep_n);
  tree->SetBranchAddress("jet_n", &jet_n);
  tree->SetBranchAddress("met_et", &MET);
  tree->SetBranchAddress("met_phi", &MET_phi);

  tree->SetBranchAddress("lep_pt", &lep_pt);
  tree->SetBranchAddress("lep_eta", &lep_eta);
  tree->SetBranchAddress("lep_phi", &lep_phi);
  tree->SetBranchAddress("lep_E", &lep_E);
  tree->SetBranchAddress("lep_type", &lep_type);
  tree->SetBranchAddress("lep_ptcone30", &lep_ptcone30);
  tree->SetBranchAddress("lep_etcone20", &lep_etcone20);

  tree->SetBranchAddress("jet_pt", &jet_pt);
  tree->SetBranchAddress("jet_eta", &jet_eta);
  tree->SetBranchAddress("jet_jvf", &jet_jvf);
  tree->SetBranchAddress("jet_MV1", &jet_mv1);

  TCanvas *canvas = new TCanvas("Canvas", "", 800, 600);

  // Histogram for cuts
  TH1F *cutflow = new TH1F("Cutflow","Cutflow; Cut; Events",10,0,10);

  // Histograms for Leptons
  TH1F *hist_lep_pt = new TH1F(
    "Leptons pT", "Leptons pT; pT (GeV); Events",
    50, 0, 1000
  );

  TH1F *hist_lep_pt_cone30 = new TH1F(
    "Track isolation", "Track isolation; lep_ptcone30/lep_pt; Events",
    50, 0, 0.2
  );

  TH1F *hist_lep_et_cone20 = new TH1F(
    "Calorimeter isolation",
    "Calorimeter isolation; lep_etcone30/lep_pt; Events",
    50, 0, 0.2
  );

  TH1F *hist_lep_eta = new TH1F(
    "Leptons eta", "Leptons eta; eta; Events",
    50, -5, 5
  );

  // Histograms for Jets
  TH1F *hist_njets = new TH1F(
    "Number of jets", "n-jets; Jet multiplicity; Events",
    10,0,10
  );

  TH1F *hist_jet_pt = new TH1F(
    "Jets pT", "Jets pT; pT (GeV); Events",
    50, 0, 1000
  );

  TH1F *hist_jet_eta = new TH1F(
    "Jets eta", "Jets eta; eta; Events",
    50, -5, 5
  );

  TH1F *hist_jet_JVF = new TH1F(
    "Jets JVF", "Jets JVF; JVF; Events",
    50, 0, 1
  );

  TH1F *hist_jet_MV1 = new TH1F(
    "Jets MV1", "Jets MV1; MV1; Events",
    50, 0.5, 1
  );

  // Histogram for bJets
  TH1F *hist_nbjets = new TH1F(
    "Number of b-jets", "n-bjets; b-Jet multiplicity; Events",
    50, 0, 6
  );

  // Histogram for MET
  TH1F *hist_MET = new TH1F(
    "MET", "MET; MET (GeV); Events",
    50, 0, 200
  );

  // Histogram for mTW
  TH1F *hist_mTW = new TH1F(
    "mTW", "mTW; mTW (GeV); Events",
    50, 0, 200
  );

  int nentries, nbytes, i;
  nentries = (Int_t)tree->GetEntries();

  int cut1 = 0;
  int cut2 = 0;
  int cut3 = 0;
  int cut4 = 0;
  int cut5 = 0;
  int cut6 = 0;
  int cut7 = 0;
  int cut8 = 0;

  for (i = 0; i < nentries; i++) {
    nbytes = tree->GetEntry(i);

    //First cut: Good vertex
    if(!good_vtx) continue;
    cut1++;
    cutflow->Fill(1);

    //Second cut: Trigger
    if(!e_trig && !mu_trig) continue;
    cut2++;
    cutflow->Fill(2);

    // Preselection of good leptons
    int n_mu = 0;
    int n_el = 0;
    int n_lep = 0;
    int g_lep = 0;

    //Loop over leptons
    for (unsigned int j=0; j < lep_n; j++) {
      hist_lep_pt->Fill(lep_pt[j] / 1000);
      hist_lep_pt_cone30->Fill(lep_ptcone30[j] / lep_pt[j]);
      hist_lep_et_cone20->Fill(lep_etcone20[j] / lep_pt[j]);
      hist_lep_eta->Fill(lep_eta[j]);

      if( lep_pt[j] < 25000.0) continue;

      if( lep_ptcone30[j]/lep_pt[j] > 0.15 ) continue;

      if( lep_etcone20[j]/lep_pt[j] > 0.15 ) continue;

      if(lep_type[j] == 13 && TMath::Abs(lep_eta[j]) < 2.5){
        n_mu++;
        g_lep = j;
      }

      if (lep_type[j]==11 && TMath::Abs(lep_eta[j]) < 2.47){
        if (TMath::Abs(lep_eta[j]) < 1.37 || TMath::Abs(lep_eta[j]) > 1.52) {
          n_el++;
          g_lep = j;
        }
      }

    }
    n_lep = n_el + n_mu;

    //Select events with only 1 good lepton and fill the cutflow histogram
    if (n_lep != 1) continue;
    cutflow->Fill(3);
    cut3++;


    int n_jets = 0;
    int n_bjets = 0;

    //Number of jets distribution
    hist_njets->Fill(jet_n);

    //Fourth cut: At least 4 jets
    if (jet_n < 4) continue;
    cutflow->Fill(4);
    cut4++;

    for (unsigned int j=0; j < jet_n; j++){
      hist_jet_pt->Fill(jet_pt[j] / 1000);
      hist_jet_eta->Fill(jet_eta[j]);
      hist_jet_JVF->Fill(jet_jvf[j]);
      hist_jet_MV1->Fill(jet_mv1[j]);
      // To complete: apply jet cuts to find the good jets
      if (jet_pt[j] < 25000.) continue;

      //Eta cut
      if (jet_eta[j] > 2.5) continue;
      // JVF cleaning
      if (jet_pt[j] < 50000. && TMath::Abs(jet_eta[j]) < 2.4) {
        if (jet_jvf[j] < 0.5) continue;
      }

      n_jets ++;

      // cut on 0.7892 MV1 and count the number of b-jets
      if (jet_mv1[j] < 0.7892) continue;

      n_bjets ++;
    }


    //Fifth cut: At least 4 good jets
    if (n_jets < 4) continue;
    cutflow->Fill(5);
    cut5++;

    //Sixth cut: at least one b-jet
    hist_nbjets->Fill(n_bjets);
    if (n_bjets < 2) continue;
    cutflow->Fill(6);
    cut6++;

    //Seventh cut: MET > 30 GeV
    hist_MET->Fill(MET / 1000);
    if (MET < 30000.) continue;
    cutflow->Fill(7);
    cut7++;

    // TLorentzVector definitions
    TLorentzVector Lepton  = TLorentzVector();
    TLorentzVector  MeT  = TLorentzVector();

    //To complete: Lorentz vectors for the lepton and MET. Use SetPtEtaPhiE().
    Lepton.SetPtEtaPhiE(
      lep_pt[g_lep],
      lep_eta[g_lep],
      lep_phi[g_lep],
      lep_E[g_lep]
    );

    MeT.SetPtEtaPhiE(
      MET,
      0,
      MET_phi,
      MET
    );

    //Calculation of the mTW using TLorentz vectors
    float mTW = sqrt(2*Lepton.Pt()*MeT.Et()*(1-cos(Lepton.DeltaPhi(MeT))));

    //Eight cut: mTW > 30 GeV
    hist_mTW->Fill(mTW / 1000);
    if (mTW < 30000.) continue;
    cutflow->Fill(8);
    cut8++;
  }

  std::cout << "Done!" << std::endl;
  std::cout << "All events:" << nentries << std::endl;
  std::cout << "Cut1:" << cut1 << std::endl;
  std::cout << "Cut2:" << cut2 << std::endl;
  std::cout << "Cut3:" << cut3 << std::endl;
  std::cout << "Cut4:" << cut4 << std::endl;
  std::cout << "Cut5:" << cut5 << std::endl;
  std::cout << "Cut6:" << cut6 << std::endl;
  std::cout << "Cut7:" << cut7 << std::endl;
  std::cout << "Cut8:" << cut8 << std::endl;

  std::cout << '\n';

  std::cout << '\n';
  std::cout << "Saving Histograms" << '\n';
  std::cout << '\n';

  std::cout << '\n';
  std::cout << "Histogram for Cuts" << '\n';
  std::cout << '\n';

  canvas->SetLogy();
  cutflow->SetFillColor(kCyan);
  cutflow->Draw("");
  canvas->Print("cutflow.svg");
  canvas->Clear();

  std::cout << '\n';
  std::cout << "Histogram for Leptons" << '\n';
  std::cout << '\n';

  hist_lep_pt->SetFillColor(kCyan);
  hist_lep_pt->Draw();
  canvas->Print("hist_lep_pt.svg");
  canvas->Clear();

  hist_lep_pt_cone30->SetFillColor(kCyan);
  hist_lep_pt_cone30->Draw();
  canvas->Print("hist_lep_pt_cone30.svg");
  canvas->Clear();

  hist_lep_et_cone20->SetFillColor(kCyan);
  hist_lep_et_cone20->Draw();
  canvas->Print("hist_lep_et_cone20.svg");
  canvas->Clear();

  hist_lep_eta->SetFillColor(kCyan);
  hist_lep_eta->Draw();
  canvas->Print("hist_lep_eta.svg");
  canvas->Clear();

  std::cout << '\n';
  std::cout << "Histogram for Jets" << '\n';
  std::cout << '\n';

  hist_njets->SetFillColor(kCyan);
  hist_njets->Draw();
  canvas->Print("hist_njets.svg");
  canvas->Clear();

  hist_jet_pt->SetFillColor(kCyan);
  hist_jet_pt->Draw();
  canvas->Print("hist_jet_pt.svg");
  canvas->Clear();

  hist_jet_eta->SetFillColor(kCyan);
  hist_jet_eta->Draw();
  canvas->Print("hist_jet_eta.svg");
  canvas->Clear();

  hist_jet_JVF->SetFillColor(kCyan);
  hist_jet_JVF->Draw();
  canvas->Print("hist_jet_JVF.svg");
  canvas->Clear();

  hist_jet_MV1->SetFillColor(kCyan);
  hist_jet_MV1->Draw();
  canvas->Print("hist_jet_MV1.svg");
  canvas->Clear();

  std::cout << '\n';
  std::cout << "Histogram for b-Jets" << '\n';
  std::cout << '\n';

  hist_nbjets->SetFillColor(kCyan);
  hist_nbjets->Draw();
  canvas->Print("hist_nbjets.svg");
  canvas->Clear();
  canvas->SetLogy(false);

  std::cout << '\n';
  std::cout << "Histogram for MET" << '\n';
  std::cout << '\n';

  hist_MET->SetFillColor(kCyan);
  hist_MET->Draw();
  canvas->Print("hist_MET.svg");
  canvas->Clear();

  std::cout << '\n';
  std::cout << "Histogram for mTW" << '\n';
  std::cout << '\n';

  hist_mTW->SetFillColor(kCyan);
  hist_mTW->Draw();
  canvas->Print("hist_mTW.svg");
  canvas->Clear();

  std::cout << '\n';
  std::cout << '\n';
  std::cout << "Done" << '\n';
}
