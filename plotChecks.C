TFile *fa=0;

TCanvas *c[30];

void plotChecks()
{
  //gStyle->SetPalette(1,0);
  gStyle->SetOptStat(0);
  //gStyle->SetOptStat(1111);
  //gStyle->SetOptFit(1111);

  TString rootHistFnameA;
  TString rootHistFnameB;
  TString rootHistFnameC;
  TString rootHistFnameD;
  TString rootHistFnameE;

  // For Res Plots
  
  rootHistFnameA="test.hist.root";

  fa=new TFile(rootHistFnameA);
  assert(fa->IsOpen());
}


void plotRecoChecks()
{
  c[0]=new TCanvas("c0","Reco Jet Internal Consistency",800,600);
  c[1]=new TCanvas("c1","Reco-Local Comparison: Number of Jets",800,600);
  c[2]=new TCanvas("c2","Reco-Local Comparison: Number of Jet Particles",800,600);
  c[3]=new TCanvas("c3","Reco-Local Comparison: Kinematics",800,600);

  for(int i=0; i<4; i++)
    {
      c[i]->Clear();
    }

  c[0]->Divide(1,1);
  c[1]->Divide(1,1);
  c[2]->Divide(1,1);
  c[3]->Divide(2,2);

  TH2D *hA=(TH2D *)fa->Get("recoJetEvsPartESum");
  TH2D *hB=(TH2D *)fa->Get("numLocalVsRecoJets");
  TH2D *hC=(TH2D *)fa->Get("numLocalVsRecoJetParts");

  TH2D *hD[4];
  hD[0]=(TH2D *)fa->Get("localVsRecoJetE");
  hD[1]=(TH2D *)fa->Get("localVsRecoJetPt");
  hD[2]=(TH2D *)fa->Get("localVsRecoJetEta");
  hD[3]=(TH2D *)fa->Get("localVsRecoJetPhi");

  // Draw Histos
  c[0]->cd(1);
  hA->Draw("COLZ");
  hA->SetTitle("Reco Consistency Check: Sum Particle E Vs Jet E;Jet E;Sum Particle E");

  c[1]->cd(1);
  hB->Draw("COLZ");
  hB->SetTitle("Reco-Local Comp: Number of Jets;Num Reco Jets;Num Local Jets");

  c[2]->cd(1);
  hC->Draw("COLZ");
  hC->GetXaxis()->SetRangeUser(0,20);
  hC->GetYaxis()->SetRangeUser(0,20);
  hC->SetTitle("Reco-Local Comp: Number of Jet Particles;Num Reco Parts;Num Local Parts");

  for(int i=0; i<4; i++)
    {
      c[3]->cd(i+1);
      hD[i]->Draw("COLZ");
      if(i==0) hD[i]->SetTitle("Reco-Local Comp: Jet Energy;Reco E;Local E");
      if(i==1) hD[i]->SetTitle("Reco-Local Comp: Jet Pt;Reco Pt;Local Pt");
      if(i==2) hD[i]->SetTitle("Reco-Local Comp: Jet Eta;Reco Eta;Local Eta");
      if(i==3) hD[i]->SetTitle("Reco-Local Comp: Jet Phi;Reco Phi;Local Phi");
    }
}


void plotGenChecks()
{
  c[0]=new TCanvas("c0","Gen Jet Internal Consistency",800,600);
  c[1]=new TCanvas("c1","Gen-Local Comparison: Number of Jets",800,600);
  c[2]=new TCanvas("c2","Gen-Local Comparison: Number of Jet Particles",800,600);
  c[3]=new TCanvas("c3","Gen-Local Comparison: Kinematics",800,600);

  for(int i=0; i<4; i++)
    {
      c[i]->Clear();
    }

  c[0]->Divide(1,1);
  c[1]->Divide(1,1);
  c[2]->Divide(1,1);
  c[3]->Divide(2,2);

  TH2D *hA=(TH2D *)fa->Get("genJetEvsPartESum");
  TH2D *hB=(TH2D *)fa->Get("numLocalVsGenJets");
  TH2D *hC=(TH2D *)fa->Get("numLocalVsGenJetParts");

  TH2D *hD[4];
  hD[0]=(TH2D *)fa->Get("localVsGenJetE");
  hD[1]=(TH2D *)fa->Get("localVsGenJetPt");
  hD[2]=(TH2D *)fa->Get("localVsGenJetEta");
  hD[3]=(TH2D *)fa->Get("localVsGenJetPhi");

  // Draw Histos
  c[0]->cd(1);
  hA->Draw("COLZ");
  hA->SetTitle("Gen Consistency Check: Sum Particle E Vs Jet E;Jet E;Sum Particle E");

  c[1]->cd(1);
  hB->Draw("COLZ");
  hB->SetTitle("Gen-Local Comp: Number of Jets;Num Gen Jets;Num Local Jets");

  c[2]->cd(1);
  hC->Draw("COLZ");
  hC->GetXaxis()->SetRangeUser(0,20);
  hC->GetYaxis()->SetRangeUser(0,20);
  hC->SetTitle("Gen-Local Comp: Number of Jet Particles;Num Gen Parts;Num Local Parts");

  for(int i=0; i<4; i++)
    {
      c[3]->cd(i+1);
      hD[i]->Draw("COLZ");
      if(i==0) hD[i]->SetTitle("Gen-Local Comp: Jet Energy;Gen E;Local E");
      if(i==1) hD[i]->SetTitle("Gen-Local Comp: Jet Pt;Gen Pt;Local Pt");
      if(i==2) hD[i]->SetTitle("Gen-Local Comp: Jet Eta;Gen Eta;Local Eta");
      if(i==3) hD[i]->SetTitle("Gen-Local Comp: Jet Phi;Gen Phi;Local Phi");
    }
}


void plotRecoTruth()
{
  c[0]=new TCanvas("c0","Reco-Truth Delta R",800,600);
  c[1]=new TCanvas("c1","Duplicate Tracks",800,600);
  c[2]=new TCanvas("c2","Reco-Truth Comparison: Kinematics",800,600);
  c[3]=new TCanvas("c3","Reco-Truth Comparison: E Res and Duplicate Tracks",800,600);

  for(int i=0; i<4; i++)
    {
      c[i]->Clear();
    }

  c[0]->Divide(1,1);
  c[1]->Divide(1,1);
  c[2]->Divide(2,2);
  c[3]->Divide(1,1);

  TH1D *hA=(TH1D *)fa->Get("localGenRecoDetlaR");

  TH1D *hB[2];
  hB[0]=(TH1D *)fa->Get("localPartCombDR");
  hB[1]=(TH1D *)fa->Get("localGenPartCombDR");

  TH2D *hC[3];
  hC[0]=(TH2D *)fa->Get("localGenVsRecoJetE");
  hC[1]=(TH2D *)fa->Get("localGenVsRecoJetEta");
  hC[2]=(TH2D *)fa->Get("localGenVsRecoJetPhi");

  TH1D *hD[2];
  hD[0]=(TH1D *)fa->Get("localGenRecoJetEDiff");
  hD[1]=(TH1D *)fa->Get("localGenRecoJetEDiffNoDupTrack");
  

  // Draw Histos
  c[0]->cd(1);
  hA->Draw("HIST");
  hA->SetTitle("Delta R Between Closest Truth and Reco Jets;Delta R");
  gPad->SetLogy();

  c[1]->cd(1);
  hB[0]->Draw("HIST");
  hB[1]->Draw("HISTSAME");
  hB[1]->SetLineColor(kRed);
  hB[0]->GetXaxis()->SetRangeUser(0,2);
  hB[0]->SetTitle("Pairwise Distance Between Particles in a Jet;Delta R");

  TLegend *leg1 = new TLegend(0.50,0.88,0.88,0.65);
  leg1->SetBorderSize(0);
  leg1->AddEntry(hB[0],"Reconstructed Jets","l");
  leg1->AddEntry(hB[1],"Generated Jets","l");
  leg1->Draw();

  for(int i=0; i<3; i++)
    {
      c[2]->cd(i+1);
      hC[i]->Draw("COLZ");
      if(i==0){
	hC[i]->GetXaxis()->SetRangeUser(0,150);
	hC[i]->GetYaxis()->SetRangeUser(0,150);
      }
      if(i==0) hC[i]->SetTitle("Matched Generated Vs Reconstructed Jet Energy;Reco E;Gen E");
      if(i==1) hC[i]->SetTitle("Matched Generated Vs Reconstructed Jet Eta;Reco Eta;Gen Eta");
      if(i==2) hC[i]->SetTitle("Matched Generated Vs Reconstructed Jet Phi;Reco Phi;Gen Phi");
    }

  c[3]->cd(1);
  hD[0]->Draw("HIST");
  hD[1]->Draw("HISTSAME");
  hD[1]->SetLineColor(kRed);
  hD[0]->GetXaxis()->SetRangeUser(-1,3);
  hD[0]->SetTitle("(Reco - Truth)/Truth Jet Energy;Jet Res");
  gPad->SetLogy();

  TLegend *leg3 = new TLegend(0.50,0.88,0.88,0.65);
  leg3->SetBorderSize(0);
  leg3->AddEntry(hD[0],"All Jets","l");
  leg3->AddEntry(hD[1],"Jets w/o Duplicate Tracks","l");
  leg3->Draw();
}

