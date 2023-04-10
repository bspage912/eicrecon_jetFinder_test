// Read edm4eic Trees

#include "fastjet/ClusterSequence.hh"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace fastjet;
using namespace std;


int main(int argc, char* argv[]) {

  //const int geometry = atoi(argv[1]);
  const char* infile = argv[1];

  TChain *mychain = new TChain("events");
  mychain->Add(infile);

  const char* output="test.hist.root";
  TFile *ofile = TFile::Open(output,"recreate");

  //TTreeReader tree_reader(tree);
  TTreeReader tree_reader(mychain);

  // Reco Jets
  TTreeReaderArray<int> recoType = {tree_reader, "ReconstructedJets.type"};
  TTreeReaderArray<float> recoNRG = {tree_reader, "ReconstructedJets.energy"};
  TTreeReaderArray<int> recoPDG = {tree_reader, "ReconstructedJets.PDG"};
  TTreeReaderArray<float> recoMomX = {tree_reader, "ReconstructedJets.momentum.x"};
  TTreeReaderArray<float> recoMomY = {tree_reader, "ReconstructedJets.momentum.y"};
  TTreeReaderArray<float> recoMomZ = {tree_reader, "ReconstructedJets.momentum.z"};
  TTreeReaderArray<float> recoM = {tree_reader, "ReconstructedJets.mass"};

  //TTreeReaderArray<int> recoType = {tree_reader, "GeneratedJets.type"};
  //TTreeReaderArray<float> recoNRG = {tree_reader, "GeneratedJets.energy"};
  //TTreeReaderArray<int> recoPDG = {tree_reader, "GeneratedJets.PDG"};
  //TTreeReaderArray<float> recoMomX = {tree_reader, "GeneratedJets.momentum.x"};
  //TTreeReaderArray<float> recoMomY = {tree_reader, "GeneratedJets.momentum.y"};
  //TTreeReaderArray<float> recoMomZ = {tree_reader, "GeneratedJets.momentum.z"};
  
  // MC
  TTreeReaderArray<int> mcGenStat = {tree_reader, "MCParticles.generatorStatus"};
  TTreeReaderArray<float> mcMomX = {tree_reader, "MCParticles.momentum.x"};
  TTreeReaderArray<float> mcMomY = {tree_reader, "MCParticles.momentum.y"};
  TTreeReaderArray<float> mcMomZ = {tree_reader, "MCParticles.momentum.z"};
  TTreeReaderArray<double> mcM = {tree_reader, "MCParticles.mass"};
  TTreeReaderArray<int> pdg = {tree_reader, "MCParticles.PDG"};

  // Reconstructed Particles
  TTreeReaderArray<float> recoPartMomX = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
  TTreeReaderArray<float> recoPartMomY = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
  TTreeReaderArray<float> recoPartMomZ = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
  TTreeReaderArray<float> recoPartM = {tree_reader, "ReconstructedChargedParticles.mass"};

  // Histograms
  // Internal Consistency Checks and Distributions
  TH1D *recoJetTypeHist = new TH1D("recoJetType","",5,0.,5);

  TH1D *numRecoJetsEventHist = new TH1D("numRecoJetsEvent","",20,0.,20.);
  TH2D *recoJetEvsEtaHist = new TH2D("recoJetEvsEta","",100,-5.,5.,300,0.,300.);
  TH2D *recoJetPhiVsEtaHist = new TH2D("recoJetPhiVsEta","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH1D *numRecoJetPartsHist = new TH1D("numRecoJetParts","",20,0.,20.);
  TH2D *recoJetPartEvsEtaHist = new TH2D("recoJetPartEvsEta","",100,-5.,5.,300,0.,300.);
  TH2D *recoJetPartPhiVsEtaHist = new TH2D("recoJetPartPhiVsEta","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH2D *recoJetEvsPartESumHist = new TH2D("recoJetEvsPartESum","",3000,0.,300.,3000,0.,300.);
  TH1D *recoJetEDiffHist = new TH1D("recoJetEDiff","",500,-10.,10.);

  TH2D *recoJetEvsEtaBadHist = new TH2D("recoJetEvsEtaBad","",100,-5.,5.,300,0.,300.);
  TH2D *recoJetPhiVsEtaBadHist = new TH2D("recoJetPhiVsEtaBad","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  // Comparison to Local Jets
  TH2D *numLocalVsRecoJetsHist = new TH2D("numLocalVsRecoJets","",20,0.,20.,20,0.,20.);
  TH2D *numLocalVsRecoJetPartsHist = new TH2D("numLocalVsRecoJetParts","",50,0.,50.,50,0.,50.);

  TH2D *localVsRecoJetPtHist = new TH2D("localVsRecoJetPt","",500,0.,50.,500,0.,50.);
  TH2D *localVsRecoJetEHist = new TH2D("localVsRecoJetE","",3000,0.,300.,3000,0.,300.);
  TH2D *localVsRecoJetEtaHist = new TH2D("localVsRecoJetEta","",1000,-5.,5.,1000,-5.,5.);
  TH2D *localVsRecoJetPhiHist = new TH2D("localVsRecoJetPhi","",500,-TMath::Pi(),TMath::Pi(),500,-TMath::Pi(),TMath::Pi());

  
  // Begin Loop Over Events
  int NEVENTS = 0;
  int TOTPARTS = 0;
  int MISMATCHJETE = 0;
  while(tree_reader.Next()) {

    if(NEVENTS%10000 == 0) cout << "Events Processed: " << NEVENTS << endl;

    // Reconstructed Particles for Clustering
    vector<PseudoJet> particles;

    int numJets = 0;
    cout << "Event " << NEVENTS << " " << recoType.GetSize() << endl;
    for(unsigned int i=0; i<recoType.GetSize(); i++) // Loop over Entries in ReconstructedJet
      {
	recoJetTypeHist->Fill(recoType[i]);

	int numParts = 0;
	TLorentzVector totVec;
	if(recoType[i] == 0) // Select Jets
	  {
	    TVector3 jetMom(recoMomX[i],recoMomY[i],recoMomZ[i]);
	    
	    cout << "Jet " << i << " Energy = " << recoNRG[i] << " Eta = " << jetMom.PseudoRapidity() << " Phi = " << jetMom.Phi() << " Number = " << recoPDG[i] << endl;

	    recoJetEvsEtaHist->Fill(jetMom.PseudoRapidity(),recoNRG[i]);
	    recoJetPhiVsEtaHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());

	    numJets++;

	    for(unsigned int j=0; j<recoType.GetSize(); j++) // For each jet, loop through all entries again
	      {
		if((recoType[j] == 1) && (recoPDG[j] == recoPDG[i])) // Select particles associated with jet
		  {
		    TLorentzVector local;
		    local.SetPxPyPzE(recoMomX[j],recoMomY[j],recoMomZ[j],recoNRG[j]);
		    totVec += local;

		    //cout << "    Particle " << j << " Energy = " << recoNRG[j] << " Eta = " << local.PseudoRapidity() << " Phi = " << local.Phi() << " Mass = " << recoM[j] << endl;

		    recoJetPartEvsEtaHist->Fill(local.PseudoRapidity(),recoNRG[j]);
		    recoJetPartPhiVsEtaHist->Fill(local.PseudoRapidity(),local.Phi());

		    numParts++;
		    TOTPARTS++;
		  }
	      }
	    cout << "    Number Of Particles = " << numParts << " EnergySum = " << totVec.E() << " Energy Diff = " << totVec.E() - recoNRG[i] << endl;

	    numRecoJetPartsHist->Fill(numParts);
	    recoJetEvsPartESumHist->Fill(recoNRG[i],totVec.E());
	    recoJetEDiffHist->Fill(recoNRG[i]-totVec.E());

	    if(TMath::Abs(totVec.E() - recoNRG[i]) > 0.00001)
	      {
		recoJetEvsEtaBadHist->Fill(jetMom.PseudoRapidity(),recoNRG[i]);
		recoJetPhiVsEtaBadHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());
	      }
	  }
      }
    cout << endl;

    numRecoJetsEventHist->Fill(numJets);


    // Print Reconstructed and Generated Particle Information
    for(unsigned int i=0; i<recoPartMomX.GetSize(); i++)
      {
	TVector3 partMom(recoPartMomX[i],recoPartMomY[i],recoPartMomZ[i]);
	double E = std::sqrt(recoPartMomX[i]*recoPartMomX[i]+recoPartMomY[i]*recoPartMomY[i]+recoPartMomZ[i]*recoPartMomZ[i]+recoPartM[i]*recoPartM[i]);
	
	//cout << "Reco Particle " << i << " Energy = " << E << " Eta = " << partMom.PseudoRapidity() << " Phi = " << partMom.Phi() << " Mass = " << recoPartM[i] << " Pt = " << partMom.Perp() << endl;

	if(1) //partMom.Perp() > 0.2 && partMom.Perp() < 100.0)
	  {
	    particles.push_back( PseudoJet(partMom.Px(),partMom.Py(),partMom.Pz(),E) );
	  }
      }
    //cout << endl;
    for(unsigned int i=0; i<mcMomX.GetSize(); i++)
      {
	if(mcGenStat[i] == 1)
	  {
	    TVector3 jetMom(mcMomX[i],mcMomY[i],mcMomZ[i]);
	    double E = std::sqrt(mcMomX[i]*mcMomX[i] + mcMomY[i]*mcMomY[i] + mcMomZ[i]*mcMomZ[i] + mcM[i]*mcM[i]);
	    
	    //cout << "Particle " << i << " Energy = " << E << " Eta = " << jetMom.PseudoRapidity() << " Phi = " << jetMom.Phi() << " PDG = " << pdg[i] << endl;
	  }
      }
    //cout << endl;


    // Find Jets from Reconstructed Particles
    // Define Jet
    JetDefinition jet_def_akt(antikt_algorithm,1.0);

    // Cluster in Beam Frame
    ClusterSequence cs_akt_all_beam(particles, jet_def_akt);

    // Set Min Jet Pt
    //double ptmin = 1.0;

    // Get Jets
    vector<PseudoJet> jets_akt_all_beam = sorted_by_pt(cs_akt_all_beam.inclusive_jets());

    // Compare to Jets from EICrecon
    numLocalVsRecoJetsHist->Fill(numJets,jets_akt_all_beam.size());

    // Loop Over Reco Jets
    for(unsigned int i=0; i<recoType.GetSize(); i++)
      {
	int numParts = 0;
	TLorentzVector totVec;
	if(recoType[i] == 0) // Select Jets
	  {
	    TVector3 jetMom(recoMomX[i],recoMomY[i],recoMomZ[i]);

	    for(unsigned int j=0; j<recoType.GetSize(); j++) // For each jet, loop through all entries again
	      {
		if((recoType[j] == 1) && (recoPDG[j] == recoPDG[i])) // Select particles associated with jet
		  {
		    TLorentzVector local;
		    local.SetPxPyPzE(recoMomX[j],recoMomY[j],recoMomZ[j],recoNRG[j]);
		    totVec += local;
		    
		    numParts++;
		  }
	      }

	    double minDist = 10000.0;
	    int matchIndex = -1;
	    for(unsigned int jn=0; jn<jets_akt_all_beam.size(); jn++)
	      {
		double dEta = jetMom.PseudoRapidity() - jets_akt_all_beam[jn].eta();
		double dPhi = TVector2::Phi_mpi_pi(jetMom.Phi() - jets_akt_all_beam[jn].phi());
		double dR = std::sqrt(dEta*dEta + dPhi*dPhi);

		if(dR < minDist)
		  {
		    minDist = dR;
		    matchIndex = jn;
		  }
	      }

	    if(matchIndex > -1)
	      {
		localVsRecoJetPtHist->Fill(jetMom.Perp(),jets_akt_all_beam[matchIndex].pt());
		localVsRecoJetEHist->Fill(recoNRG[i],jets_akt_all_beam[matchIndex].e());
		localVsRecoJetEtaHist->Fill(jetMom.PseudoRapidity(),jets_akt_all_beam[matchIndex].eta());
		localVsRecoJetPhiHist->Fill(jetMom.Phi(),jets_akt_all_beam[matchIndex].phi());

		// Compare Constituents
		vector<PseudoJet> cons = jets_akt_all_beam[matchIndex].constituents();

		int numLocalJetPart = 0;
		for(unsigned int jc=0; jc<cons.size(); jc++)
		  {
		    if(cons[jc].pt() > 0.2 && cons[jc].pt() < 100.0)
		      numLocalJetPart++;
		  }

		numLocalVsRecoJetPartsHist->Fill(numParts,numLocalJetPart);

		if(std::abs(recoNRG[i] - jets_akt_all_beam[matchIndex].e()) > 10e-6) 
		  {
		    cout << "Event " << NEVENTS << " Jet " << i << " has Mismatched Energy = " << recoNRG[i] - jets_akt_all_beam[matchIndex].e() << endl;
		    MISMATCHJETE++;
		  }
	      }
	  }
      }

    NEVENTS++;
  } // End Loop Over Events

  cout << "Total Particles = " << TOTPARTS << endl;
  cout << "Total Reco/Local Jet Pairs with Mismatched Energy = " << MISMATCHJETE << endl;

  //mcEta->Write();
  //mcNRG->Write();

  ofile->Write();
  ofile->Close();
  return 0;
}
