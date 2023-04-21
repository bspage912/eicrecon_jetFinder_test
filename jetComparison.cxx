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

  // Generated Jets
  TTreeReaderArray<int> genType = {tree_reader, "GeneratedJets.type"};
  TTreeReaderArray<float> genNRG = {tree_reader, "GeneratedJets.energy"};
  TTreeReaderArray<int> genPDG = {tree_reader, "GeneratedJets.PDG"};
  TTreeReaderArray<float> genMomX = {tree_reader, "GeneratedJets.momentum.x"};
  TTreeReaderArray<float> genMomY = {tree_reader, "GeneratedJets.momentum.y"};
  TTreeReaderArray<float> genMomZ = {tree_reader, "GeneratedJets.momentum.z"};
  TTreeReaderArray<float> genM = {tree_reader, "GeneratedJets.mass"};
  
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

  // Internal Consistency Checks and Distributions (Generated Jets)
  TH1D *numGenJetsEventHist = new TH1D("numGenJetsEvent","",20,0.,20.);
  TH2D *genJetEvsEtaHist = new TH2D("genJetEvsEta","",100,-5.,5.,300,0.,300.);
  TH2D *genJetPhiVsEtaHist = new TH2D("genJetPhiVsEta","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH1D *numGenJetPartsHist = new TH1D("numGenJetParts","",20,0.,20.);
  TH2D *genJetPartEvsEtaHist = new TH2D("genJetPartEvsEta","",100,-5.,5.,300,0.,300.);
  TH2D *genJetPartPhiVsEtaHist = new TH2D("genJetPartPhiVsEta","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH2D *genJetEvsPartESumHist = new TH2D("genJetEvsPartESum","",3000,0.,300.,3000,0.,300.);
  TH1D *genJetEDiffHist = new TH1D("genJetEDiff","",500,-10.,10.);

  TH2D *genJetEvsEtaBadHist = new TH2D("genJetEvsEtaBad","",100,-5.,5.,300,0.,300.);
  TH2D *genJetPhiVsEtaBadHist = new TH2D("genJetPhiVsEtaBad","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  // Comparison to Local Jets
  TH2D *numLocalVsRecoJetsHist = new TH2D("numLocalVsRecoJets","",20,0.,20.,20,0.,20.);
  TH2D *numLocalVsRecoJetPartsHist = new TH2D("numLocalVsRecoJetParts","",50,0.,50.,50,0.,50.);

  TH2D *localVsRecoJetPtHist = new TH2D("localVsRecoJetPt","",500,0.,50.,500,0.,50.);
  TH2D *localVsRecoJetEHist = new TH2D("localVsRecoJetE","",3000,0.,300.,3000,0.,300.);
  TH2D *localVsRecoJetEtaHist = new TH2D("localVsRecoJetEta","",1000,-5.,5.,1000,-5.,5.);
  TH2D *localVsRecoJetPhiHist = new TH2D("localVsRecoJetPhi","",500,-TMath::Pi(),TMath::Pi(),500,-TMath::Pi(),TMath::Pi());

  TH1D *recoPartMassHist = new TH1D("recoPartMass","",10000,0.,1.);

  TH1D *localPartCombDRHist = new TH1D("localPartCombDR","",1000,0.,10.);

  TH2D *numLocalVsGenJetsHist = new TH2D("numLocalVsGenJets","",20,0.,20.,20,0.,20.);
  TH2D *numLocalVsGenJetPartsHist = new TH2D("numLocalVsGenJetParts","",50,0.,50.,50,0.,50.);

  TH2D *localVsGenJetPtHist = new TH2D("localVsGenJetPt","",500,0.,50.,500,0.,50.);
  TH2D *localVsGenJetEHist = new TH2D("localVsGenJetE","",3000,0.,300.,3000,0.,300.);
  TH2D *localVsGenJetEtaHist = new TH2D("localVsGenJetEta","",1000,-5.,5.,1000,-5.,5.);
  TH2D *localVsGenJetPhiHist = new TH2D("localVsGenJetPhi","",500,-TMath::Pi(),TMath::Pi(),500,-TMath::Pi(),TMath::Pi());

  // Comparison between Local Generated and Reco Jets
  TH1D *localGenRecoDetlaRHist = new TH1D("localGenRecoDetlaR","",500,0.,5.);

  TH2D *localGenVsRecoJetEHist = new TH2D("localGenVsRecoJetE","",300,0.,300.,300,0.,300.);
  TH2D *localGenVsRecoJetEtaHist = new TH2D("localGenVsRecoJetEta","",100,-5.,5.,100,-5.,5.);
  TH2D *localGenVsRecoJetPhiHist = new TH2D("localGenVsRecoJetPhi","",100,-TMath::Pi(),TMath::Pi(),500,-TMath::Pi(),TMath::Pi());

  TH1D *localGenRecoJetEDiffHist = new TH1D("localGenRecoJetEDiff","",200,-10.,10.);
  TH2D *localGenRecoJetEDiffVsRecoEHist = new TH2D("localGenRecoJetEDiffVsRecoE","",300,0.,300.,200,-10.,10.);
  TH1D *localGenRecoJetEDiffDupTrackHist = new TH1D("localGenRecoJetEDiffDupTrack","",200,-10.,10.);
  TH1D *localGenRecoJetEDiffNoDupTrackHist = new TH1D("localGenRecoJetEDiffNoDupTrack","",200,-10.,10.);

  TH2D *localGenRecoJetEDiffVsDeltaRHist = new TH2D("localGenRecoJetEDiffVsDeltaR","",500,0.,5.,200,-10.,10.);

  TH1D *localGenPartCombDRHist = new TH1D("localGenPartCombDR","",1000,0.,10.);

  
  // Begin Loop Over Events
  int NEVENTS = 0;
  int TOTPARTS = 0;
  int MISMATCHJETE = 0;
  int MISMATCHGENJETE = 0;
  vector<int> dupTrackEvents;
  while(tree_reader.Next()) {

    if(NEVENTS%10000 == 0) cout << "Events Processed: " << NEVENTS << endl;

    // Reconstructed Particles for Clustering
    vector<PseudoJet> particles;
    vector<PseudoJet> particlesGenLoc;
    vector<PseudoJet> particlesGen;

    // Look at ReconstructedJets
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

		    cout << "    Particle " << j << " Energy = " << recoNRG[j] << " Eta = " << local.PseudoRapidity() << " Phi = " << local.Phi() << " Mass = " << recoM[j] << " Pt = " << local.Perp() << endl;

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


    // Look at Generated Jets
    int numGenJets = 0;
    for(unsigned int i=0; i<genType.GetSize(); i++) // Loop over Entries in GeneratedJet
      {
	int numParts = 0;
	TLorentzVector totVec;
	if(genType[i] == 0) // Select Jets
	  {
	    TVector3 jetMom(genMomX[i],genMomY[i],genMomZ[i]);

	    genJetEvsEtaHist->Fill(jetMom.PseudoRapidity(),genNRG[i]);
	    genJetPhiVsEtaHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());

	    numGenJets++;

	    for(unsigned int j=0; j<genType.GetSize(); j++) // For each jet, loop through all entries again
	      {
		if((genType[j] == 1) && (genPDG[j] == genPDG[i])) // Select particles associated with jet
		  {
		    TLorentzVector local;
		    local.SetPxPyPzE(genMomX[j],genMomY[j],genMomZ[j],genNRG[j]);
		    totVec += local;

		    genJetPartEvsEtaHist->Fill(local.PseudoRapidity(),genNRG[j]);
		    genJetPartPhiVsEtaHist->Fill(local.PseudoRapidity(),local.Phi());

		    numParts++;
		  }
	      }

	    numGenJetPartsHist->Fill(numParts);
	    genJetEvsPartESumHist->Fill(genNRG[i],totVec.E());
	    genJetEDiffHist->Fill(genNRG[i]-totVec.E());

	    if(TMath::Abs(totVec.E() - genNRG[i]) > 0.00001)
	      {
		genJetEvsEtaBadHist->Fill(jetMom.PseudoRapidity(),genNRG[i]);
		genJetPhiVsEtaBadHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());
	      }
	  }
      }

    numGenJetsEventHist->Fill(numGenJets);


    // Print Reconstructed and Generated Particle Information
    for(unsigned int i=0; i<recoPartMomX.GetSize(); i++)
      {
	TVector3 partMom(recoPartMomX[i],recoPartMomY[i],recoPartMomZ[i]);
	double E = std::sqrt(recoPartMomX[i]*recoPartMomX[i]+recoPartMomY[i]*recoPartMomY[i]+recoPartMomZ[i]*recoPartMomZ[i]+recoPartM[i]*recoPartM[i]);

	recoPartMassHist->Fill(recoPartM[i]);
	
	cout << "Reco Particle " << i << " Energy = " << E << " Eta = " << partMom.PseudoRapidity() << " Phi = " << partMom.Phi() << " Mass = " << recoPartM[i] << " Pt = " << partMom.Perp() << endl;

	// Fill Local Jet Using Same Cuts as ReconstructedJet
	if(partMom.Perp() > 0.2 && partMom.Perp() < 100.0)
	  {
	    particles.push_back( PseudoJet(partMom.Px(),partMom.Py(),partMom.Pz(),E) );
	  }
      }
    cout << endl;
    for(unsigned int i=0; i<mcMomX.GetSize(); i++)
      {
	if(mcGenStat[i] == 1)
	  {
	    TVector3 partMom(mcMomX[i],mcMomY[i],mcMomZ[i]);
	    double E = std::sqrt(mcMomX[i]*mcMomX[i] + mcMomY[i]*mcMomY[i] + mcMomZ[i]*mcMomZ[i] + mcM[i]*mcM[i]);
	    int id = std::abs(pdg[i]);
	    
	    cout << "Particle " << i << " Energy = " << E << " Eta = " << partMom.PseudoRapidity() << " Phi = " << partMom.Phi() << " PDG = " << pdg[i] << endl;

	    // MC Jet Particles for GeneratedJets / Local Comparison
	    if(partMom.Perp() > 0.2 && partMom.Perp() < 100.0)
	      {
		particlesGenLoc.push_back( PseudoJet(partMom.Px(),partMom.Py(),partMom.Pz(),E) );
	      }

	    // MC Jet Particles for Reco/Truth Comparison: No min pT cut and only Charged Particles
	    if(id==11 || id==13 || id==211 || id==321 || id==2212)
	      {
		particlesGen.push_back( PseudoJet(partMom.Px(),partMom.Py(),partMom.Pz(),E) );
	      }
	  }
      }
    cout << endl;


    // Find Jets from Reconstructed Particles
    // Define Jet
    JetDefinition jet_def_akt(antikt_algorithm,1.0);

    // Cluster in Beam Frame
    ClusterSequence cs_akt_reco(particles, jet_def_akt);
    ClusterSequence cs_akt_gen_loc(particlesGenLoc, jet_def_akt);
    ClusterSequence cs_akt_gen(particlesGen, jet_def_akt);

    // Set Min Jet Pt
    //double ptmin = 1.0;

    // Get Jets
    vector<PseudoJet> jets_akt_reco = sorted_by_pt(cs_akt_reco.inclusive_jets(1.0));
    vector<PseudoJet> jets_akt_gen_loc = sorted_by_pt(cs_akt_gen_loc.inclusive_jets(1.0));
    vector<PseudoJet> jets_akt_gen = sorted_by_pt(cs_akt_gen.inclusive_jets(1.0));

    // Compare to Jets from EICrecon
    numLocalVsRecoJetsHist->Fill(numJets,jets_akt_reco.size());
    numLocalVsGenJetsHist->Fill(numGenJets,jets_akt_gen_loc.size());

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

	    // Match reco Jet to Local Jet
	    double minDist = 10000.0;
	    int matchIndex = -1;
	    for(unsigned int jn=0; jn<jets_akt_reco.size(); jn++)
	      {
		double dEta = jetMom.PseudoRapidity() - jets_akt_reco[jn].eta();
		double dPhi = TVector2::Phi_mpi_pi(jetMom.Phi() - jets_akt_reco[jn].phi());
		double dR = std::sqrt(dEta*dEta + dPhi*dPhi);

		if(dR < minDist)
		  {
		    minDist = dR;
		    matchIndex = jn;
		  }
	      }

	    // Match reco jet to local Gen jet
	    double minDistGen = 10000.0;
	    int matchIndexGen = -1;
	    for(unsigned int jn=0; jn<jets_akt_gen.size(); jn++)
	      {
		double dEta = jetMom.PseudoRapidity() - jets_akt_gen[jn].eta();
		double dPhi = TVector2::Phi_mpi_pi(jetMom.Phi() - jets_akt_gen[jn].phi());
		double dR = std::sqrt(dEta*dEta + dPhi*dPhi);

		if(dR < minDistGen)
		  {
		    minDistGen = dR;
		    matchIndexGen = jn;
		  }
	      }

	    // Tag Local Reco Jets With Duplicate Tracks
	    bool dupTrack = false;

	    // Look at Matched Reco / Local jets
	    if(matchIndex > -1)
	      {
		localVsRecoJetPtHist->Fill(jetMom.Perp(),jets_akt_reco[matchIndex].pt());
		localVsRecoJetEHist->Fill(recoNRG[i],jets_akt_reco[matchIndex].e());
		localVsRecoJetEtaHist->Fill(jetMom.PseudoRapidity(),jets_akt_reco[matchIndex].eta());
		localVsRecoJetPhiHist->Fill(TVector2::Phi_mpi_pi(jetMom.Phi()),TVector2::Phi_mpi_pi(jets_akt_reco[matchIndex].phi()));

		// Compare Constituents
		vector<PseudoJet> cons = jets_akt_reco[matchIndex].constituents();

		int numLocalJetPart = 0;
		for(unsigned int jc=0; jc<cons.size(); jc++)
		  {
		    if(cons[jc].pt() > 0.2 && cons[jc].pt() < 100.0)
		      numLocalJetPart++;

		    if(jc < cons.size() -1)
		      {
			for(unsigned int jc1 = jc+1; jc1<cons.size(); jc1++)
			  {
			    double dEta = cons[jc].eta() - cons[jc1].eta();
			    double dPhi = TVector2::Phi_mpi_pi(cons[jc].phi() - cons[jc1].phi());
			    double dR = std::sqrt(dEta*dEta + dPhi*dPhi);

			    localPartCombDRHist->Fill(dR);

			    if(dR < 0.05) 
			      {
				dupTrack = true;
				dupTrackEvents.push_back(NEVENTS);
			      }
			  }
		      }
		  }

		numLocalVsRecoJetPartsHist->Fill(numParts,numLocalJetPart);

		if(std::abs(recoNRG[i] - jets_akt_reco[matchIndex].e()) > 10e-6) 
		  {
		    //cout << "Event " << NEVENTS << " Jet " << i << " has Mismatched Energy = " << recoNRG[i] - jets_akt_all_beam[matchIndex].e() << endl;
		    MISMATCHJETE++;
		  }
	      }

	    // Look at matched local Gen vs Reco jets
	    localGenRecoDetlaRHist->Fill(minDistGen);
	    if(matchIndexGen > -1 && minDistGen < 0.275)
	      {
		localGenVsRecoJetEHist->Fill(recoNRG[i],jets_akt_gen[matchIndexGen].e());
		localGenVsRecoJetEtaHist->Fill(jetMom.PseudoRapidity(),jets_akt_gen[matchIndexGen].eta());
		localGenVsRecoJetPhiHist->Fill(TVector2::Phi_mpi_pi(jetMom.Phi()),TVector2::Phi_mpi_pi(jets_akt_gen[matchIndexGen].phi()));

		localGenRecoJetEDiffHist->Fill((recoNRG[i]-jets_akt_gen[matchIndexGen].e())/jets_akt_gen[matchIndexGen].e());
		localGenRecoJetEDiffVsRecoEHist->Fill(recoNRG[i],(recoNRG[i]-jets_akt_gen[matchIndexGen].e())/jets_akt_gen[matchIndexGen].e());
		if(dupTrack) localGenRecoJetEDiffDupTrackHist->Fill((recoNRG[i]-jets_akt_gen[matchIndexGen].e())/jets_akt_gen[matchIndexGen].e());
		if(!dupTrack) localGenRecoJetEDiffNoDupTrackHist->Fill((recoNRG[i]-jets_akt_gen[matchIndexGen].e())/jets_akt_gen[matchIndexGen].e());

		localGenRecoJetEDiffVsDeltaRHist->Fill(minDistGen,(recoNRG[i]-jets_akt_gen[matchIndexGen].e())/jets_akt_gen[matchIndexGen].e());

		cout << "Gen Jet (Reco = " << recoPDG[i] << ") Energy = " << jets_akt_gen[matchIndexGen].e() << " Eta = " << jets_akt_gen[matchIndexGen].eta() << " Phi = " << TVector2::Phi_mpi_pi(jets_akt_gen[matchIndexGen].phi()) << " Delta R = " << minDistGen << " Delta E/E = " << (recoNRG[i]-jets_akt_gen[matchIndexGen].e())/jets_akt_gen[matchIndexGen].e() << endl;

		// Compare Constituents
		vector<PseudoJet> cons = jets_akt_gen[matchIndexGen].constituents();
		for(unsigned int jc=0; jc<cons.size(); jc++)
		  {
		    cout << "    Particle " << jc << " Energy = " << cons[jc].e() << " Eta = " << cons[jc].eta() << " Phi = " << TVector2::Phi_mpi_pi(cons[jc].phi()) << " Pt = " << cons[jc].pt() << endl;

		    if(jc < cons.size() -1)
		      {
			for(unsigned int jc1 = jc+1; jc1<cons.size(); jc1++)
			  {
			    double dEta = cons[jc].eta() - cons[jc1].eta();
			    double dPhi = TVector2::Phi_mpi_pi(cons[jc].phi() - cons[jc1].phi());
			    double dR = std::sqrt(dEta*dEta + dPhi*dPhi);

			    localGenPartCombDRHist->Fill(dR);
			  }
		      }
		  }
	      }
	    cout << endl;
	  }
      }


    // Loop Over Gen Jets
    for(unsigned int i=0; i<genType.GetSize(); i++)
      {
	int numParts = 0;
	TLorentzVector totVec;
	if(genType[i] == 0) // Select Jets
	  {
	    TVector3 jetMom(genMomX[i],genMomY[i],genMomZ[i]);

	    for(unsigned int j=0; j<genType.GetSize(); j++) // For each jet, loop through all entries again
	      {
		if((genType[j] == 1) && (genPDG[j] == genPDG[i])) // Select particles associated with jet
		  {
		    TLorentzVector local;
		    local.SetPxPyPzE(genMomX[j],genMomY[j],genMomZ[j],genNRG[j]);
		    totVec += local;
		    
		    numParts++;
		  }
	      }

	    // Match Gen Jet to Local Jet
	    double minDist = 10000.0;
	    int matchIndex = -1;
	    for(unsigned int jn=0; jn<jets_akt_gen_loc.size(); jn++)
	      {
		double dEta = jetMom.PseudoRapidity() - jets_akt_gen_loc[jn].eta();
		double dPhi = TVector2::Phi_mpi_pi(jetMom.Phi() - jets_akt_gen_loc[jn].phi());
		double dR = std::sqrt(dEta*dEta + dPhi*dPhi);

		if(dR < minDist)
		  {
		    minDist = dR;
		    matchIndex = jn;
		  }
	      }

	    // Look at Matched Gen / Local jets
	    if(matchIndex > -1)
	      {
		localVsGenJetPtHist->Fill(jetMom.Perp(),jets_akt_gen_loc[matchIndex].pt());
		localVsGenJetEHist->Fill(genNRG[i],jets_akt_gen_loc[matchIndex].e());
		localVsGenJetEtaHist->Fill(jetMom.PseudoRapidity(),jets_akt_gen_loc[matchIndex].eta());
		localVsGenJetPhiHist->Fill(TVector2::Phi_mpi_pi(jetMom.Phi()),TVector2::Phi_mpi_pi(jets_akt_gen_loc[matchIndex].phi()));

		// Compare Constituents
		vector<PseudoJet> cons = jets_akt_gen_loc[matchIndex].constituents();

		int numLocalJetPart = 0;
		for(unsigned int jc=0; jc<cons.size(); jc++)
		  {
		    if(cons[jc].pt() > 0.2 && cons[jc].pt() < 100.0)
		      numLocalJetPart++;
		    /*
		    if(jc < cons.size() -1)
		      {
			for(unsigned int jc1 = jc+1; jc1<cons.size(); jc1++)
			  {
			    double dEta = cons[jc].eta() - cons[jc1].eta();
			    double dPhi = TVector2::Phi_mpi_pi(cons[jc].phi() - cons[jc1].phi());
			    double dR = std::sqrt(dEta*dEta + dPhi*dPhi);

			    localGenPartCombDRHist->Fill(dR);
			  }
		      }
		    */
		  }

		numLocalVsGenJetPartsHist->Fill(numParts,numLocalJetPart);

		if(std::abs(genNRG[i] - jets_akt_gen_loc[matchIndex].e()) > 10e-6) 
		  {
		    //cout << "Event " << NEVENTS << " Jet " << i << " has Mismatched Energy = " << recoNRG[i] - jets_akt_all_beam[matchIndex].e() << endl;
		    MISMATCHGENJETE++;
		  }
	      }
	  }
      }

    NEVENTS++;
  } // End Loop Over Events

  cout << "Total Particles = " << TOTPARTS << endl;
  cout << "Total Reco/Local Jet Pairs with Mismatched Energy = " << MISMATCHJETE << endl;
  cout << "Total Gen/Local Jet Pairs with Mismatched Energy = " << MISMATCHGENJETE << endl;
  //cout << "Events with Duplicate Tracks" << endl;
  for(unsigned int i=0; i<dupTrackEvents.size(); i++)
    {
      //cout << dupTrackEvents[i] << endl;
    }

  //mcEta->Write();
  //mcNRG->Write();

  ofile->Write();
  ofile->Close();
  return 0;
}
