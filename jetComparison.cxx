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

  /*
  auto file = new TFile(infile);
  auto tree = (TTree *)file->Get("events");
  //tree->Print();
  */

  //if(geometry == 0) cout << "Using Arches Geometry" << endl;
  //if(geometry == 1) cout << "Using BryceCanyon Geometry" << endl;
  //cout << infile << endl;

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


  /*
  TH1D *mcGen = new TH1D("mcGen","",100,0.,100.);
  TH1D *mcPt = new TH1D("mcPt","",100,0.,25.);
  TH1D *mcEta = new TH1D("mcEta","",100,-5.,5.);
  TH1D *mcPhi = new TH1D("mcPhi","",100,-TMath::Pi(),TMath::Pi());
  TH1D *mcNRG = new TH1D("mcNRG","",200,0.,100.);
  TH2D *mcPhiVsEta = new TH2D("mcPhiVsEta","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH1D *mcNRGCal = new TH1D("mcNRGCal","",200,0.,100.);
  TH1D *mcNRGElec = new TH1D("mcNRGElec","",200,0.,100.);
  TH1D *mcNRGGamma = new TH1D("mcNRGGamma","",200,0.,100.);
  //TH2D *mcNRGVsEtaElec = new TH2D("mcNRGVsEtaElec","",100,-5.,5.,200,0.,100.);
  TH2D *mcNRGVsEtaGamma = new TH2D("mcNRGVsEtaGamma","",100,-5.,5.,200,0.,100.);

  TH2D *mcPhiVsEtaBeam = new TH2D("mcPhiVsEtaBeam","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH2D *mcPhiVsEtaBeamTest = new TH2D("mcPhiVsEtaBeamTest","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());
  TH2D *mcPhiVsEtaTest = new TH2D("mcPhiVsEtaTest","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  //TH2D *testMCBarrelScatElecEVsEta = new TH2D("testMCBarrelScatElecEVsEta","",100,-5.,5.,200,0.,100.);
  TH2D *mcNRGVsEtaElecHist = new TH2D("mcNRGVsEtaElec","",100,-5.,5.,100,0.,100.);
  TH2D *mcPhiVsEtaElecHist = new TH2D("mcPhiVsEtaElec","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  // Reconstructed Particles
  TTreeReaderArray<float> recoMomX = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
  TTreeReaderArray<float> recoMomY = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
  TTreeReaderArray<float> recoMomZ = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
  TTreeReaderArray<float> recoM = {tree_reader, "ReconstructedChargedParticles.mass"};

  TH1D *recoPt = new TH1D("recoPt","",100,0.,25.);
  TH1D *recoEta = new TH1D("recoEta","",100,-5.,5.);
  TH1D *recoPhi = new TH1D("recoPhi","",100,-TMath::Pi(),TMath::Pi());
  TH1D *recoNRG = new TH1D("recoNRG","",200,0.,100.);
  TH2D *recoPhiVsEta = new TH2D("recoPhiVsEta","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH2D *recoPhiVsEtaBeam = new TH2D("recoPhiVsEtaBeam","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH2D *recoPhiVsEtaBeamTest = new TH2D("recoPhiVsEtaBeamTest","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());
  TH2D *recoPhiVsEtaTest = new TH2D("recoPhiVsEtaTest","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());
  */
  
  /*
  // Jets
  char jetNames[5][20];
  sprintf(jetNames[0],"allMC_beam");
  sprintf(jetNames[1],"chargeMC_beam");
  sprintf(jetNames[2],"reco_beam");
  sprintf(jetNames[3],"recoECal_beam");
  sprintf(jetNames[4],"recoECalTruth_beam");

  TH1D *jetPtHist[5];
  TH1D *jetEHist[5];
  TH2D *jetPtVsEtaHist[5];
  TH2D *jetEVsEtaHist[5];
  TH2D *jetPhiVsEtaHist[5];

  for(int i=0; i<5; i++)
    {
      jetPtHist[i] = new TH1D(Form("jetPt_%s",jetNames[i]),"",100,0,50.);
      jetEHist[i] = new TH1D(Form("jetE_%s",jetNames[i]),"",100,0,100.);
      jetPtVsEtaHist[i] = new TH2D(Form("jetPtVsEta_%s",jetNames[i]),"",100,-5.,5.,100,0.,50.);
      jetEVsEtaHist[i] = new TH2D(Form("jetEVsEta_%s",jetNames[i]),"",100,-5.,5.,100,0.,100.);
      jetPhiVsEtaHist[i] = new TH2D(Form("jetPhiVsEta_%s",jetNames[i]),"",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());
    }

  char jetComp[4][30];
  sprintf(jetComp[0],"allMCReco_loopReco");
  sprintf(jetComp[1],"chargeMCReco_loopReco");
  sprintf(jetComp[2],"allMCRecoECal_loopReco");
  sprintf(jetComp[3],"allMCRecoECalTruth_loopReco");

  TH1D *jetDeltaRHist[4];
  TH2D *jetMissingPtVsEtaHist[4];
  TH2D *jetMissingEVsEtaHist[4];
  TH2D *jetMissingPhiVsEtaHist[4];
  
  TH2D *jetRecoVsPartPtHist[4];
  TH2D *jetRecoVsPartEHist[4];
  
  TH2D *jetPtResVsPtEtaBin0Hist[4];
  TH2D *jetEResVsEEtaBin0Hist[4];
  TH2D *jetPtResVsPtEtaBin1Hist[4];
  TH2D *jetEResVsEEtaBin1Hist[4];
  TH2D *jetPtResVsPtEtaBin2Hist[4];
  TH2D *jetEResVsEEtaBin2Hist[4];
  
  for(int i=0; i<4; i++)
    {
      jetDeltaRHist[i] = new TH1D(Form("jetDeltaR_%s",jetComp[i]),"",1000,0.,10.);

      jetMissingPtVsEtaHist[i] = new TH2D(Form("jetMissingPtVsEta_%s",jetComp[i]),"",100,-5.,5.,100,0.,50.);
      jetMissingEVsEtaHist[i] = new TH2D(Form("jetMissingEVsEta_%s",jetComp[i]),"",100,-5.,5.,100,0.,100.); 
      jetMissingPhiVsEtaHist[i] = new TH2D(Form("jetMissingPhiVsEta_%s",jetComp[i]),"",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());
      
      jetRecoVsPartPtHist[i] = new TH2D(Form("jetRecoVsPartPt_%s",jetComp[i]),"",100,0.,50.,100,0.,50.);
      jetRecoVsPartEHist[i] = new TH2D(Form("jetRecoVsPartE_%s",jetComp[i]),"",100,0.,100.,100,0.,100.);
      
      jetPtResVsPtEtaBin0Hist[i] = new TH2D(Form("jetPtResVsPtEtaBin0_%s",jetComp[i]),"",100,0.,50.,10000,-10.,10.);
      jetEResVsEEtaBin0Hist[i] = new TH2D(Form("jetEResVsEEtaBin0_%s",jetComp[i]),"",100,0.,100.,10000,-10.,10.);
      jetPtResVsPtEtaBin1Hist[i] = new TH2D(Form("jetPtResVsPtEtaBin1_%s",jetComp[i]),"",100,0.,50.,10000,-10.,10.);
      jetEResVsEEtaBin1Hist[i] = new TH2D(Form("jetEResVsEEtaBin1_%s",jetComp[i]),"",100,0.,100.,10000,-10.,10.);
      jetPtResVsPtEtaBin2Hist[i] = new TH2D(Form("jetPtResVsPtEtaBin2_%s",jetComp[i]),"",100,0.,50.,10000,-10.,10.);
      jetEResVsEEtaBin2Hist[i] = new TH2D(Form("jetEResVsEEtaBin2_%s",jetComp[i]),"",100,0.,100.,10000,-10.,10.);
    }
  */
  

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
	    
	    


    //cout << "Number of Particles: " << mcMomZ.GetSize() << endl;

    //numHcalEndNHits->Fill(hcalEndNHitE.GetSize());
    //numClust->Fill(hcalEndNClustE.GetSize());

    //vector<PseudoJet> particles;
    //vector<PseudoJet> particlesCharged;
    //vector<PseudoJet> particlesReco;
    //vector<PseudoJet> particlesECal;
    //vector<PseudoJet> particlesRecoECal;
    //vector<PseudoJet> particlesECalTruth;
    //vector<PseudoJet> particlesRecoECalTruth;

    /*
    // Define MC Beam Kinematics
    TLorentzVector hBeamMC;
    if(pdg[0] != 2212 || mcGenStat[0] != 4)
      {
	cout << "PROBLEM GETTING MC HADRON BEAM" << endl;
	break;
      }
    */

    //hBeamMC.SetPxPyPzE(mcMomX[0],mcMomY[0],mcMomZ[0],std::sqrt(mcMomX[0]*mcMomX[0] + mcMomY[0]*mcMomY[0] + mcMomZ[0]*mcMomZ[0] + mcM[0]*mcM[0]));
    //hBeamMC.SetPxPyPzE(-6.8742839,0.0,274.91407,275.0);	

    /*
    // Identify Electron in Barrel
    Double_t scatElecEta = -99.;
    Double_t scatElecPhi = -99.;
    Double_t scatElecE = -99.;
    //bool firstElec = true;
    Double_t maxE = 0.0;

    // Loop Over MC Particles
    for(unsigned int i=0; i<mcMomX.GetSize(); i++)
      {
	mcGen->Fill(mcGenStat[i]);

	if(i < 10)
	  {
	    //cout << i << " " << pdg[i] << " " << mcGenStat[i] << " " << mcMomX[i] << " " << mcMomY[i] << " " << mcMomZ[i] << " " << std::sqrt(mcMomX[i]*mcMomX[i] + mcMomY[i]*mcMomY[i] + mcMomZ[i]*mcMomZ[i] + mcM[i]*mcM[i]) << endl;
	  }

	if(mcGenStat[i] != 1) continue; // Only look at thrown particles

	TVector3 mom(mcMomX[i],mcMomY[i],mcMomZ[i]);
	double E = std::sqrt(mcMomX[i]*mcMomX[i] + mcMomY[i]*mcMomY[i] + mcMomZ[i]*mcMomZ[i] + mcM[i]*mcM[i]);
	mcPt->Fill(mom.Perp());
	mcEta->Fill(mom.PseudoRapidity());
	mcPhi->Fill(mom.Phi());
	mcNRG->Fill(E);
	mcPhiVsEta->Fill(mom.PseudoRapidity(),mom.Phi());

	//cout << i << " " << pdg[i] << " " << mom.PseudoRapidity() << " " << E << endl;
	//if(mom.PseudoRapidity() > -1.8 && mom.PseudoRapidity() < 1.4)
	//if(mom.PseudoRapidity() > -1.8 && mom.PseudoRapidity() < 1.4 && firstElec)
	if(1)
	  {
	    if(pdg[i] == 11)
	      {
		if(E > maxE) // E > maxE
		  {
		    scatElecEta = mom.PseudoRapidity();
		    scatElecPhi = mom.Phi();
		    scatElecE = E;
		    //firstElec = false;
		    maxE = E;
		  }
	      }
	  }

	if(mom.PseudoRapidity() > -1.8 && mom.PseudoRapidity() < 1.4)
	  {
	    if(TMath::Abs(pdg[i]) == 11 || TMath::Abs(pdg[i]) == 22)
	      mcNRGCal->Fill(E);
	    
	    if(TMath::Abs(pdg[i]) == 11)
	      {
		mcNRGElec->Fill(E);
		//mcNRGVsEtaElec->Fill(mom.PseudoRapidity(),E);
	      }

	    if(TMath::Abs(pdg[i]) == 22)
	      {
		mcNRGGamma->Fill(E);
		mcNRGVsEtaGamma->Fill(mom.PseudoRapidity(),E);
	      }
	  }
	    
	// Calculate Kinematics W.R.T. Hadron Beam Line
	TLorentzVector partWRThBeam;
	partWRThBeam.SetPxPyPzE(mcMomX[i],mcMomY[i],mcMomZ[i],E);
        
	// Calc Track Pt wrt HBeam
	TVector3 longUnitB = hBeamMC.Vect().Unit();
	TVector3 normUnitB = TVector3(0,0,1).Cross(longUnitB).Unit();
	TVector3 sideUnitB = longUnitB.Cross(normUnitB);
	TVector3 momB = partWRThBeam.Vect();
	TVector3 localMomB(momB.Dot(sideUnitB),momB.Dot(normUnitB),momB.Dot(longUnitB));
	
	mcPhiVsEtaBeam->Fill(localMomB.PseudoRapidity(),localMomB.Phi());

	if(mom.PseudoRapidity() >= 1.4)
	  {
	    mcPhiVsEtaBeamTest->Fill(localMomB.PseudoRapidity(),localMomB.Phi());
	  }
	if(mom.PseudoRapidity() < 1.4)
	  mcPhiVsEtaTest->Fill(mom.PseudoRapidity(),mom.Phi());
	
	// Fill FastJet
	if(localMomB.PseudoRapidity() > -5.0 && localMomB.PseudoRapidity() < 5.0)
	  {
	    int apdg = TMath::Abs(pdg[i]);

	    if(mom.PseudoRapidity() < 1.4)
	      {
		fastjet::PseudoJet p(mom.Px(),mom.Py(),mom.Pz(),E);
		p.set_user_index(i);
		if(apdg != 11) particles.push_back(p);
		//particles.push_back(p);

		if(apdg == 211 || apdg == 321 || apdg == 2212)
		  particlesCharged.push_back(p);
	      }

	    if(mom.PseudoRapidity() >= 1.4 && localMomB.PseudoRapidity() >= 1.4)
	      {
		fastjet::PseudoJet p(localMomB.Px(),localMomB.Py(),localMomB.Pz(),E);
		p.set_user_index(i);
		if(apdg != 11) particles.push_back(p);
		//particles.push_back(p);

		if(apdg == 211 || apdg == 321 || apdg == 2212)
		  particlesCharged.push_back(p);
	      }
	  }
      }

    mcNRGVsEtaElecHist->Fill(scatElecEta,scatElecE);
    mcPhiVsEtaElecHist->Fill(scatElecEta,scatElecPhi);
    */

    /*
    // Loop Over Reconstructed Charged Particles
    for(unsigned int i=0; i<recoMomX.GetSize(); i++)
      {
	TVector3 mom(recoMomX[i],recoMomY[i],recoMomZ[i]);
	double E = std::sqrt(recoMomX[i]*recoMomX[i] + recoMomY[i]*recoMomY[i] + recoMomZ[i]*recoMomZ[i] + recoM[i]*recoM[i]);
	recoPt->Fill(mom.Perp());
	recoEta->Fill(mom.PseudoRapidity());
	recoPhi->Fill(mom.Phi());
	recoNRG->Fill(E);
	recoPhiVsEta->Fill(mom.PseudoRapidity(),mom.Phi());

	// Calculate Kinematics W.R.T. Hadron Beam Line
	TLorentzVector partWRThBeam;
	partWRThBeam.SetPxPyPzE(recoMomX[i],recoMomY[i],recoMomZ[i],E);
        
	// Calc Track Pt wrt HBeam
	TVector3 longUnitB = hBeamMC.Vect().Unit();
	TVector3 normUnitB = TVector3(0,0,1).Cross(longUnitB).Unit();
	TVector3 sideUnitB = longUnitB.Cross(normUnitB);
	TVector3 momB = partWRThBeam.Vect();
	TVector3 localMomB(momB.Dot(sideUnitB),momB.Dot(normUnitB),momB.Dot(longUnitB));
	
	recoPhiVsEtaBeam->Fill(localMomB.PseudoRapidity(),localMomB.Phi());

	if(mom.PseudoRapidity() >= 1.4)
	  recoPhiVsEtaBeamTest->Fill(localMomB.PseudoRapidity(),localMomB.Phi());

	if(mom.PseudoRapidity() < 1.4)
	  recoPhiVsEtaTest->Fill(mom.PseudoRapidity(),mom.Phi());

	// Fill FastJet
	if(localMomB.PseudoRapidity() > -3.2 && localMomB.PseudoRapidity() < 3.2)
	  {
	    if(mom.PseudoRapidity() < 1.4)
	      {
		fastjet::PseudoJet p(mom.Px(),mom.Py(),mom.Pz(),E);
		p.set_user_index(i);
		if(recoM[i] > 0.001) particlesReco.push_back(p);
		if(recoM[i] > 0.001) particlesRecoECal.push_back(p);
		if(recoM[i] > 0.001) particlesRecoECalTruth.push_back(p);
		//particlesReco.push_back(p);
		//particlesRecoECal.push_back(p);
	      }

	    if(mom.PseudoRapidity() >= 1.4 && localMomB.PseudoRapidity() >= 1.4)
	      {
		fastjet::PseudoJet p(localMomB.Px(),localMomB.Py(),localMomB.Pz(),E);
		p.set_user_index(i);
		if(recoM[i] > 0.001) particlesReco.push_back(p);
		if(recoM[i] > 0.001) particlesRecoECal.push_back(p);
		if(recoM[i] > 0.001) particlesRecoECalTruth.push_back(p);
		//particlesReco.push_back(p);
		//particlesRecoECal.push_back(p);
	      }
	  }

      }
    */

    /*
    // Loop Over Merged nECal Clusters
    //cout << "Event " << NEVENTS << endl;
    for(unsigned int i=0; i<nECalPosX.GetSize(); i++)
      {
	TVector3 pos(nECalPosX[i],nECalPosY[i],nECalPosZ[i]);

	nECalPosYvsXHist->Fill(nECalPosX[i],nECalPosY[i]);
	nECalPosZHist->Fill(TMath::Abs(nECalPosZ[i]));

	if(nECalPosX[i] > -100.0 && nECalPosX[i] < 100.0 && nECalPosY[i] > -100.0 && nECalPosY[i] < 100.0) // Kill Clusters in Gap
	  continue;

	Double_t px = nECalE[i]*TMath::Sin(pos.Theta())*TMath::Cos(pos.Phi());
	Double_t py = nECalE[i]*TMath::Sin(pos.Theta())*TMath::Sin(pos.Phi());
	Double_t pz = nECalE[i]*TMath::Cos(pos.Theta());

	TVector3 mom(px,py,pz);

	nECalNRGHist->Fill(nECalE[i]);
	nECalPhiVsEtaHist->Fill(pos.PseudoRapidity(),TVector2::Phi_mpi_pi(pos.Phi()));

	int clusterPID = pdg[nECalSimID[i]];

	if(clusterPID != 11)
	  {
	    nECalNRGNoElecHist->Fill(nECalE[i]);
	    nECalPhiVsEtaNoElecHist->Fill(pos.PseudoRapidity(),TVector2::Phi_mpi_pi(pos.Phi()));
	  }
	if(clusterPID == 11)
	  nECalNRGVsEtaElecHist->Fill(pos.PseudoRapidity(),nECalE[i]);

	// Calc Track Pt wrt HBeam
	TVector3 longUnitB = hBeamMC.Vect().Unit();
	TVector3 normUnitB = TVector3(0,0,1).Cross(longUnitB).Unit();
	TVector3 sideUnitB = longUnitB.Cross(normUnitB);
	//TVector3 momB = partWRThBeam.Vect();
	TVector3 localMomB(mom.Dot(sideUnitB),mom.Dot(normUnitB),mom.Dot(longUnitB));
	TVector3 localPosB(pos.Dot(sideUnitB),pos.Dot(normUnitB),pos.Dot(longUnitB));

	nECalPhiVsEtaBeamHist->Fill(localMomB.PseudoRapidity(),TVector2::Phi_mpi_pi(localMomB.Phi()));
	nECalPosYvsXBeamHist->Fill(localPosB.Px(),localPosB.Py());
	if(clusterPID != 11)
	  nECalPhiVsEtaNoElecBeamHist->Fill(localMomB.PseudoRapidity(),TVector2::Phi_mpi_pi(localMomB.Phi()));

	//Fill FastJet
	if(localMomB.PseudoRapidity() > -4.0 && localMomB.PseudoRapidity() < -1.8)
	  {
	    fastjet::PseudoJet p(mom.Px(),mom.Py(),mom.Pz(),nECalE[i]);
	    p.set_user_index(i);
	    if(clusterPID != 11) // clusterPID != 11
	      {
		particlesECal.push_back(p);
		particlesRecoECal.push_back(p);
	      }
	  }

	//cout << "Cluster " << i << " " << nECalPosX[i] << " " << nECalPosY[i] << " " << nECalPosZ[i] << " " << nECalE[i] << " " << pidBegin[i] << " " << pidEnd[i] << " " << pos.PseudoRapidity() << " " << pos.Phi() << endl;
      }
    */

    /*
    // Loop Over Truth nECal Clusters
    for(unsigned int i=0; i<nECalTruthPosX.GetSize(); i++)
      {
	TVector3 pos(nECalTruthPosX[i],nECalTruthPosY[i],nECalTruthPosZ[i]);

	nECalTruthPosYvsXHist->Fill(nECalTruthPosX[i],nECalTruthPosY[i]);
	nECalTruthPosZHist->Fill(TMath::Abs(nECalTruthPosZ[i]));

	if(nECalTruthPosX[i] > -100.0 && nECalTruthPosX[i] < 100.0 && nECalTruthPosY[i] > -100.0 && nECalTruthPosY[i] < 100.0) // Kill Clusters in Gap
	  continue;

	Double_t px = nECalTruthE[i]*TMath::Sin(pos.Theta())*TMath::Cos(pos.Phi());
	Double_t py = nECalTruthE[i]*TMath::Sin(pos.Theta())*TMath::Sin(pos.Phi());
	Double_t pz = nECalTruthE[i]*TMath::Cos(pos.Theta());

	TVector3 mom(px,py,pz);

	nECalTruthNRGHist->Fill(nECalTruthE[i]);
	nECalTruthPhiVsEtaHist->Fill(pos.PseudoRapidity(),TVector2::Phi_mpi_pi(pos.Phi()));

	int clusterPID = pdg[nECalTruthSimID[i]];

	// Calc Track Pt wrt HBeam
	TVector3 longUnitB = hBeamMC.Vect().Unit();
	TVector3 normUnitB = TVector3(0,0,1).Cross(longUnitB).Unit();
	TVector3 sideUnitB = longUnitB.Cross(normUnitB);
	//TVector3 momB = partWRThBeam.Vect();
	TVector3 localMomB(mom.Dot(sideUnitB),mom.Dot(normUnitB),mom.Dot(longUnitB));
	TVector3 localPosB(pos.Dot(sideUnitB),pos.Dot(normUnitB),pos.Dot(longUnitB));

	nECalTruthPhiVsEtaBeamHist->Fill(localMomB.PseudoRapidity(),TVector2::Phi_mpi_pi(localMomB.Phi()));
	nECalTruthPosYvsXBeamHist->Fill(localPosB.Px(),localPosB.Py());

	//Fill FastJet
	if(mom.PseudoRapidity() > -4.0 && mom.PseudoRapidity() < -1.8)
	  {
	    fastjet::PseudoJet p(localMomB.Px(),localMomB.Py(),localMomB.Pz(),nECalTruthE[i]);
	    p.set_user_index(i);
	    if(clusterPID != 11) // clusterPID != 11
	      {
		particlesECalTruth.push_back(p);
		particlesRecoECalTruth.push_back(p);
	      }
	  }
      }

    for(unsigned int i=0; i<nECalSimID.GetSize(); i++)
      {
	TVector3 mom(mcMomX[nECalSimID[i]],mcMomY[nECalSimID[i]],mcMomZ[nECalSimID[i]]);
	double E = std::sqrt(mcMomX[nECalSimID[i]]*mcMomX[nECalSimID[i]] + mcMomY[nECalSimID[i]]*mcMomY[nECalSimID[i]] + mcMomZ[nECalSimID[i]]*mcMomZ[nECalSimID[i]] + mcM[nECalSimID[i]]*mcM[nECalSimID[i]]);

	if(pdg[nECalSimID[i]] == 11)
	  {
	    testElecERes->Fill((nECalE[i]-E)/E);
	  }
	if(pdg[nECalSimID[i]] == 22)
	  {
	    testGammaERes->Fill((nECalE[i]-E)/E);
	  }

	//cout << i << " " << simID[i] << " " << recID[i] << endl;
	//cout << pdg[simID[i]] << " " << E << " " << mom.PseudoRapidity() << " " << mom.Phi() << endl;
      }
    */

    /*
    // Loop bECal Clusters
    //cout << "Event " << NEVENTS << endl;
    int numElecClusters = 0;
    Double_t totElecClusterE = 0;
    for(unsigned int i=0; i<bECalPosX.GetSize(); i++)
      {
	TVector3 pos(bECalPosX[i],bECalPosY[i],bECalPosZ[i]);

	bECalPosYvsXHist->Fill(bECalPosX[i],bECalPosY[i]);
	bECalPosZHist->Fill(bECalPosZ[i]);

	Double_t px = bECalE[i]*TMath::Sin(pos.Theta())*TMath::Cos(pos.Phi());
	Double_t py = bECalE[i]*TMath::Sin(pos.Theta())*TMath::Sin(pos.Phi());
	Double_t pz = bECalE[i]*TMath::Cos(pos.Theta());

	TVector3 mom(px,py,pz);

	bECalNRGHist->Fill(bECalE[i]);
	bECalPhiVsEtaHist->Fill(pos.PseudoRapidity(),TVector2::Phi_mpi_pi(pos.Phi()));

	// Find Electron Cluster
	bool elecFlag = false;
	Double_t elecDeltaEta = mom.PseudoRapidity() - scatElecEta;
	Double_t elecDeltaPhi = TVector2::Phi_mpi_pi(mom.Phi() - scatElecPhi);
	Double_t elecDeltaR = TMath::Sqrt(elecDeltaEta*elecDeltaEta + elecDeltaPhi*elecDeltaPhi);
	if(elecDeltaR < 0.1) elecFlag = true;

	if(mom.PseudoRapidity() > -1.8 && mom.PseudoRapidity() < 1.4)
	  {
	    if(!elecFlag)
	      {
		bECalNRGNoElecHist->Fill(bECalE[i]);
		bECalNRGVsEtaNoElecHist->Fill(mom.PseudoRapidity(),bECalE[i]);
	      }
	    if(elecFlag)
	      {
		bECalNRGVsEtaElecHist->Fill(mom.PseudoRapidity(),bECalE[i]);
		bECalPhiVsEtaElecHist->Fill(mom.PseudoRapidity(),mom.Phi());
		bECalNRGVsMCNRGElecHist->Fill(scatElecE,bECalE[i]);
		numElecClusters++;
		totElecClusterE += bECalE[i];
	      }
	  }

	//int clusterPID = pdg[nECalSimID[i]];
	

	// Calc Track Pt wrt HBeam
	TVector3 longUnitB = hBeamMC.Vect().Unit();
	TVector3 normUnitB = TVector3(0,0,1).Cross(longUnitB).Unit();
	TVector3 sideUnitB = longUnitB.Cross(normUnitB);
	//TVector3 momB = partWRThBeam.Vect();
	TVector3 localMomB(mom.Dot(sideUnitB),mom.Dot(normUnitB),mom.Dot(longUnitB));
	TVector3 localPosB(pos.Dot(sideUnitB),pos.Dot(normUnitB),pos.Dot(longUnitB));

	bECalPhiVsEtaBeamHist->Fill(localMomB.PseudoRapidity(),TVector2::Phi_mpi_pi(localMomB.Phi()));
	bECalPosYvsXBeamHist->Fill(localPosB.Px(),localPosB.Py());
	//if(clusterPID != 11)
	//nECalPhiVsEtaNoElecBeamHist->Fill(localMomB.PseudoRapidity(),TVector2::Phi_mpi_pi(localMomB.Phi()));

	//Fill FastJet
	if(mom.PseudoRapidity() > -1.8 && mom.PseudoRapidity() < 1.4)
	  {
	    fastjet::PseudoJet p(mom.Px(),mom.Py(),mom.Pz(),bECalE[i]);
	    p.set_user_index(i);
	    if(!elecFlag) 
	      {
		particlesECal.push_back(p);
		particlesRecoECal.push_back(p);
	      }
	  }
      }

    bECalTotClusterEVsNumHist->Fill(numElecClusters,totElecClusterE);

    // Loop Over Truth Barrel Clusters
    for(unsigned int i=0; i<bECalTruthPosX.GetSize(); i++)
      {
	TVector3 pos(bECalTruthPosX[i],bECalTruthPosY[i],bECalTruthPosZ[i]);

	bECalTruthPosYvsXHist->Fill(bECalTruthPosX[i],bECalTruthPosY[i]);
	bECalTruthPosZHist->Fill(TMath::Abs(bECalTruthPosZ[i]));

	Double_t px = bECalTruthE[i]*TMath::Sin(pos.Theta())*TMath::Cos(pos.Phi());
	Double_t py = bECalTruthE[i]*TMath::Sin(pos.Theta())*TMath::Sin(pos.Phi());
	Double_t pz = bECalTruthE[i]*TMath::Cos(pos.Theta());

	TVector3 mom(px,py,pz);

	bECalTruthNRGHist->Fill(bECalTruthE[i]);
	bECalTruthPhiVsEtaHist->Fill(pos.PseudoRapidity(),TVector2::Phi_mpi_pi(pos.Phi()));

	// Find Electron Cluster
	bool elecFlag = false;
	Double_t elecDeltaEta = mom.PseudoRapidity() - scatElecEta;
	Double_t elecDeltaPhi = TVector2::Phi_mpi_pi(mom.Phi() - scatElecPhi);
	Double_t elecDeltaR = TMath::Sqrt(elecDeltaEta*elecDeltaEta + elecDeltaPhi*elecDeltaPhi);
	if(elecDeltaR < 0.1) elecFlag = true;

	// Calc Track Pt wrt HBeam
	TVector3 longUnitB = hBeamMC.Vect().Unit();
	TVector3 normUnitB = TVector3(0,0,1).Cross(longUnitB).Unit();
	TVector3 sideUnitB = longUnitB.Cross(normUnitB);
	//TVector3 momB = partWRThBeam.Vect();
	TVector3 localMomB(mom.Dot(sideUnitB),mom.Dot(normUnitB),mom.Dot(longUnitB));
	TVector3 localPosB(pos.Dot(sideUnitB),pos.Dot(normUnitB),pos.Dot(longUnitB));

	bECalTruthPhiVsEtaBeamHist->Fill(localMomB.PseudoRapidity(),TVector2::Phi_mpi_pi(localMomB.Phi()));
	bECalTruthPosYvsXBeamHist->Fill(localPosB.Px(),localPosB.Py());

	//Fill FastJet
	if(mom.PseudoRapidity() > -1.8 && mom.PseudoRapidity() < 1.4)
	  {
	    fastjet::PseudoJet p(mom.Px(),mom.Py(),mom.Pz(),bECalTruthE[i]);
	    p.set_user_index(i);
	    if(!elecFlag) 
	      {
		particlesECalTruth.push_back(p);
		particlesRecoECalTruth.push_back(p);
	      }
	  }
      }

    // Loop Over Merged pECal Clusters
    for(unsigned int i=0; i<pECalPosX.GetSize(); i++)
      {
	TVector3 pos(pECalPosX[i],pECalPosY[i],pECalPosZ[i]);

	pECalPosYvsXHist->Fill(pECalPosX[i],pECalPosY[i]);
	pECalPosZHist->Fill(TMath::Abs(pECalPosZ[i]));

	if(pECalPosX[i] > -380.0 && pECalPosX[i] < 180.0 && pECalPosY[i] > -280.0 && pECalPosY[i] < 280.0) // Kill Clusters in Gap
	  continue;

	Double_t px = pECalE[i]*TMath::Sin(pos.Theta())*TMath::Cos(pos.Phi());
	Double_t py = pECalE[i]*TMath::Sin(pos.Theta())*TMath::Sin(pos.Phi());
	Double_t pz = pECalE[i]*TMath::Cos(pos.Theta());

	//posYvsXHist->Fill(px,py);
	//posZHist->Fill(pz);

	TVector3 mom(px,py,pz);

	pECalNRGHist->Fill(pECalE[i]);
	pECalPhiVsEtaHist->Fill(pos.PseudoRapidity(),TVector2::Phi_mpi_pi(pos.Phi()));

	int clusterPID = pdg[pECalSimID[i]];

	if(clusterPID != 11)
	  {
	    pECalNRGNoElecHist->Fill(pECalE[i]);
	    pECalPhiVsEtaNoElecHist->Fill(pos.PseudoRapidity(),TVector2::Phi_mpi_pi(pos.Phi()));
	  }
	if(clusterPID == 11)
	  pECalNRGVsEtaElecHist->Fill(pos.PseudoRapidity(),pECalE[i]);

	// Calc Track Pt wrt HBeam
	TVector3 longUnitB = hBeamMC.Vect().Unit();
	TVector3 normUnitB = TVector3(0,0,1).Cross(longUnitB).Unit();
	TVector3 sideUnitB = longUnitB.Cross(normUnitB);
	//TVector3 momB = partWRThBeam.Vect();
	TVector3 localMomB(mom.Dot(sideUnitB),mom.Dot(normUnitB),mom.Dot(longUnitB));

	// Test
	TVector3 localPosB(pos.Dot(sideUnitB),pos.Dot(normUnitB),pos.Dot(longUnitB));
	if(clusterPID == 22) testPosPhiVsEtaBeamHist->Fill(localPosB.PseudoRapidity(),TVector2::Phi_mpi_pi(localPosB.Phi()));
	if(clusterPID == 22) testPosYvsXBeamHist->Fill(localPosB.Px(),localPosB.Py());

	pECalPhiVsEtaBeamHist->Fill(localMomB.PseudoRapidity(),TVector2::Phi_mpi_pi(localMomB.Phi()));
	pECalPosYvsXBeamHist->Fill(localPosB.Px(),localPosB.Py());
	if(clusterPID != 11)
	  pECalPhiVsEtaNoElecBeamHist->Fill(localMomB.PseudoRapidity(),TVector2::Phi_mpi_pi(localMomB.Phi()));

	//Fill FastJet
	if(mom.PseudoRapidity() > 1.4 && mom.PseudoRapidity() < 4.0)
	  {
	    fastjet::PseudoJet p(localMomB.Px(),localMomB.Py(),localMomB.Pz(),pECalE[i]);
	    p.set_user_index(i);
	    if(clusterPID != 11) // clusterPID != 11
	      {
		particlesECal.push_back(p);
		particlesRecoECal.push_back(p);
	      }
	  }

	//cout << "Cluster " << i << " " << nECalPosX[i] << " " << nECalPosY[i] << " " << nECalPosZ[i] << " " << nECalE[i] << " " << pidBegin[i] << " " << pidEnd[i] << " " << pos.PseudoRapidity() << " " << pos.Phi() << endl;
      }

    // Loop Over Truth pECal Clusters
    for(unsigned int i=0; i<pECalTruthPosX.GetSize(); i++)
      {
	TVector3 pos(pECalTruthPosX[i],pECalTruthPosY[i],pECalTruthPosZ[i]);

	pECalTruthPosYvsXHist->Fill(pECalTruthPosX[i],pECalTruthPosY[i]);
	pECalTruthPosZHist->Fill(TMath::Abs(pECalTruthPosZ[i]));

	if(pECalTruthPosX[i] > -380.0 && pECalTruthPosX[i] < 180.0 && pECalTruthPosY[i] > -280.0 && pECalTruthPosY[i] < 280.0) // Kill Clusters in Gap
	  continue;

	Double_t px = pECalTruthE[i]*TMath::Sin(pos.Theta())*TMath::Cos(pos.Phi());
	Double_t py = pECalTruthE[i]*TMath::Sin(pos.Theta())*TMath::Sin(pos.Phi());
	Double_t pz = pECalTruthE[i]*TMath::Cos(pos.Theta());

	TVector3 mom(px,py,pz);

	pECalTruthNRGHist->Fill(pECalTruthE[i]);
	pECalTruthPhiVsEtaHist->Fill(pos.PseudoRapidity(),TVector2::Phi_mpi_pi(pos.Phi()));

	int clusterPID = pdg[pECalTruthSimID[i]];

	// Calc Track Pt wrt HBeam
	TVector3 longUnitB = hBeamMC.Vect().Unit();
	TVector3 normUnitB = TVector3(0,0,1).Cross(longUnitB).Unit();
	TVector3 sideUnitB = longUnitB.Cross(normUnitB);
	//TVector3 momB = partWRThBeam.Vect();
	TVector3 localMomB(mom.Dot(sideUnitB),mom.Dot(normUnitB),mom.Dot(longUnitB));
	TVector3 localPosB(pos.Dot(sideUnitB),pos.Dot(normUnitB),pos.Dot(longUnitB));

	pECalTruthPhiVsEtaBeamHist->Fill(localMomB.PseudoRapidity(),TVector2::Phi_mpi_pi(localMomB.Phi()));
	pECalTruthPosYvsXBeamHist->Fill(localPosB.Px(),localPosB.Py());

	//Fill FastJet
	if(mom.PseudoRapidity() > 1.4 && mom.PseudoRapidity() < 4.0)
	  {
	    fastjet::PseudoJet p(localMomB.Px(),localMomB.Py(),localMomB.Pz(),pECalTruthE[i]);
	    p.set_user_index(i);
	    if(clusterPID != 11) // clusterPID != 11
	      {
		particlesECalTruth.push_back(p);
		particlesRecoECalTruth.push_back(p);
	      }
	  }
      }

    // Loop Over pECal Reconstructed Hits
    for(unsigned int i=0; i<pECalRecPosX.GetSize(); i++)
      {
	TVector3 pos(pECalRecPosX[i],pECalRecPosY[i],pECalRecPosZ[i]);

	pECalRecPosYvsXHist->Fill(pECalRecPosX[i],pECalRecPosY[i]);
	pECalRecPosZHist->Fill(TMath::Abs(pECalRecPosZ[i]));

	Double_t px = pECalRecE[i]*TMath::Sin(pos.Theta())*TMath::Cos(pos.Phi());
	Double_t py = pECalRecE[i]*TMath::Sin(pos.Theta())*TMath::Sin(pos.Phi());
	Double_t pz = pECalRecE[i]*TMath::Cos(pos.Theta());

	TVector3 mom(px,py,pz);

	pECalRecNRGHist->Fill(pECalRecE[i]);
	pECalRecPhiVsEtaHist->Fill(pos.PseudoRapidity(),TVector2::Phi_mpi_pi(pos.Phi()));

	// Calc Track Pt wrt HBeam
	TVector3 longUnitB = hBeamMC.Vect().Unit();
	TVector3 normUnitB = TVector3(0,0,1).Cross(longUnitB).Unit();
	TVector3 sideUnitB = longUnitB.Cross(normUnitB);
	//TVector3 momB = partWRThBeam.Vect();
	TVector3 localMomB(mom.Dot(sideUnitB),mom.Dot(normUnitB),mom.Dot(longUnitB));
	TVector3 localPosB(pos.Dot(sideUnitB),pos.Dot(normUnitB),pos.Dot(longUnitB));

	pECalRecPhiVsEtaBeamHist->Fill(localMomB.PseudoRapidity(),TVector2::Phi_mpi_pi(localMomB.Phi()));
	pECalRecPosYvsXBeamHist->Fill(localPosB.Px(),localPosB.Py());
      }
    */
  
    /*
    // Cluster and Find Jets
    double R[2] = {1.0, 0.5};
    for(int rad=0; rad<1; rad++)
      {
	// Define Jet
	JetDefinition jet_def_akt(antikt_algorithm,R[rad]);

	// Cluster in Beam Frame
	ClusterSequence cs_akt_all_beam(particles, jet_def_akt);
	ClusterSequence cs_akt_charge_beam(particlesCharged, jet_def_akt);
	ClusterSequence cs_akt_reco_beam(particlesReco, jet_def_akt);
	ClusterSequence cs_akt_recoECal_beam(particlesRecoECal, jet_def_akt);
	ClusterSequence cs_akt_recoECalTruth_beam(particlesRecoECalTruth, jet_def_akt);

	// Set Min Jet Pt
	double ptmin = 1.0;

	// Get Jets
	vector<vector<PseudoJet> > jetVector;
	vector<PseudoJet> jets_akt_all_beam = sorted_by_pt(cs_akt_all_beam.inclusive_jets(ptmin));
	jetVector.push_back(jets_akt_all_beam);
	vector<PseudoJet> jets_akt_charge_beam = sorted_by_pt(cs_akt_charge_beam.inclusive_jets(ptmin));
	jetVector.push_back(jets_akt_charge_beam);
	vector<PseudoJet> jets_akt_reco_beam = sorted_by_pt(cs_akt_reco_beam.inclusive_jets(ptmin));
	jetVector.push_back(jets_akt_reco_beam);
	vector<PseudoJet> jets_akt_recoECal_beam = sorted_by_pt(cs_akt_recoECal_beam.inclusive_jets(ptmin));
	jetVector.push_back(jets_akt_recoECal_beam);
	vector<PseudoJet> jets_akt_recoECalTruth_beam = sorted_by_pt(cs_akt_recoECalTruth_beam.inclusive_jets(ptmin));
	jetVector.push_back(jets_akt_recoECalTruth_beam);

	// Loop Through MC Jets
	for(int jT=0; jT<5; jT++)
	  {
	    for(unsigned int jn=0; jn<jetVector[jT].size(); jn++)
	      {
		jetPtHist[jT]->Fill(jetVector[jT][jn].pt());
		jetEHist[jT]->Fill(jetVector[jT][jn].e());
		jetPtVsEtaHist[jT]->Fill(jetVector[jT][jn].eta(),jetVector[jT][jn].pt());
		jetEVsEtaHist[jT]->Fill(jetVector[jT][jn].eta(),jetVector[jT][jn].e());
		jetPhiVsEtaHist[jT]->Fill(jetVector[jT][jn].eta(),TVector2::Phi_mpi_pi(jetVector[jT][jn].phi()));
	      }
	  }


	// Loop Through All Particles MC Jets and Find Matching Reco
	for(unsigned int jn=0; jn<jetVector[0].size(); jn++)
	  {
	    Double_t partPt = jetVector[0][jn].pt();
	    Double_t partEta = jetVector[0][jn].eta();
	    Double_t partPhi = jetVector[0][jn].phi();
	    Double_t partE = jetVector[0][jn].e();

	    if(TMath::Abs(partEta) > 2.2) continue;

	    // Loop Through Reco Jets and Find Match
	    Double_t minDR = 100000.0;
	    Int_t minIndex = -1;
	    for(unsigned int jnn=0; jnn<jetVector[2].size(); jnn++)
	      {
		//Double_t recoPt = jetVector[2][jnn].pt();
		Double_t recoEta = jetVector[2][jnn].eta();
		Double_t recoPhi = jetVector[2][jnn].phi();
		//Double_t recoE = jetVector[2][jnn].e();

		Double_t deltaEta = partEta - recoEta;
		Double_t deltaPhi = TVector2::Phi_mpi_pi(partPhi - recoPhi);
		Double_t deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

		if(deltaR < minDR)
		  {
		    minDR = deltaR;
		    minIndex = jnn;
		  }
	      }

	    // Plot With Matching Jet
	    jetDeltaRHist[0]->Fill(minDR);

	    if(minIndex == -1) 
	      {
		jetMissingPtVsEtaHist[0]->Fill(partEta,partPt);
		jetMissingEVsEtaHist[0]->Fill(partEta,partE);
		jetMissingPhiVsEtaHist[0]->Fill(partEta,TVector2::Phi_mpi_pi(partPhi));
		continue;
	      }

	    if(minDR > 0.5) continue;

	    Double_t recoPt = jetVector[2][minIndex].pt();
	    Double_t recoE = jetVector[2][minIndex].e();

	    jetRecoVsPartPtHist[0]->Fill(partPt,recoPt);
	    jetRecoVsPartEHist[0]->Fill(partE,recoE);

	    if(partEta < -1.0)
	      {
		jetPtResVsPtEtaBin0Hist[0]->Fill(partPt,(recoPt-partPt)/partPt);
		jetEResVsEEtaBin0Hist[0]->Fill(partE,(recoE-partE)/partE);
	      }
	    if(partEta >= -1.0 && partEta < 1.0)
	      {
		jetPtResVsPtEtaBin1Hist[0]->Fill(partPt,(recoPt-partPt)/partPt);
		jetEResVsEEtaBin1Hist[0]->Fill(partE,(recoE-partE)/partE);
	      }
	    if(partEta >= 1.0)
	      {
		jetPtResVsPtEtaBin2Hist[0]->Fill(partPt,(recoPt-partPt)/partPt);
		jetEResVsEEtaBin2Hist[0]->Fill(partE,(recoE-partE)/partE);
	      }
	  }


	// Loop Through Charged MC Jets and Find Matching Reco
	for(unsigned int jn=0; jn<jetVector[1].size(); jn++)
	  {
	    Double_t partPt = jetVector[1][jn].pt();
	    Double_t partEta = jetVector[1][jn].eta();
	    Double_t partPhi = jetVector[1][jn].phi();
	    Double_t partE = jetVector[1][jn].e();

	    if(TMath::Abs(partEta) > 2.2) continue;

	    // Loop Through Reco Jets and Find Match
	    Double_t minDR = 100000.0;
	    Int_t minIndex = -1;
	    for(unsigned int jnn=0; jnn<jetVector[2].size(); jnn++)
	      {
		//Double_t recoPt = jetVector[2][jnn].pt();
		Double_t recoEta = jetVector[2][jnn].eta();
		Double_t recoPhi = jetVector[2][jnn].phi();
		//Double_t recoE = jetVector[2][jnn].e();

		Double_t deltaEta = partEta - recoEta;
		Double_t deltaPhi = TVector2::Phi_mpi_pi(partPhi - recoPhi);
		Double_t deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

		if(deltaR < minDR)
		  {
		    minDR = deltaR;
		    minIndex = jnn;
		  }
	      }

	    // Plot With Matching Jet
	    jetDeltaRHist[1]->Fill(minDR);

	    if(minIndex == -1) 
	      {
		jetMissingPtVsEtaHist[1]->Fill(partEta,partPt);
		jetMissingEVsEtaHist[1]->Fill(partEta,partE);
		jetMissingPhiVsEtaHist[1]->Fill(partEta,TVector2::Phi_mpi_pi(partPhi));
		continue;
	      }

	    if(minDR > 0.5) continue;

	    Double_t recoPt = jetVector[2][minIndex].pt();
	    Double_t recoE = jetVector[2][minIndex].e();

	    jetRecoVsPartPtHist[1]->Fill(partPt,recoPt);
	    jetRecoVsPartEHist[1]->Fill(partE,recoE);

	    if(partEta < -1.0)
	      {
		jetPtResVsPtEtaBin0Hist[1]->Fill(partPt,(recoPt-partPt)/partPt);
		jetEResVsEEtaBin0Hist[1]->Fill(partE,(recoE-partE)/partE);
	      }
	    if(partEta >= -1.0 && partEta < 1.0)
	      {
		jetPtResVsPtEtaBin1Hist[1]->Fill(partPt,(recoPt-partPt)/partPt);
		jetEResVsEEtaBin1Hist[1]->Fill(partE,(recoE-partE)/partE);
	      }
	    if(partEta >= 1.0)
	      {
		jetPtResVsPtEtaBin2Hist[1]->Fill(partPt,(recoPt-partPt)/partPt);
		jetEResVsEEtaBin2Hist[1]->Fill(partE,(recoE-partE)/partE);
	      }
	  }


	// Loop Through All MC Jets and Find Matching Reco+ECal
	for(unsigned int jn=0; jn<jetVector[0].size(); jn++)
	  {
	    Double_t partPt = jetVector[0][jn].pt();
	    Double_t partEta = jetVector[0][jn].eta();
	    Double_t partPhi = jetVector[0][jn].phi();
	    Double_t partE = jetVector[0][jn].e();

	    if(TMath::Abs(partEta) > 2.2) continue;

	    // Loop Through Reco Jets and Find Match
	    Double_t minDR = 100000.0;
	    Int_t minIndex = -1;
	    for(unsigned int jnn=0; jnn<jetVector[3].size(); jnn++)
	      {
		//Double_t recoPt = jetVector[3][jnn].pt();
		Double_t recoEta = jetVector[3][jnn].eta();
		Double_t recoPhi = jetVector[3][jnn].phi();
		//Double_t recoE = jetVector[3][jnn].e();

		Double_t deltaEta = partEta - recoEta;
		Double_t deltaPhi = TVector2::Phi_mpi_pi(partPhi - recoPhi);
		Double_t deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

		if(deltaR < minDR)
		  {
		    minDR = deltaR;
		    minIndex = jnn;
		  }
	      }

	    // Plot With Matching Jet
	    jetDeltaRHist[2]->Fill(minDR);

	    if(minIndex == -1) 
	      {
		jetMissingPtVsEtaHist[2]->Fill(partEta,partPt);
		jetMissingEVsEtaHist[2]->Fill(partEta,partE);
		jetMissingPhiVsEtaHist[2]->Fill(partEta,TVector2::Phi_mpi_pi(partPhi));
		continue;
	      }

	    if(minDR > 0.5) continue;

	    Double_t recoPt = jetVector[3][minIndex].pt();
	    Double_t recoE = jetVector[3][minIndex].e();

	    jetRecoVsPartPtHist[2]->Fill(partPt,recoPt);
	    jetRecoVsPartEHist[2]->Fill(partE,recoE);

	    if(partEta < -1.0)
	      {
		jetPtResVsPtEtaBin0Hist[2]->Fill(partPt,(recoPt-partPt)/partPt);
		jetEResVsEEtaBin0Hist[2]->Fill(partE,(recoE-partE)/partE);
	      }
	    if(partEta >= -1.0 && partEta < 1.0)
	      {
		jetPtResVsPtEtaBin1Hist[2]->Fill(partPt,(recoPt-partPt)/partPt);
		jetEResVsEEtaBin1Hist[2]->Fill(partE,(recoE-partE)/partE);
	      }
	    if(partEta >= 1.0)
	      {
		jetPtResVsPtEtaBin2Hist[2]->Fill(partPt,(recoPt-partPt)/partPt);
		jetEResVsEEtaBin2Hist[2]->Fill(partE,(recoE-partE)/partE);
	      }
	  }

	// Loop Through All MC Jets and Find Matching Reco+ECalTruth
	for(unsigned int jn=0; jn<jetVector[0].size(); jn++)
	  {
	    Double_t partPt = jetVector[0][jn].pt();
	    Double_t partEta = jetVector[0][jn].eta();
	    Double_t partPhi = jetVector[0][jn].phi();
	    Double_t partE = jetVector[0][jn].e();

	    if(TMath::Abs(partEta) > 2.2) continue;

	    // Loop Through Reco Jets and Find Match
	    Double_t minDR = 100000.0;
	    Int_t minIndex = -1;
	    for(unsigned int jnn=0; jnn<jetVector[4].size(); jnn++)
	      {
		//Double_t recoPt = jetVector[4][jnn].pt();
		Double_t recoEta = jetVector[4][jnn].eta();
		Double_t recoPhi = jetVector[4][jnn].phi();
		//Double_t recoE = jetVector[4][jnn].e();

		Double_t deltaEta = partEta - recoEta;
		Double_t deltaPhi = TVector2::Phi_mpi_pi(partPhi - recoPhi);
		Double_t deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

		if(deltaR < minDR)
		  {
		    minDR = deltaR;
		    minIndex = jnn;
		  }
	      }

	    // Plot With Matching Jet
	    jetDeltaRHist[3]->Fill(minDR);

	    if(minIndex == -1) 
	      {
		jetMissingPtVsEtaHist[3]->Fill(partEta,partPt);
		jetMissingEVsEtaHist[3]->Fill(partEta,partE);
		jetMissingPhiVsEtaHist[3]->Fill(partEta,TVector2::Phi_mpi_pi(partPhi));
		continue;
	      }

	    if(minDR > 0.5) continue;

	    Double_t recoPt = jetVector[4][minIndex].pt();
	    Double_t recoE = jetVector[4][minIndex].e();

	    jetRecoVsPartPtHist[3]->Fill(partPt,recoPt);
	    jetRecoVsPartEHist[3]->Fill(partE,recoE);

	    if(partEta < -1.0)
	      {
		jetPtResVsPtEtaBin0Hist[3]->Fill(partPt,(recoPt-partPt)/partPt);
		jetEResVsEEtaBin0Hist[3]->Fill(partE,(recoE-partE)/partE);
	      }
	    if(partEta >= -1.0 && partEta < 1.0)
	      {
		jetPtResVsPtEtaBin1Hist[3]->Fill(partPt,(recoPt-partPt)/partPt);
		jetEResVsEEtaBin1Hist[3]->Fill(partE,(recoE-partE)/partE);
	      }
	    if(partEta >= 1.0)
	      {
		jetPtResVsPtEtaBin2Hist[3]->Fill(partPt,(recoPt-partPt)/partPt);
		jetEResVsEEtaBin2Hist[3]->Fill(partE,(recoE-partE)/partE);
	      }
	  }
      } // rad
    */


    /*
    // Loop over Negative HCal
    Double_t nhcalHitE = 0.;
    for(unsigned int i=0; i<hcalEndNHitE.GetSize(); i++)
      {
	hcalEndNHitPosXY->Fill(hcalEndNHitPosX[i],hcalEndNHitPosY[i]);
	nhcalHitE += hcalEndNHitE[i];
      }

    hcalEndNHitTotE->Fill(nhcalHitE);

    // Loop over Barrel HCal
    hcalBarrelNumClust->Fill(hcalBarrelClustE.GetSize());
    hcalBarrelNumRecHits->Fill(hcalBarrelRecHitsE.GetSize());

    Double_t bhcalClustE = 0.;
    Double_t bhcalRecHitE = 0.;
    for(unsigned int i=0; i<hcalBarrelClustE.GetSize(); i++)
      {
	Double_t thetaToEta = (-1.0)*TMath::Log(TMath::Tan(0.5*hcalBarrelClustInTheta[i]));

	//cout << hcalBarrelClustInTheta[i] << " " << thetaToEta << endl;

	hcalBarrelNumHits->Fill(hcalBarrelClustNHits[i]);
	hcalBarrelClustNRG->Fill(hcalBarrelClustE[i]);
	hcalBarrelClustYvsX->Fill(hcalBarrelClustPosX[i],hcalBarrelClustPosY[i]);
	hcalBarrelClustPhiVsEta->Fill(thetaToEta,hcalBarrelClustInPhi[i]);
	hcalBarrelClustZ->Fill(hcalBarrelClustPosZ[i]);
	
	bhcalClustE += hcalBarrelClustE[i];
      }

    for(unsigned int i=0; i<hcalBarrelRecHitsE.GetSize(); i++)
      {
	bhcalRecHitE += hcalBarrelRecHitsE[i];
      }

    hcalBarrelClustTotE->Fill(bhcalClustE);
    hcalBarrelRecHitsTotE->Fill(bhcalRecHitE);

    numReco->Fill(recoMomX.GetSize());
    numRecoAssoc->Fill(assocRecID.GetSize());

    // Loop Over Reconstructed Charged Particles
    for(unsigned int i=0; i<recoMomX.GetSize(); i++)
      {
	TVector3 mom(recoMomX[i],recoMomY[i],recoMomZ[i]);
	double E = std::sqrt(recoMomX[i]*recoMomX[i] + recoMomY[i]*recoMomY[i] + recoMomZ[i]*recoMomZ[i] + recoM[i]*recoM[i]);
	recoPt->Fill(mom.Perp());
	recoEta->Fill(mom.PseudoRapidity());
	recoPhi->Fill(mom.Phi());
	recoNRG->Fill(E);
	recoPhiVsEta->Fill(mom.PseudoRapidity(),mom.Phi());
      }

    // Loop Over Reconstructed Associations
    for(unsigned int i=0; i<assocRecID.GetSize(); i++)
      {
	auto rI = assocRecID[i];
	auto sI = assocSimID[i];

	TVector3 momRec(recoMomX[rI],recoMomY[rI],recoMomZ[rI]);
	TVector3 momSim(mcMomX[sI],mcMomY[sI],mcMomZ[sI]);

	assocDiffPt->Fill(momRec.Perp() - momSim.Perp());
	assocResPt->Fill((momRec.Perp() - momSim.Perp())/momSim.Perp());

	// Angular Res
	double deltaEta = momRec.PseudoRapidity() - momSim.PseudoRapidity();
	double deltaPhi = TVector2::Phi_mpi_pi(momRec.Phi() - momSim.Phi());
	double deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

	assocDiffDeltaEta->Fill(deltaEta);
	assocDiffDeltaPhi->Fill(deltaPhi);
	assocDiffDeltaR->Fill(deltaR);
      }
    */

    NEVENTS++;
  }

  cout << "Total Particles = " << TOTPARTS << endl;
  cout << "Total Reco/Local Jet Pairs with Mismatched Energy = " << MISMATCHJETE << endl;

  //mcEta->Write();
  //mcNRG->Write();

  ofile->Write();
  ofile->Close();
  return 0;
}
