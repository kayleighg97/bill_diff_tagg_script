#include "diff_tagg_ana.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

#include <stdio.h>

#include <fun4all/Fun4AllHistoManager.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <g4eval/CaloEvalStack.h>
//#include <calobase/RawCluster.h>
#include <g4eval/CaloRawClusterEval.h>
//#include <calobase/RawClusterContainer.h>
#include <fun4all/SubsysReco.h>
#include <calobase/RawTowerContainer.h>

// G4Hits includes
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

// Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TMath.h>

#include <cassert>
#include <sstream>
#include <string>
#include <iostream>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

/// Tracking includes
#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <g4eval/SvtxEvalStack.h>

// G4Cells includes
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>

/// HEPMC truth includes
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

/// Fun4All includes
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Reco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <pdbcalbase/PdbParameterMap.h>
#include <phparameter/PHParameters.h>
#include <pdbcalbase/PdbParameterMapContainer.h>

// Tower includes
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

/// HEPMC truth includes
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>



//#include <coresoftware/blob/master/simulation/g4simulation/g4eval/SvtxEvalStack.h>
//#include "/cvmfs/eic.opensciencegrid.org/ecce/gcc-8.3/release/release_new/new.1/include/g4eval/SvtxEvalStack.h"

using namespace std;

diff_tagg_ana::diff_tagg_ana(const std::string &name, const std::string& filename):
 SubsysReco(name)
 , outfilename(filename)
{
  std::cout << "Diff_Tagg_example::Diff_Tagg_example(const std::string &name) Calling ctor" << std::endl;
	 
  _caloevalstackFEMC=nullptr;
  _caloevalstackEEMC=nullptr;
  _caloevalstackBECAL=nullptr;
  _svtxEvalStack = nullptr;

  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(m_RandomGenerator, seed);

  initializeVariables();
  initializeTrees();

}




//____________________________________________________________________________..
diff_tagg_ana::~diff_tagg_ana()
{
  delete tree;
  if(_caloevalstackBECAL) delete _caloevalstackBECAL;
  if(_caloevalstackFEMC) delete _caloevalstackFEMC;
  if(_caloevalstackEEMC) delete _caloevalstackEEMC;
  if(_svtxEvalStack) delete _svtxEvalStack;
  gsl_rng_free(m_RandomGenerator);

  std::cout << "diff_tagg_ana::~diff_tagg_ana() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int diff_tagg_ana::Init(PHCompositeNode *topNode)
{

  static_event_counter = 0;

  hm = new Fun4AllHistoManager(Name());
  // create and register your histos (all types) here
  // TH1 *h1 = new TH1F("h1",....)
  // hm->registerHisto(h1);
  outfile = new TFile(outfilename.c_str(), "RECREATE");
  g4hitntuple = new TNtuple("hitntup", "G4Hits", "x0:y0:z0:x1:y1:z1:edep");

  std::cout << "diff_tagg_ana::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  event_itt = 0;

  //**************
  // ZDC

  gDirectory->mkdir("ZDC");
  gDirectory->cd("ZDC");

  h2_ZDC_XY_g = new TH2F("ZDC_XY_g", "ZDC XY", 200, -200, 200, 200, -50, 50);
  h2_ZDC_XY_g_double = new TH2F("ZDC_XY_g_double", "ZDC XY Double gamma", 200, -200, 200, 200, -50, 50);

  h2_ZDC_XY_l = new TH2F("ZDC_XY_l", "ZDC XY_l", 200, -50, 50, 200, -50, 50);
  h2_ZDC_XY_l_double = new TH2F("ZDC_XY_l_double", "ZDC XY Double gamma", 200, -50, -50, 200, -50, 50);

  h1_E_dep  = new TH1F("E_dep", "E Dependence", 120, 0.0, 60.0);

  h1_E_dep_smeared = new TH1F("E_dep_smeared", "E Dependence Smeared", 120, 0.0, 60.0);

  gDirectory->cd("/");

  //**************
  // RP

  // This is a test for the git complict
  // This is a test for the git complict
  // This is a test for the git complict
	
  gDirectory->mkdir("RP") ;
  gDirectory->cd("RP") ;

  // This is a test for the git complict

//  h2_RP_XY_g = new TH2F("RP_XY_global", "RP_XY_global", 100, -500, 500, 100, -500, 500); 
//  h2_RP_XY_l = new TH2F("RP_XY_local", "RP_XY_local", 100, -50, 50, 100, -50, 50); 
  h2_RP_XY_signal = new TH2F("RP_XY_signal", "RP_XY_signal", 100, -50, 50, 100, -50, 50); 

  h2_RP_XY_g = new TH2F("RP_XY_global", "RP_XY_global", 200, -150, 150, 100, -100, 100); 
//  h2_RP_XY_l = new TH2F("RP_XY_local", "RP_XY_local", 200, -50, 50, 200, -50, 50); 
  h2_RP_XY_l = new TH2F("RP_XY_local", "RP_XY_local", 200, -5, 5, 200, -5, 5); 

  gDirectory->cd("/");

  //**************
  // Low Q2
  gDirectory->mkdir("LowQ2");
  gDirectory->cd("LowQ2");

  h2_lowQ2_XY = new TH2F("h2_lowQ2_XY", "h2_lowQ2_XY", 200, -80, -20, 200, -20, 20); 
  h_Q2_truth = new TH1F("h_Q2_truth", "h_Q2_truth", 200, 1e-9, 1); 
  h_Q2_truth_LowQ2tag = new TH1F("h_Q2_truth_LowQ2tag", "h_Q2_truth_LowQ2tag", 200, 0, 5); 
  h2_Q2_pos = new TH2F("h2_Q2_truth_pos", "h_Q2_truth_pos", 200, -80, -20, 200, 0, 5); 
//  h2_Q2_pos = new TH2F("h2_Q2_truth_pos", "h_Q2_truth_pos", 20, 0, 100, 20, 0, 5); 

  h2_Q2_truth_E = new TH2F("h2_Q2_truth_E", "h2_Q2_truth_E", 200, 0, 5, 200, 0, 18); 
  h2_pos_mom = new TH2F("h2_pos_mom", "h2_pos_mom", 200, -80, -20, 200, 0, 0.00005); 
//  h2_Q2_theta = new TH2F("h2_Q2_truth_theta", "h_Q2_truth_theta", 200, 0, 3.14, 200, 0, 0.00005); 

  
  // ----------------------------------
  // Low Q2 tagger


   const Int_t nbins = 150;
   Double_t xmin = 1e-9;
   Double_t xmax = 1e0;
   Double_t logxmin = log10(xmin);
   Double_t logxmax = log10(xmax);
   Double_t binwidth = (logxmax-logxmin)/nbins;
   Double_t xbins[nbins+1];
   xbins[0] = xmin;
   for (Int_t i=1;i<=nbins;i++) {
      xbins[i] = xmin + TMath::Power(10,logxmin+i*binwidth);
   }
 
//   TH1F *h = new TH1F("h","hist with log x axis",nbins,xbins);
 
  h_Q2 = new TH1F("h_Q2", "h2_Q2", nbins, xbins); 


//  h_log_Q2 = new TH1F("h_log_Q2", "h2_log_Q2", 200, -10, 0); 
//  h_log_Q2_LowQ2tag = new TH1F("h_log_Q2_LowQ2tag", "h2_log_LowQ2tag", 200, -10, 0); 



  h_log_Q2 = new TH1F("h_log_Q2", "h2_log_Q2", nbins, xbins); 
  h_log_Q2_LowQ2tag = new TH1F("h_log_Q2_LowQ2tag", "h2_log_LowQ2tag", nbins, xbins); 

  h_E = new TH1F("h_E", "h_E", 200, 4, 24); 
  h_E_LowQ2tag = new TH1F("h_E_LowQ2tag", "h_E_LowQ2tag", 200, 4, 24); 

  h_eta = new TH1F("h_eta", "h_eta", 200, -14, -4); 
  h_eta_LowQ2tag = new TH1F("h_eta_LowQ2tag", "h_eta_LowQ2tag", 200, -14, -4); 

  h_polar = new TH1F("h_polar", "h_polar", 100, 0, 15); 
  h_polar_LowQ2tag = new TH1F("h_polar_LowQ2tag", "h_polar_LowQ2tag", 100, 0, 15); 

  h2_E_Q2 = new TH2F("h2_E_Q2", "h2_E_Q2", 100, 0, 20, nbins, xbins);
  h2_E_Q2_LowQ2tag = new TH2F("h2_E_Q2_LowQ2tag", "h2_E_Q2_LowQ2tag", 100, 0, 20,  nbins, xbins); 
 

  gDirectory->cd("/");

  //**************
  // B0

  gDirectory->mkdir("B0");
  gDirectory->cd("B0");

  h2_B0_XY = new TH2F("B0_XY_global", "B0_XY_global", 50, -50, 0, 50, -25, 25); 
  h2_B0_XY_l = new TH2F("B0_XY_local", "B0_XY_local", 50, -25, 25, 50, -25, 25); 

  gDirectory->cd("/");

  //***********************8

  m_mpi = -99;
  m_process_id = -99;
  m_truthenergy = -99;
  m_trutheta = -99;
  m_truthphi = -99;
  m_truthp = -99;
  m_truthpx = -99;
  m_truthpy = -99;
  m_truthpz = -99;
  m_truthpt = -99;
  m_numparticlesinevent = -99;
  m_truthpid = -99;



  ///**********************************/
  // Parameter definition

  mElec = 0.000510998950;
  mProt = 0.93827208816;

  // Define beam 4 vectors
  e_beam_energy = 18;
  e_beam_pmag = sqrt(pow(e_beam_energy,2)-pow(mElec,2));
  ion_beam_energy = 275;
  ion_beam_pmag = sqrt((pow(ion_beam_energy,2)-pow(mProt,2)));
  crossing_angle = 0.025; 

  //Double_t Pi = TMath::ACos(-1);
  eBeam4Vect.SetPxPyPzE(0,0,-1*e_beam_pmag,e_beam_energy);
  pBeam4Vect.SetPxPyPzE(-ion_beam_pmag*TMath::Sin(crossing_angle),0,ion_beam_pmag*TMath::Cos(crossing_angle),ion_beam_energy);

  /**********************************/

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int diff_tagg_ana::InitRun(PHCompositeNode *topNode)
{
  if( static_event_counter == 0) {
  
  	encloseure_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_hFarFwdBeamLineEnclosure_0");

        encloseure_nodeparams->Print();

  	if (encloseure_nodeparams)
  	{
  	   Enclosure_params.FillFrom(encloseure_nodeparams, 0);
  	} else {
  	   cerr << "There is a issue finding the detector paramter node!" << endl;
  	}


  	beamlinemagnet_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_BEAMLINEMAGNET");

//	beamlinemagnet_nodeparams->print();

  	if (encloseure_nodeparams)
  	{
  	   BeamLineMagnet_params.FillFrom(beamlinemagnet_nodeparams, 0);
  	} else {
  	   cerr << "There is a issue finding the detector paramter node!" << endl;
  	}

//	exit(0);

	std::cout << "diff_tagg_ana::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
	rp_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_rpTruth");
	rp2_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_rpTruth2");
	b0_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_b0Truth_0");
	zdc_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_ZDCsurrogate");


  	static_event_counter++;

       /// Determining which IP design
	if (zdc_nodeparams) {
  	   if (rp2_nodeparams) {
		IP_design = "IP8";
	   } else {
		IP_design = "IP6";
           }
	} else {
	   IP_design = "UNKNOWN";
        }
 	std::cout << "IP design: " << IP_design << std::endl;
  }

  cout << " END initialization" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int diff_tagg_ana::process_event(PHCompositeNode *topNode)
{
  nHits=0;
  hitsEEMC=0;
  hitsFEMC=0;
  hitsBECAL=0;
  RP1 = 0;
  RP2 =0;
  RPhits=0;
  B0hits=0;
  ntr=0;
  cout<<" event = "<<event_itt<<endl;
  ZDC_hit = 0;

  event_itt++; 
 
  if(event_itt%100 == 0)
     std::cout << "Event Processing Counter: " << event_itt << endl;
	
  PHG4TruthInfoContainer* m_TruthInfoContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  PHG4Particle* primary = m_TruthInfoContainer->GetParticle(4);
  // std::cout << "pid 1 = "<< primary->get_pid() << std::endl;
  //PHG4TruthInfoContainer::ConstRange particles = m_TruthInfoContainer->GetPrimaryParticleRange();
  Epx = primary->get_px();
  Epy = primary->get_py();
  Epz = primary->get_pz();

  //if (m_TruthInfoContainer->GetParticle(0)!=0){
  primary = m_TruthInfoContainer->GetParticle(3);
  //std::cout << "pid 2 = "<< primary->get_pid() << std::endl;
  Ppx = primary->get_px();
  Ppy = primary->get_py();
  Ppz = primary->get_pz();
  // }
  primary = m_TruthInfoContainer->GetParticle(2); //PID 11
  //std::cout<< "pid 3 = " << primary->get_pid()<< std::endl;
  Gpx = primary->get_px(); //GpX = decay electron, Pos_ will be positron
  Gpy = primary->get_py();
  Gpz = primary->get_pz();

  primary = m_TruthInfoContainer->GetParticle(1);
  Pos_px = primary->get_px(); //GpX = decay electron, Pos_ will be positron                                                                                                                               
  Pos_py = primary->get_py();
  Pos_pz = primary->get_pz();


  if(!_caloevalstackFEMC){
    _caloevalstackFEMC = new CaloEvalStack(topNode, "FEMC");
    _caloevalstackFEMC->set_strict(true);
  }
  else{
    _caloevalstackFEMC->next_event(topNode);
  }

  if(!_caloevalstackEEMC){
    _caloevalstackEEMC = new CaloEvalStack(topNode, "EEMC");
    _caloevalstackEEMC->set_strict(true);
  }
  else{
    _caloevalstackEEMC->next_event(topNode);
  }


  if(!_caloevalstackBECAL){
    _caloevalstackBECAL = new CaloEvalStack(topNode, "BECAL");
    _caloevalstackBECAL->set_strict(true);
  }
  else{
    _caloevalstackBECAL->next_event(topNode);
  }

  std::string caloName;
  caloName="CLUSTER_FEMC";
  process_ClusterCalo(topNode,caloName);
  caloName="CLUSTER_EEMC";
  process_ClusterCalo(topNode,caloName);
  caloName="CLUSTER_BECAL";
  process_ClusterCalo(topNode,caloName);
  // process_g4hits(topNode);


  process_tracks(topNode);
  process_RomanPots(topNode);
  process_B0(topNode);
  process_Q2tagger(topNode);
  process_ZDC(topNode);
	
  process_PHG4Truth_Primary_Particles(topNode);

  process_PHG4Truth(topNode);
  
  process_g4hits_LowQ2Tagger(topNode);

  tree->Fill();
	
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int diff_tagg_ana::ResetEvent(PHCompositeNode *topNode)
{
//  std::cout << "diff_tagg_ana::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
//
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int diff_tagg_ana::EndRun(const int runnumber)
{
  std::cout << "diff_tagg_ana::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int diff_tagg_ana::End(PHCompositeNode *topNode)
{
  std::cout << "diff_tagg_ana::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  if(_caloevalstackBECAL) delete _caloevalstackBECAL;
  if(_caloevalstackFEMC) delete _caloevalstackFEMC;
  if(_caloevalstackEEMC) delete _caloevalstackEEMC;

//  h2_ZDC_XY->Write();
//  h2_ZDC_XY_double->Write();
//  
//  h1_E_dep->Write();
//  h1_E_dep_smeared->Write();

  outfile->cd();
	
  tree->Write();
	
  g4hitntuple->Write();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int diff_tagg_ana::Reset(PHCompositeNode *topNode)
{
 std::cout << "diff_tagg_ana::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void diff_tagg_ana::Print(const std::string &what) const
{
  std::cout << "diff_tagg_ana::Print(const std::string &what) const Printing info for " << what << std::endl;
}


//***************************************************
//


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int diff_tagg_ana::process_PHG4Truth_Primary_Particles(PHCompositeNode* topNode) {

 PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  if (!truthinfo)
  {
    cout << PHWHERE
         << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles"
         << endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  /// Get the primary particle range
  ///PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();

  /// Loop over the G4 truth (stable) particles
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second;
       ++iter)
  {
    /// Get this truth particle
    const PHG4Particle *truth = iter->second;

    /// Get this particles momentum, etc.
    m_truthpx = truth->get_px();
    m_truthpy = truth->get_py();
    m_truthpz = truth->get_pz();
    m_truthp = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy + m_truthpz * m_truthpz);
    m_truthenergy = truth->get_e();

    m_truthpt = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy);

    m_truthphi = atan(m_truthpy / m_truthpx);

    float m_trutheta = atanh(m_truthpz / m_truthenergy);
    /// Check for nans
    if (m_trutheta != m_trutheta)
      m_trutheta = -99;

    float m_truthpid = truth->get_pid();

	m_truthpid = m_truthpid;

    cout << setprecision(10) << "truth: " << m_truthpid << "  " << m_truthpx << "  " << m_truthpy  << "  " << m_truthpz << endl;

    /// Fill the g4 truth tree
//    m_truthtree->Fill();
  }

  return Fun4AllReturnCodes::EVENT_OK;

}



//***************************************************

int diff_tagg_ana::process_PHG4Truth(PHCompositeNode* topNode) {

 PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  if (!truthinfo)
  {
    cout << PHWHERE
         << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles"
         << endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  /// Get the primary particle range
  ///PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
  PHG4TruthInfoContainer::Range range = truthinfo->GetParticleRange();

  /// Loop over the G4 truth (stable) particles
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second;
       ++iter)
  {
    /// Get this truth particle
    const PHG4Particle *truth = iter->second;

    /// Get this particles momentum, etc.
    m_truthpx = truth->get_px();
    m_truthpy = truth->get_py();
    m_truthpz = truth->get_pz();
    m_truthp = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy + m_truthpz * m_truthpz);
    m_truthenergy = truth->get_e();

    m_truthpt = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy);

    m_truthphi = atan(m_truthpy / m_truthpx);

    float m_trutheta = atanh(m_truthpz / m_truthenergy);
    /// Check for nans
    if (m_trutheta != m_trutheta)
      m_trutheta = -99;
    float m_truthpid = truth->get_pid();

	m_truthpid = m_truthpid;

//    cout << "truth: " << m_truthpid << "  " << m_truthpx << "  " << m_truthpy 
//         << "  " << m_truthpz << endl;

    /// Fill the g4 truth tree
//    m_truthtree->Fill();
  }

  return Fun4AllReturnCodes::EVENT_OK;

}

//***************************************************

int diff_tagg_ana::process_ZDC(PHCompositeNode* topNode)
{

   ostringstream nodename;
 
   nodename.str("");
   nodename << "G4HIT_" << "ZDC";
 
  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());

  //float smeared_E;

  if (hits) {
//    // this returns an iterator to the beginning and the end of our G4Hits
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++){
//	HIT_IN_ZDC=true;
	ZDC_hit++;
      	g4hitntuple->Fill(hit_iter->second->get_x(0),
                        hit_iter->second->get_y(0),
                        hit_iter->second->get_z(0),
                        hit_iter->second->get_x(1),
                        hit_iter->second->get_y(1),
                        hit_iter->second->get_z(1),
                        hit_iter->second->get_edep());

//       h2_ZDC_XY->Fill(hit_iter->second->get_x(0), hit_iter->second->get_y(0)); 
// //
//       //smeared_E = ZDC_Energy_Smear_EMCAL(hit_iter->second->get_edep());
// //
//       if (ZDC_hit == 2 ) {

// //      cout << hit_iter->second->get_x(0)-90 << "   " << hit_iter->second->get_y(0) << endl;

//         h2_ZDC_XY_g_double->Fill(hit_iter->second->get_x(0), hit_iter->second->get_y(0)); 
// //      h1_E_dep->Fill(hit_iter->second->get_edep()); 
// //
//         h1_E_dep->Fill(hit_iter->second->get_edep()); 
//         h1_E_dep_smeared->Fill(smeared_E); 
// //
//       }
    }
  }

//       //*******************************************/
//       // Method 1 
// //
// //      float det_xCent = Enclosure_params.get_double_param("place_x") + ZDC_params.get_double_param("place_x");
// //      float det_zCent = Enclosure_params.get_double_param("place_z") + ZDC_params.get_double_param("place_z");
// //      float det_tilt = ZDC_params.get_double_param("rot_y")/180. * TMath::Pi(); // in Rad
// //
// //      float det_rot = atan( det_xCent / det_zCent);  // in Rad
// //
// //      float local_x = Get_Local_X(hit_iter->second->get_x(0), hit_iter->second->get_y(0), hit_iter->second->get_z(0), det_tilt, det_rot) ;
// //      float local_y = Get_Local_Y(hit_iter->second->get_x(0), hit_iter->second->get_y(0), hit_iter->second->get_z(0), det_tilt, det_rot) ;



//       //*******************************************/
//       // Method 2    
//       //

//       PHParameters ZDC_params{"PHG4RP"};
            
//       if (zdc_nodeparams)
//       {
//          ZDC_params.FillFrom(zdc_nodeparams, 0);
//       } else {
//          cerr << "There is a issue finding the detector paramter node!" << endl;
//       }

// //      cout << "z original: " << Enclosure_params.get_double_param("place_z")  + ZDC_params.get_double_param("place_z") << endl;

//       float det_x_pos = Enclosure_params.get_double_param("place_x")  + ZDC_params.get_double_param("place_x");
//       float det_z_pos = Enclosure_params.get_double_param("place_z")  + ZDC_params.get_double_param("place_z");

//       ZDC_params.set_double_param("place_x", det_x_pos); 
//       ZDC_params.set_double_param("place_z", det_z_pos); 

//       float local_x = Get_Local_X(hit_iter->second->get_x(0), hit_iter->second->get_y(0), hit_iter->second->get_z(0), ZDC_params);
//       float local_y = hit_iter->second->get_y(0);

// //      cout << "x: " << hit_iter->second->get_x(0) << "  " << det_x_pos << endl;
// //      cout << "z: " << hit_iter->second->get_z(0) << "  " << det_z_pos << endl;

//       h2_ZDC_XY_l->Fill(local_x, local_y); 

//     }
//   }

// //  cout << "BB" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//***************************************************
// Getting the RomanPots hits

int diff_tagg_ana::process_RomanPots(PHCompositeNode* topNode)
{
  ostringstream nodename;

  // loop over the G4Hits
  nodename.str("");
  nodename << "G4HIT_" << "rpTruth_VirtSheet";

  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  if (!truthinfo)
  {
    cout << PHWHERE
         << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles"
         << endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

//   /// Get the primary particle range
//   PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
  if (hits) {
//    // this returns an iterator to the beginning and the end of our G4Hits
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
    //int counter = 0;
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++) {
      PHParameters RP_1_params{"PHRP_1"};
      if (rp_nodeparams)
	{
	  RP_1_params.FillFrom(rp_nodeparams, 0);
	  /*
	  cout << RP_1_params.get_double_param("Layer1_pos_x") << endl;
	  cout << RP_1_params.get_double_param("Layer1_pos_z") << endl;
	  cout << RP_1_params.get_double_param("Layer1_rot_y") << endl;
	  cout << RP_1_params.get_double_param("Layer2_pos_x") << endl;
	  cout << RP_1_params.get_double_param("Layer2_pos_z") << endl;
	  cout << RP_1_params.get_double_param("Layer2_rot_y") << endl;
	  cout << endl;
	  */	
	} else {
	cerr << "There is a issue finding the detector paramter node!" << endl;
      }
      
      PHParameters RP_2_params{"PHRP_2"};
      if (IP_design == "IP8"){
	if (rp2_nodeparams)
	  {
	    RP_2_params.FillFrom(rp2_nodeparams, 0);
	    cout << RP_2_params.get_double_param("Layer3_pos_x") << endl;
	    cout << RP_2_params.get_double_param("Layer3_pos_z") << endl;
	    cout << RP_2_params.get_double_param("Layer3_rot_y") << endl;
	    cout << RP_2_params.get_double_param("Layer4_pos_x") << endl;
	    cout << RP_2_params.get_double_param("Layer4_pos_z") << endl;
	    cout << RP_2_params.get_double_param("Layer4_rot_y") << endl;
	    cout << endl;
	  } else {
	  cerr << "There is a issue finding the detector paramter node!" << endl;
	}
      }
      
      //Converting to the global coordinate, where the forward vacuum encloseure must be taken into account  
      //float det_x_pos = Enclosure_params.get_double_param("place_x")  + RP_1_params.get_double_param("Layer1_pos_x");
      //float det_z_pos = Enclosure_params.get_double_param("place_z")  + RP_1_params.get_double_param("Layer1_pos_z");
      
      //RP_1_params.set_double_param("Layer1_pos_x", det_x_pos); 
      //RP_1_params.set_double_param("Layer1_pos_z", det_z_pos);
      
      
      if (RP_1_params.get_double_param("Layer1_rot_y") > 1.5) {
	float deg_to_rad = RP_1_params.get_double_param("rot_y") * TMath::Pi() / 180.;
	RP_1_params.set_double_param("Layer1_rot_y", deg_to_rad); 
      }
      
      float local_x = hit_iter->second->get_x(0);
      
      //if (hit_iter->second->get_trkid() == 1){
      if(TMath::Abs( hit_iter->second->get_z(0) - RP_1_params.get_double_param("Layer1_pos_z" ))<50){
	RPind[RPhits ] = 1;
	local_x = hit_iter->second->get_x(0) - RP_1_params.get_double_param("Layer1_pos_x");
      }else if(TMath::Abs( hit_iter->second->get_z(0) - RP_1_params.get_double_param("Layer2_pos_z"))<50){ 
	RPind[RPhits ] = 2;
	local_x = hit_iter->second->get_x(0) - RP_1_params.get_double_param("Layer2_pos_x");
      }else if(TMath::Abs( hit_iter->second->get_z(0) - RP_2_params.get_double_param("Layer1_pos_z" ))<50){
	RPind[RPhits ] = 3;
	local_x = hit_iter->second->get_x(0) - RP_2_params.get_double_param("Layer1_pos_x");
      }else if(TMath::Abs( hit_iter->second->get_z(0) - RP_2_params.get_double_param("Layer2_pos_z"))<50){
	RPind[RPhits ] = 4;
	local_x = hit_iter->second->get_x(0) - RP_2_params.get_double_param("Layer2_pos_x");
      }
      //float local_x = Get_Local_X(hit_iter->second->get_x(0), hit_iter->second->get_y(0), hit_iter->second->get_z(0), RP_1_params);
     
      float local_y = hit_iter->second->get_y(0);
      float local_z = hit_iter->second->get_z(0);
      
      //cout << local_x << " " << local_y << " " << local_z << endl;
      h2_RP_XY_g->Fill(hit_iter->second->get_x(0), hit_iter->second->get_y(0));
      h2_RP_XY_l->Fill(local_x, local_y);
	
      RPx[RPhits] = hit_iter->second->get_x(0);
      RPy[RPhits] = hit_iter->second->get_y(0);
      RPz[RPhits] = hit_iter->second->get_z(0);
	
      //RPxloc[RPhits] = hit_iter->second->get_local_x(0);
      //RPyloc[RPhits] = hit_iter->second->get_local_y(0);
      //RPzloc[RPhits] = hit_iter->second->get_local_z(0);
	
      RPxloc[RPhits] = local_x;
      RPyloc[RPhits] = local_y;
      RPzloc[RPhits] = local_z;
	
      RP_px[RPhits] = hit_iter->second->get_px(0);
      RP_py[RPhits] = hit_iter->second->get_py(0);
      RP_pz[RPhits] = hit_iter->second->get_pz(0);
      RP_edep[RPhits] = hit_iter->second->get_edep();
      
      RPtrkid[RPhits] = hit_iter->second->get_trkid();
      
      RPhits++;
      //}//if hit id loop
    }//for hit loop
  }//if hits loop
// 	PHG4Particle* g4particle = truthinfo->GetParticle(hit_iter->second->get_trkid());

// 	g4particle->get_px();

// 	//---------------------------------

//  	 /// Loop over the G4 truth (stable) particles
//  	 for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
//  	      iter != range.second;
//  	      ++iter)
//  	 {
//  	   /// Get this truth particle
//  	   const PHG4Particle *truth = iter->second;


//            float m_trutheta = atanh(m_truthpz / m_truthenergy);
//     	   /// Check for nans
//     	   if (m_trutheta != m_trutheta)
//       	   m_trutheta = -99;

//           truth->get_pid();
// //	   cout << "Particle in Roman Pot: "<< truth->get_pid() << endl;
// //	   cout << "Particle barcode: "<< truth->get_barcode() << endl;
// //	   cout << "Particle primary ID: "<< truth->get_primary_id() << endl;
// //   	   float m_truthpid = truth->get_pid();



// 	   /// Generic filling algorithm
// 	   if (hit_iter->second->get_hit_id() == 3 || hit_iter->second->get_hit_id() == 4294967297) {

// 	   PHParameters RP_1_params{"PHRP_1"};
	   
// 	   if (rp_nodeparams)
// 	   {
// 	      RP_1_params.FillFrom(rp_nodeparams, 0);
// 	   } else {
// 	      cerr << "There is a issue finding the detector paramter node!" << endl;
// 	   }


// 	   if (hit_iter->second->get_z(0) > Enclosure_params.get_double_param("place_z") + RP_1_params.get_double_param("Layer1_pos_z") - 50 &&    hit_iter->second->get_z(0) < Enclosure_params.get_double_param("place_z") + RP_1_params.get_double_param("Layer1_pos_z") + 50 ) {
 
// //           return 0;

//            h2_RP_XY_g->Fill(hit_iter->second->get_x(0), hit_iter->second->get_y(0));

// //	   float local_x = Get_Local_X(hit_iter->second->get_x(0), hit_iter->second->get_y(0), hit_iter->second->get_z(0), rp_nodeparams);


// 	    PHParameters B0_1_params{"PHB0_1"};
	
// 	    if (b0_nodeparams)
// 	    {
// 	       B0_1_params.FillFrom(b0_nodeparams, 0);
// 	    } else {
// 	       cerr << "There is a issue finding the detector paramter node!" << endl;
// 	    }

// 	   //******************************
// 	   /// Converting to the global coordinate, where the forward vacuum encloseure must be taken into account  


// 	   float det_x_pos = Enclosure_params.get_double_param("place_x")  + RP_1_params.get_double_param("place_x");
// 	   float det_z_pos = Enclosure_params.get_double_param("place_z")  + RP_1_params.get_double_param("place_z");

// 	   RP_1_params.set_double_param("place_x", det_x_pos); 
// 	   RP_1_params.set_double_param("place_z", det_z_pos); 

// 	   //******************************
// 	   /// Converting to the degree to radian. This is due to the fact that RP tilt was defined in degrees instead of Rad. 
// 	   /// Note that the RP is the only detector having this issue.
// 	   /// To minimize future confusion and an additional check is in place to make sure RP tile in rad is less than 1.5.
// 	   /// If it is larger, than it have to be in degrees

// 	   if (RP_1_params.get_double_param("rot_y") > 1.5) {
// 	  	float deg_to_rad = RP_1_params.get_double_param("rot_y") * TMath::Pi() / 180.;
// 	   	RP_1_params.set_double_param("rot_y", deg_to_rad); 
// 	   }

// 	   float local_x = Get_Local_X(hit_iter->second->get_x(0), hit_iter->second->get_y(0), hit_iter->second->get_z(0), RP_1_params);
//            float local_y = hit_iter->second->get_y(0);

// //	   cout << local_x << "    " << RP_1_params.get_double_param("rot_y") << "    "  << RP_1_params.get_double_param("rot_y")*TMath::Pi()/180 << endl;
// //	   exit(0);
// 	   RPx[RPhits] = local_x;
// 	   RPy[RPhits] = local_y;
// 	   RPz[RPhits] = hit_iter->second->get_z(0);
// 	   RP_px[RPhits] = hit_iter->second->get_px(0);
// 	   RP_py[RPhits] = hit_iter->second->get_py(0);
// 	   RP_pz[RPhits] = hit_iter->second->get_pz(0);
// 	   RP_edep[RPhits] = hit_iter->second->get_edep();
//            h2_RP_XY_l->Fill(local_x, local_y);
	   
// 	   if(TMath::Abs(RPz[RPhits] - 2600.0)<50) RPind[RPhits ] = 1;
// 	   if(TMath::Abs(RPz[RPhits] - 2800.0)<50) RPind[RPhits ] = 2;
// 	   RPhits++;
// 	   counter++;
//                //---------------------------------------------
//                // Standarized Roman pot cut
//                //
//               // if (local_x > -5. && local_x < 5. && local_y > -1.0 && local_y < 1.0)  {
               
//                	//// This is beam contribution here
               
//                //} else {
               
//                	  ////This is where you signal is!
//                	  //
//                  // if (local_x > -12.5 && local_x < 12.5 && local_y > -5.0 && local_y < 5.0) {
//                    //    h2_RP_XY_signal->Fill(local_x, local_y);
               
//                   //}   
//                }
//             }
// 	 }

//       }
//     }

   return Fun4AllReturnCodes::EVENT_OK;
 }


//***************************************************
// Getting the B0 hits

int diff_tagg_ana::process_B0(PHCompositeNode* topNode)
{
 ostringstream nodename;

  // loop over the G4Hits
  nodename.str("");
  nodename << "G4HIT_" << "b0Truth_0";

  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
//   PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

//   if (!truthinfo)
//     {
//       cout << PHWHERE
// 	   << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles"
// 	   << endl;
//       return Fun4AllReturnCodes::EVENT_OK;
//     }

//   //Get the primary particle range
//   PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();

  Int_t layer=-1;
  if (hits) {
    //this returns an iterator to the beginning and the end of our G4Hits
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
   for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++) {
      //cout << "b0 hit id: " <<  hit_iter->second->get_trkid() << " b0 pz: " << hit_iter->second->get_pz(0) << endl;
      
      h2_B0_XY->Fill(hit_iter->second->get_x(0), hit_iter->second->get_y(0));
      
      //if(hit_iter->second->get_trkid() == 1){
      B0x[B0hits] = (Float_t)hit_iter->second->get_x(0);
      B0y[B0hits] = (Float_t)hit_iter->second->get_y(0);
      B0z[B0hits] = (Float_t)hit_iter->second->get_z(0);
	
      B0xloc[B0hits] = (Float_t)hit_iter->second->get_local_x(0);
      B0yloc[B0hits] = (Float_t)hit_iter->second->get_local_y(0);
      B0zloc[B0hits] = (Float_t)hit_iter->second->get_local_z(0);
	
      B0_px[B0hits] = (Float_t)hit_iter->second->get_px(0);
      B0_py[B0hits] = (Float_t)hit_iter->second->get_py(0);
      B0_pz[B0hits] = (Float_t)hit_iter->second->get_pz(0);
      B0_edep[B0hits] = (Float_t)hit_iter->second->get_edep();
	    
      if(TMath::Abs(B0z[B0hits] - 591.0)<5) {layer = 1;}
      if(TMath::Abs(B0z[B0hits] - 615.0)<5) {layer = 2;}
      if(TMath::Abs(B0z[B0hits] - 639.0)<5) {layer = 3;}
      if(TMath::Abs(B0z[B0hits] - 663.0)<5) {layer = 4;}
	    
      B0ind[B0hits] = layer;
      B0trkid[B0hits] = hit_iter->second->get_trkid();
      
      B0hits++;
      //}//if hit id loop
    }//for hits loop
  }//if hits loop
  return Fun4AllReturnCodes::EVENT_OK;
}
//     for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++) {
//       for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter){

// 	  /// Get this truth particle
// 	  const PHG4Particle *truth = iter->second;

// 	  float m_truthpz = truth->get_pz();
// 	  float m_truthenergy = truth->get_e();
// 	  float m_trutheta = atanh(m_truthpz / m_truthenergy);
// 	  /// Check for nans
// 	  if (m_trutheta != m_trutheta)
// 	    m_trutheta = -99;

// 	  if ( truth->get_pid() == 2212 ){
// 	    if( hit_iter->second->get_hit_id() == 3 || hit_iter->second->get_hit_id() == 4294967297 || hit_iter->second->get_hit_id() == 8589934593 || hit_iter->second->get_hit_id() == 12884901889 ) {

//         PHParameters B0_1_params{"PHB0_1"};

// 	      if (b0_nodeparams){
// 	           B0_1_params.FillFrom(b0_nodeparams, 0);
// 	       }
//         else {
// 	         cerr << "There is a issue finding the detector paramter node!" << endl;
// 	      }
//         float det_x_pos = B0_1_params.get_double_param("place_x") + Enclosure_params.get_double_param("place_x")  + BeamLineMagnet_params.get_double_param("place_x");
//         float det_z_pos = B0_1_params.get_double_param("place_z") + Enclosure_params.get_double_param("place_z")  + BeamLineMagnet_params.get_double_param("place_z");
//         ///	cout << hit_iter->second->get_z(0) << "    " <<  BeamLineMagnet_params.get_double_param("place_z") + Enclosure_params.get_double_param("place_z") << "    " <<  B0_1_params.get_double_param("place_z") << "    " << B0_1_params.get_double_param("length")/(b0DetNr + 1) * (0 - b0DetNr / 2) <<  "    " << BeamLineMagnet_params.get_double_param("length") << "    " << b0DetNr << "    " << z_pos << endl;
//         B0_1_params.set_double_param("place_x", det_x_pos);
//         B0_1_params.set_double_param("place_z", det_z_pos);
//         float local_x = Get_Local_X(hit_iter->second->get_x(0), hit_iter->second->get_y(0), hit_iter->second->get_z(0), B0_1_params);
//         float local_y = hit_iter->second->get_y(0);

// 	      B0x[B0hits] = local_x;
// 	      B0y[B0hits] = local_y;
// 	      B0z[B0hits] = (Float_t)hit_iter->second->get_z(0);

// 	      B0xloc[B0hits] = (Float_t)hit_iter->second->get_local_x(0);
// 	      B0yloc[B0hits] = (Float_t)hit_iter->second->get_local_y(0);
// 	      B0zloc[B0hits] = (Float_t)hit_iter->second->get_local_z(0);

// 	      B0_px[B0hits] = (Float_t)hit_iter->second->get_px(0);
// 	      B0_py[B0hits] = (Float_t)hit_iter->second->get_py(0);
// 	      B0_pz[B0hits] = (Float_t)hit_iter->second->get_pz(0);
// 	      B0_edep[B0hits] = (Float_t)hit_iter->second->get_edep();
// 	      B0truth_px[B0hits] = (Float_t)truth->get_px();
// 	      B0truth_py[B0hits] = (Float_t)truth->get_py();
// 	      B0truth_pz[B0hits] = (Float_t)truth->get_pz();
// 	      B0truth_E[B0hits] = (Float_t)truth->get_e();

// 	      if(TMath::Abs(B0z[B0hits] - 591.0)<5) {layer = 1;}
// 	      if(TMath::Abs(B0z[B0hits] - 615.0)<5) {layer = 2;}
// 	      if(TMath::Abs(B0z[B0hits] - 639.0)<5) {layer = 3;}
// 	      if(TMath::Abs(B0z[B0hits] - 663.0)<5) {layer = 4;}

// 	      //if(TMath::Abs(B0z[B0hits] - 541.0)<50) {layer = 1;}
// 	      //if(TMath::Abs(B0z[B0hits] - 565.0)<50) {layer = 2;}
// 	      //if(TMath::Abs(B0z[B0hits] - 589.0)<50) {layer = 3;}
// 	      //if(TMath::Abs(B0z[B0hits] - 613.0)<50) {layer = 4;}

// 	      B0ind[B0hits ] = layer;
// 	      B0hits++;
// 	    }//if hit_iter
// 	  }//if truth conditions
// 	}//for truth container
//     }//for hit container
//   }//if hits
//   return Fun4AllReturnCodes::EVENT_OK;

// }



//***************************************************

int diff_tagg_ana::process_g4hits_LowQ2Tagger(PHCompositeNode* topNode)
{

   ostringstream nodename;
 
   // loop over the G4Hits
   nodename.str("");
 //  nodename << "G4HIT_" << detector;
 //  nodename << "G4HIT_" << "ZDC";
 //  nodename << "G4HIT_" << "EICG4ZDC";
   nodename << "G4HIT_" << "backLowQ2Tag_1";
 //  nodename << "G4HIT_" << "EEMC";
 
  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());

  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

//  cout << ">>>>>>>>>????????????" << endl;

  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();


  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second;
       ++iter)
  {
      // Get this truth particle
      const PHG4Particle *truth = iter->second;
      if ( truth->get_pid() == 11){ // PDG 11 -> Scattered electron
	e4VectTruth.SetPxPyPzE(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());
 	virtphoton4VectTruth = eBeam4Vect - e4VectTruth;
        Q2_truth = -1*(virtphoton4VectTruth.Mag2());
  
	h_Q2_truth->Fill(Q2_truth);
	h_Q2->Fill(Q2_truth);

	h_log_Q2->Fill(Q2_truth);



	e_E_truth = e4VectTruth.E();
	e_eta_truth = e4VectTruth.Eta();
        e_Phi_truth = e4VectTruth.Phi();

	e_Polar_truth = (TMath::Pi() - e4VectTruth.Theta()) * 1000.;


	h_E->Fill(e_E_truth);
	h_eta->Fill(e_eta_truth);
        h_polar->Fill(e_Polar_truth);
	h2_E_Q2->Fill(e_E_truth, Q2_truth);



//        cout  << "Q2: "<<  Q2_truth << endl;
//
	
  
      }
   
  }

  if (!truthinfo)
  {
    cout << PHWHERE
         << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles"
         << endl;

    return Fun4AllReturnCodes::EVENT_OK;
  }


//  /// Get the primary particle range
//  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();


  if (hits) {
//    // this returns an iterator to the beginning and the end of our G4Hits
    PHG4HitContainer::ConstRange hit_range = hits->getHits();

    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++) {

	cout << "Low Q2 Tagger hits? " ; 
	cout << "This is where you can fill your loop " << endl;

	
//	PHG4Particle* g4particle = truthinfo->GetParticle(hit_iter->second->get_trkid());

//	cout << Q2_truth << endl;
//      exit(0);

	////************************************************************************
	//// From hit to particle

	PHG4Particle* g4particle_hit = truthinfo->GetParticle(hit_iter->second->get_trkid());

//	cout << g4particle_hit->get_px() << "    " << hit_iter->second->get_px(0) << endl;
//        exit(0);

	float r = sqrt(hit_iter->second->get_x(0) * hit_iter->second->get_x(0) + hit_iter->second->get_y(0)* hit_iter->second->get_y(0));

	float theta = fabs(atan(r/hit_iter->second->get_z(0)));

	cout << theta << endl;


	if (Q2_truth < 0.1 && e_eta_truth < -4.5 && e_E_truth < 17.7) {

//	h_Q2_truth->Fill(Q2_truth);
	h2_Q2_pos->Fill(hit_iter->second->get_x(0), Q2_truth);

	h2_lowQ2_XY->Fill(hit_iter->second->get_x(0), hit_iter->second->get_y(0));

//	h2_Q2_theta->Fill(theta, Q2_truth);

	h_Q2_truth_LowQ2tag->Fill(Q2_truth);

	h_log_Q2_LowQ2tag->Fill(Q2_truth);

//	cout << hit_iter->second->get_edep() << endl;
	h2_Q2_truth_E->Fill(Q2_truth,  g4particle_hit->get_e());


//	h_E->Fill(e_E_truth);
//	h_eta->Fill(e_eta_truth);
//      h_polar->Fill(e_Phi_truth);
//	h2_E_Q2->Fill(e_E_truth, Q2_truth);

	h_E_LowQ2tag->Fill(e_E_truth);
	h_eta_LowQ2tag->Fill(e_eta_truth);
        h_polar_LowQ2tag->Fill(e_Polar_truth);
	h2_E_Q2_LowQ2tag->Fill(e_E_truth, Q2_truth);

//	cout << r << "  "<< Q2_truth << endl;

//	exit(0);

	}


    }
  
  }



  return Fun4AllReturnCodes::EVENT_OK;

}

int diff_tagg_ana::process_ClusterCalo(PHCompositeNode* topNode, string caloName)
{

  if(caloName == "CLUSTER_FEMC"){
    CaloRawClusterEval* clusterevalFEMC = _caloevalstackFEMC->get_rawcluster_eval();
    RawClusterContainer* clusters = findNode::getClass<RawClusterContainer>(topNode, caloName.c_str());
    if(clusters){
      RawClusterContainer::ConstRange cluster_range = clusters->getClusters();
      for(RawClusterContainer::ConstIterator cluster_iter = cluster_range.first; cluster_iter != cluster_range.second; cluster_iter++){
	if(!cluster_iter->second) continue;
	float clX = cluster_iter->second->get_x();
	float clY = cluster_iter->second->get_y();
	float clZ = cluster_iter->second->get_z();
	//float clPX = cluster_iter->second->get_px();
	//float clPY = cluster_iter->second->get_py();
	//float clPZ = cluster_iter->second->get_pz();
	float clEn = cluster_iter->second->get_energy();

	hitX[nHits] = clX;
	hitY[nHits] = clY;
	hitZ[nHits] = clZ;
	//hitPX[nHits] = clPX;
	//hitPY[nHits] = clPY;
	//hitPZ[nHits] = clPZ;
	hitE[nHits] = clEn;
	caloInd[nHits] = 1;
	hitsNtowers[nHits] = cluster_iter->second->getNTowers();
	if(clusterevalFEMC){
	  PHG4Particle* primary = clusterevalFEMC->max_truth_primary_particle_by_energy(cluster_iter->second);
	  if(primary){
	    hitPid[nHits] = primary->get_pid();
	  }
	}
	nHits++;
	hitsFEMC++;
      } // for loop
    } // clusters

  } // FEMC

  if(caloName == "CLUSTER_EEMC"){
    CaloRawClusterEval* clusterevalEEMC = _caloevalstackEEMC->get_rawcluster_eval();
    RawClusterContainer* clusters = findNode::getClass<RawClusterContainer>(topNode, caloName.c_str());
    if(clusters){
      RawClusterContainer::ConstRange cluster_range = clusters->getClusters();
      for(RawClusterContainer::ConstIterator cluster_iter = cluster_range.first; cluster_iter != cluster_range.second; cluster_iter++){
	if(!cluster_iter->second) continue;
	float clX = cluster_iter->second->get_x();
	float clY = cluster_iter->second->get_y();
	float clZ = cluster_iter->second->get_z();
	//float clPX = cluster_iter->second->get_px();
	//float clPY = cluster_iter->second->get_py();
	//float clPZ = cluster_iter->second->get_pz();
	float clEn = cluster_iter->second->get_energy();

	hitX[nHits] = clX;
	hitY[nHits] = clY;
	hitZ[nHits] = clZ;
	//hitPX[nHits] = clPX;
	//hitPY[nHits] = clPY;
	//hitPZ[nHits] = clPZ;
	hitE[nHits] = clEn;
	caloInd[nHits] = 2;
	hitsNtowers[nHits] = cluster_iter->second->getNTowers();

	if(clusterevalEEMC){
	  PHG4Particle* primary = clusterevalEEMC->max_truth_primary_particle_by_energy(cluster_iter->second);
	  if(primary){
	    hitPid[nHits] = primary->get_pid();
	  }
	}
	nHits++;
	hitsEEMC++;
      } // for loop
    } // clusters

  } // EEMC



  if(caloName == "CLUSTER_BECAL"){
    CaloRawClusterEval* clusterevalBECAL = _caloevalstackBECAL->get_rawcluster_eval();
    RawClusterContainer* clusters = findNode::getClass<RawClusterContainer>(topNode, caloName.c_str());
    if(clusters){
      RawClusterContainer::ConstRange cluster_range = clusters->getClusters();
      for(RawClusterContainer::ConstIterator cluster_iter = cluster_range.first; cluster_iter != cluster_range.second; cluster_iter++){
	if(!(cluster_iter->second)) continue;
	float clX = cluster_iter->second->get_x();
	float clY = cluster_iter->second->get_y();
	float clZ = cluster_iter->second->get_z();
	//float clPX = cluster_iter->second->get_px();
	//float clPY = cluster_iter->second->get_py();
	//float clPZ = cluster_iter->second->get_pz();
	float clEn = cluster_iter->second->get_energy();

	hitX[nHits] = clX;
	hitY[nHits] = clY;
	hitZ[nHits] = clZ;
	//hitPX[nHits] = clPX;
	//hitPY[nHits] = clPY;
	//hitPZ[nHits] = clPZ;
	hitE[nHits] = clEn;
	caloInd[nHits] = 3;
	hitsNtowers[nHits] = cluster_iter->second->getNTowers();

	if(clusterevalBECAL){
	  PHG4Particle* primary = clusterevalBECAL->max_truth_primary_particle_by_energy(cluster_iter->second);

	  //int val;

	  if(primary)
	    {
	      hitPid[nHits] = primary->get_pid();
	    }
	}
	nHits++;
	hitsBECAL++;
      } // for loop

    } // clusters

  } // EEMC

  return Fun4AllReturnCodes::EVENT_OK;
}

///*****************************************************
/// ZDC Energy and Poisition smearing functions
//
// Energy smearing

float diff_tagg_ana::ZDC_Energy_Smear_EMCAL(float E) {

  float resolution, E_reco;

  resolution = sqrt(.45*.45/E + 0.075*0.075);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;
}

// Energy smearing

float diff_tagg_ana::ZDC_Energy_Smear_HCAL(float E) {

  float resolution, E_reco;

  resolution = sqrt(.50*.50/E + 0.1*0.1);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;
}

// Energy smearing

float diff_tagg_ana::ZDC_Energy_Smear_PbWO4(float E) {

  float resolution, E_reco;

  resolution = sqrt(.25*.25/E + 0.04*0.04);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;

}

// Posision smearing

float diff_tagg_ana::ZDC_Position_Smear(float P) {

  float resolution, P_reco;

  resolution = 0.1;         /// Position resolution 0.1 cm
  P_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * P;

  return P_reco;

}

///*****************************************************
/// B0 tracker smearing functions

// Energy smearing

float diff_tagg_ana::B0Tracker_Energy_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.25*.25/E + 0.04*0.04);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;

}

// Posision smearing

float diff_tagg_ana::B0Tracker_Position_Smear(float P) {

  float resolution, P_reco;

  resolution = 0.1;         /// Position resolution 0.1 cm
  P_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * P;

  return P_reco;

}


///*****************************************************
/// B0 Cal smearing functions

// Energy smearing

float diff_tagg_ana::B0Cal_Energy_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.25*.25/E + 0.04*0.04);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;

}

// Posision smearing

float diff_tagg_ana::B0Cal_Position_Smear(float P) {

  float resolution, P_reco;

  resolution = 0.1;         /// Position resolution 0.1 cm
  P_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * P;

  return P_reco;

}

///*****************************************************
/// RP smearing functions

float diff_tagg_ana::RP_Energy_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.25*.25/E + 0.04*0.04);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;

}

// Posision smearing

float diff_tagg_ana::RP_Position_Smear(float P) {

  float resolution, P_reco;

  resolution = 0.1;         /// Position resolution 0.1 cm
  P_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * P;

  return P_reco;

}

///*****************************************************
/// Off momentum smearing functions

float diff_tagg_ana::Off_Mom_Energy_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.25*.25/E + 0.04*0.04);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;

}

// Posision smearing

float diff_tagg_ana::Off_Mom_Position_Smear(float P) {

  float resolution, P_reco;

  resolution = 0.1;         /// Position resolution 0.1 cm
  P_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * P;

  return P_reco;

}

int diff_tagg_ana::process_g4hits(PHCompositeNode* topNode, const std::string& detector)
{
  ostringstream nodename;

  // loop over the G4Hits
  nodename.str("");
  nodename << "G4HIT_" << detector;
  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
  if (hits)
  {
    // this returns an iterator to the beginning and the end of our G4Hits
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)

    {
      // the pointer to the G4Hit is hit_iter->second
      g4hitntuple->Fill(hit_iter->second->get_x(0),
                        hit_iter->second->get_y(0),
                        hit_iter->second->get_z(0),
                        hit_iter->second->get_x(1),
                        hit_iter->second->get_y(1),
                        hit_iter->second->get_z(1),
                        hit_iter->second->get_edep());
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//*******************************************

int diff_tagg_ana::process_g4clusters(PHCompositeNode* topNode, const string& detector)
{
  ostringstream nodename;

  // loop over the G4Hits
  nodename.str("");
  nodename << "CLUSTER_" << detector;
  RawClusterContainer* clusters = findNode::getClass<RawClusterContainer>(topNode, nodename.str().c_str());
  if (clusters)
  {
    RawClusterContainer::ConstRange cluster_range = clusters->getClusters();
    for (RawClusterContainer::ConstIterator cluster_iter = cluster_range.first; cluster_iter != cluster_range.second; cluster_iter++)
    {
      clusterntuple->Fill(cluster_iter->second->get_phi(),
                          cluster_iter->second->get_z(),
                          cluster_iter->second->get_energy(),
                          cluster_iter->second->getNTowers());
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//*******************************************
//Local Coords Functions
float diff_tagg_ana::Get_Local_Y(float global_x, float global_y, float global_z, float det_tilt, float cross_angle) {

	return global_y;

}


float diff_tagg_ana::Get_Local_X(float global_x, float global_y, float global_z, float det_tilt, float det_rot) {

  TVector3 global_cor(global_x, global_y, global_z);
  float local_x;

  global_cor.RotateY(-det_rot);
  local_x = global_cor.X()/cos(det_tilt - det_rot);
	
  return local_x;

}


float diff_tagg_ana::Get_Local_X(float global_x, float global_y, float global_z, PdbParameterMapContainer *det_nodeparams) {
  
  PHParameters Det_params{"PHDet"};
   
  if (det_nodeparams)
    {
      Det_params.FillFrom(det_nodeparams, 0);
    } else {
    cerr << "There is a issue finding the detector paramter node!" << endl;
  }
  
  float det_xCent = Enclosure_params.get_double_param("place_x") + Det_params.get_double_param("Layer1_pos_x");
  float det_zCent = Enclosure_params.get_double_param("place_z") + Det_params.get_double_param("Layer1_pos_z");
  float det_tilt = Det_params.get_double_param("Layer1_rot_y")/180 * TMath::Pi();

  float det_rot = atan( det_xCent / det_zCent);

  TVector3 global_cor(global_x, global_y, global_z);
  float local_x;

  global_cor.RotateY(-det_rot);
  
  local_x = global_cor.X()/cos(det_tilt);
  
  return local_x;
}


float diff_tagg_ana::Get_Local_X(float global_x, float global_y, float global_z, PHParameters Det_params) {

   float det_xCent = Det_params.get_double_param("Layer1_pos_x");
   float det_zCent = Det_params.get_double_param("Layer1_pos_z");

//   float det_tilt = Det_params.get_double_param("rot_y")/180. * TMath::Pi(); // in Rad
   float det_tilt = Det_params.get_double_param("Layer1_rot_y"); // in Rad

   float det_rot = atan( det_xCent / det_zCent);  // in Rad

   TVector3 global_cor(global_x, global_y, global_z);

//   float local_x;
//   global_cor.RotateY(-det_rot);
//   local_x = global_cor.X()/cos(det_tilt);


   float local_x1 = Get_Local_X(global_x, global_y, global_z, det_tilt, det_rot);
//   cout << local_x1 << "    "<< local_x << endl;
//   exit(0);

   return local_x1;

}



void diff_tagg_ana::initializeTrees()
{
  tree = new TTree("T", "A tree for TCS");

  // tree->Branch("local_x", &local_x, "local_x/F");
  //tree->Branch("local_y", &local_y, "local_y/F");

  tree->Branch("Epx", &Epx, "Epx/F");
  tree->Branch("Epy", &Epy, "Epy/F");
  tree->Branch("Epz", &Epz, "Epz/F");

  tree->Branch("Ppx", &Ppx, "Ppx/F");
  tree->Branch("Ppy", &Ppy, "Ppy/F");
  tree->Branch("Ppz", &Ppz, "Ppz/F");

  tree->Branch("Gpx", &Gpx, "Gpx/F");
  tree->Branch("Gpy", &Gpy, "Gpy/F");
  tree->Branch("Gpz", &Gpz, "Gpz/F");

  tree->Branch("Pos_px", &Pos_px, "Pos_px/F");
  tree->Branch("Pos_py", &Pos_py, "Pos_py/F");
  tree->Branch("Pos_pz", &Pos_pz, "Pos_pz/F");

  tree->Branch("nHits",&nHits,"nHits/I");
  tree->Branch("caloInd",&caloInd,"caloInd[nHits]/I");
  tree->Branch("hitX",&hitX,"hitX[nHits]/F");
  tree->Branch("hitY",&hitY,"hitY[nHits]/F");
  tree->Branch("hitZ",&hitZ,"hitZ[nHits]/F");
  //tree->Branch("hitPX",&hitPX,"hitPX[nHits]/F");
  //tree->Branch("hitPY",&hitPY,"hitPY[nHits]/F");
  //tree->Branch("hitPZ",&hitPZ,"hitPZ[nHits]/F");
  tree->Branch("hitE",&hitE,"hitE[nHits]/F");
  tree->Branch("hitPid",&hitPid,"hitPid[nHits]/I");
  tree->Branch("hitsNtowers",&hitsNtowers,"hitsNtowers[nHits]/I");
	
  tree->Branch("RPhits",&RPhits,"RPhits/I");
  tree->Branch("RPind",&RPind,"RPind[RPhits]/I");
  tree->Branch("RPx",&RPx,"RPx[RPhits]/F");
  tree->Branch("RPy",&RPy,"RPy[RPhits]/F");
  tree->Branch("RPz",&RPz,"RPz[RPhits]/F");
  tree->Branch("RPxloc",&RPxloc,"RPxloc[RPhits]/F");
  tree->Branch("RPyloc",&RPyloc,"RPyloc[RPhits]/F");
  tree->Branch("RPzloc",&RPzloc,"RPzloc[RPhits]/F");
  tree->Branch("RP_px",&RP_px,"RP_px[RPhits]/F");
  tree->Branch("RP_py",&RP_py,"RP_py[RPhits]/F");
  tree->Branch("RP_pz",&RP_pz,"RP_pz[RPhits]/F");
  tree->Branch("RP_edep",&RP_edep,"RP_edep[RPhits]/F");
  tree->Branch("RPtrkid",&RPtrkid,"RPtrkid[RPhits]/F");
	
  tree->Branch("B0hits",&B0hits,"B0hits/I");
  tree->Branch("B0x",&B0x,"B0Px[B0hits]/F");
  tree->Branch("B0y",&B0y,"B0y[B0hits]/F");
  tree->Branch("B0z",&B0z,"B0z[B0hits]/F");
  tree->Branch("B0xloc",&B0xloc,"B0xloc[B0hits]/F");
  tree->Branch("B0yloc",&B0yloc,"B0yloc[B0hits]/F");
  tree->Branch("B0zloc",&B0zloc,"B0zloc[B0hits]/F");
  tree->Branch("B0_px",&B0_px,"B0_px[B0hits]/F");
  tree->Branch("B0_py",&B0_py,"B0_py[B0hits]/F");
  tree->Branch("B0_pz",&B0_pz,"B0_pz[B0hits]/F");
  tree->Branch("B0_edep",&B0_edep,"B0_edep[B0hits]/F");
  tree->Branch("B0ind",&B0ind,"B0ind[B0hits]/I");
  tree->Branch("B0trkid",&B0trkid,"B0trkid[B0hits]/F");


  tree->Branch("ntr",&ntr,"ntr/I");
  tree->Branch("tr_px",&tr_px,"tr_px[ntr]/F");
  tree->Branch("tr_py",&tr_py,"tr_py[ntr]/F");
  tree->Branch("tr_pz",&tr_pz,"tr_pz[ntr]/F");
  tree->Branch("tr_p",&tr_p,"tr_p[ntr]/F");

  tree->Branch("tr_x",&tr_x,"tr_x[ntr]/F");
  tree->Branch("tr_y",&tr_y,"tr_y[ntr]/F");
  tree->Branch("tr_z",&tr_z,"tr_z[ntr]/F");

  tree->Branch("tr_phi",&tr_phi,"tr_phi[ntr]/F");
  tree->Branch("tr_eta",&tr_eta,"tr_eta[ntr]/F");
  tree->Branch("tr_Pid",&tr_Pid,"tr_Pid[ntr]/I");
  tree->Branch("charge",&charge,"charge[ntr]/F");


  tree->Branch("hitsEEMC",&hitsEEMC,"hitsEEMC/I");
  tree->Branch("hitsFEMC",&hitsFEMC,"hitsFEMC/I");
  tree->Branch("hitsBECAL",&hitsBECAL,"hitsBECAL/I");
}


void diff_tagg_ana::initializeVariables()
{
  Epx= -1000;
  Epy= -1000;
  Epz= -1000;

  Ppx= -1000;
  Ppy= -1000;
  Ppz= -1000;
  Gpx= -1000;
  Gpy= -1000;
  Gpz= -1000;
  Pos_px= -1000;
  Pos_py= -1000;
  Pos_pz= -1000;

  nHits=0;
  ntr=0;
  RPhits=0;
  B0hits=0;
  RP1=0;
  RP2=0;
  BRP1=0;
  BRP2=0;

  hitsEEMC=0;
  hitsFEMC=0;
  hitsBECAL=0;

  for(int i=0;i<100;i++){
    hitsNtowers[i]=0;
    caloInd[i]=-1;
    hitX[i]=-1000;
    hitY[i]=-1000;
    hitZ[i]=-1000;
    //hitPX[i]=-1000;
    //hitPY[i]=-1000;
    //hitPZ[i]=-1000;
    hitE[i]=-1000;
    hitPid[i]=0;

    tr_px[i]=-1000;
    tr_py[i]=-1000;
    tr_pz[i]=-1000;
    tr_p[i]=-1000;
    tr_phi[i]=-1000;
    tr_eta[i]=-1000;
    charge[i]=-1000;
    tr_x[i]=-1000;
    tr_y[i]=-1000;
    tr_z[i]=-1000;
    tr_Pid[i]=0;

    RPx[i]=-1000;
    RPy[i]=-1000;
    RPz[i]=-1000;
    RP_px[i]=-1000;
    RP_py[i]=-1000;
    RP_pz[i]=-1000;
//     RPtruth_px[RPhits]=-1000;
//     RPtruth_py[RPhits]=-1000;
//     RPtruth_pz[RPhits]=-1000;
//     RPtruth_E[i]=-1000;
    RPind[i]=-1;

//     BRPx[i]=-1000;
//     BRPy[i]=-1000;
//     BRPz[i]=-1000;
//     BRP_px[i]=-1000;
//     BRP_py[i]=-1000;
//     BRP_pz[i]=-1000;
//     BRPtruth_px[RPhits]=-1000;
//     BRPtruth_py[RPhits]=-1000;
//     BRPtruth_pz[RPhits]=-1000;
//     BRPtruth_E[i]=-1000;
//     BRPind[i]=-1;



    B0x[i]=-1000;
    B0y[i]=-1000;
    B0z[i]=-1000;
    B0xloc[i]=-1000;
    B0yloc[i]=-1000;
    B0zloc[i]=-1000;
    B0_px[i]=-1000;
    B0_py[i]=-1000;
    B0_pz[i]=-1000;
    B0_edep[i]=-1000;
//     B0truth_px[i]=-1000;
//     B0truth_py[i]=-1000;
//     B0truth_pz[i]=-1000;
//     B0truth_E[i]=-1000;
    B0ind[i]=-1;
  }
}


int diff_tagg_ana::process_tracks(PHCompositeNode *topNode)
{
  SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");

  if (!trackmap)
    return Fun4AllReturnCodes::EVENT_OK;

 // EvalStack for truth track matching
    if(!_svtxEvalStack)
        {
          _svtxEvalStack = new SvtxEvalStack(topNode);
           _svtxEvalStack->set_verbosity(Verbosity());
         }
    else{
           _svtxEvalStack->next_event(topNode);
         }
  
    SvtxTrackEval *trackeval = _svtxEvalStack->get_track_eval();
//    PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  
    for (SvtxTrackMap::Iter iter = trackmap->begin();iter != trackmap->end();++iter)
        {
         SvtxTrack *track = iter->second;
  
         // Get the reconstructed track info
         tr_px[ntr] = track->get_px();
         tr_py[ntr] = track->get_py();
         tr_pz[ntr] = track->get_pz();
         tr_p[ntr] = TMath::Sqrt(tr_px[ntr] * tr_px[ntr] + tr_py[ntr] * tr_py[ntr] + tr_pz[ntr] * tr_pz[ntr]);
  
         tr_phi[ntr] = track->get_phi();
         tr_eta[ntr] = track->get_eta();
  
         charge[ntr] = track->get_charge();
         tr_x[ntr] = track->get_x();
         tr_y[ntr] = track->get_y();
         tr_z[ntr] = track->get_z();
  
         // Get truth track info that matches this reconstructed track
         PHG4Particle *truthtrack = trackeval->max_truth_particle_by_nclusters(track);
         if(truthtrack) tr_Pid[ntr] = truthtrack->get_pid();
         ntr++;
        }
  return Fun4AllReturnCodes::EVENT_OK;
}
  

