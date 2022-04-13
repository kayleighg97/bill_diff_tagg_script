// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef DIFF_TAGG_ANA_H
#define DIFF_TAGG_ANA_H

#include <fun4all/SubsysReco.h>

#include <fun4all/Fun4AllServer.h>
#include <g4main/PHG4Reco.h>

#include <string>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"

#include <pdbcalbase/PdbParameterMap.h>
#include <phparameter/PHParameters.h>

class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TTree;
class TNtuple;
class CaloEvalStack;
class CaloRawClusterEval;
class RawClusterContainer;
class SvtxTrackMap;
class SvtxEvalStack;
class SvtxTrackEval;
class PHG4TruthInfoContainer;


class diff_tagg_ana : public SubsysReco
{
 public:
  CaloEvalStack *_caloevalstack;
  CaloEvalStack* _caloevalstackFEMC;
  CaloEvalStack* _caloevalstackEEMC;
  CaloEvalStack* _caloevalstackBECAL;

  SvtxEvalStack *_svtxEvalStack;
//  diff_tagg_ana(const std::string &name = "diff_tagg_ana");
  diff_tagg_ana(const std::string &name = "Diff_Tagg_ana", const std::string &fname = "MyNtuple.root");

  virtual ~diff_tagg_ana();

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  int process_ZDC(PHCompositeNode *);
  int process_RomanPots(PHCompositeNode *);//, int);
  int process_B0(PHCompositeNode *);
  int process_ClusterCalo(PHCompositeNode*, std::string);
  int process_tracks(PHCompositeNode *);
  int process_Q2tagger(PHCompositeNode *);
  
  int process_g4hits_ZDC(PHCompositeNode *);
  int process_g4hits_RomanPots(PHCompositeNode *);
  int process_g4hits_B0(PHCompositeNode *);

  int process_g4hits_LowQ2Tagger(PHCompositeNode *);
  int process_g4hits(PHCompositeNode* topNode, const std::string&);
  int process_g4clusters(PHCompositeNode *, const std::string&);

  int process_PHG4Truth(PHCompositeNode* topNode);
  int process_PHG4Truth_Primary_Particles(PHCompositeNode* topNode);
  
 private:
 
 TTree *tree;

  float Epx;
  float Epy;
  float Epz;

  float Ppx;
  float Ppy;
  float Ppz;

  float Gpx;
  float Gpy;
  float Gpz;

  float Pos_px;
  float Pos_py;
  float Pos_pz;

  int RP1;
  int RP2;
  int BRP1;
  int BRP2;
  int RPhits;
    
  float RPx[100];
  float RPy[100];
  float RPz[100];
  Float_t RPxloc[100];
  Float_t RPyloc[100];
  Float_t RPzloc[100];
  float RP_px[100];
  float RP_py[100];
  float RP_pz[100];
  float RP_edep[100];
  int RPind[100];
  float RPtrkid[100];
  float RPtruth_px[100]
  float RPtruth_py[100]
  float RPtruth_pz[100]
  float RPtruth_E[100]
  
  Int_t B0hits;
    
  Float_t B0x[100];
  Float_t B0y[100];
  Float_t B0z[100];
  Float_t B0xloc[100];
  Float_t B0yloc[100];
  Float_t B0zloc[100];
  float B0_px[100];
  float B0_py[100];
  float B0_pz[100];
  float B0_edep[100];
  Int_t B0ind[100];
  float B0trkid[100];
    
  int nHits;
  int caloInd[10000];
  float hitX[10000];
  float hitY[10000];
  float hitZ[10000];
  float hitE[10000];
  int hitPid[10000];
  int hitsNtowers[10000];
  int hitsEEMC;
  int hitsFEMC;
  int hitsBECAL;

  int ntr;
  float tr_px[10000];
  float tr_py[10000];
  float tr_pz[10000];
  float tr_p[10000];
  float tr_phi[10000];
  float tr_eta[10000];
  float charge[10000];
  float tr_x[10000];
  float tr_y[10000];
  float tr_z[10000];
  int tr_Pid[10000];

  void initializeVariables();
  void initializeTrees();
  
 protected:

  std::string detector;
  std::string outfilename;
  Fun4AllHistoManager *hm;

  TFile *outfile;
  TNtuple *g4hitntuple;
  TNtuple *clusterntuple;

  unsigned long long int event_itt;
  gsl_rng* m_RandomGenerator;

  int static_event_counter;

  //*********************************
  // ZDC Energy and Position smearing

  float ZDC_Energy_Smear_EMCAL(float E);
  float ZDC_Energy_Smear_HCAL(float E);
  float ZDC_Energy_Smear_PbWO4(float E);
  float ZDC_Position_Smear(float E);

  //*********************************
  // B0 Tracker Energy and Position smearing

  float B0Tracker_Energy_Smear(float E);
  float B0Tracker_Position_Smear(float E);

  //*********************************
  // B0 Cal Energy and Position smearing

  float B0Cal_Energy_Smear(float E);
  float B0Cal_Position_Smear(float E);

  //*********************************
  // RP Energy and Position smearing

  float RP_Energy_Smear(float E);
  float RP_Position_Smear(float E);

  //*********************************
  // Off Momentum Energy and Position smearing

  float Off_Mom_Energy_Smear(float E);
  float Off_Mom_Position_Smear(float E);

  //*********************************
  // Coordinate transformation from global to local

  float Get_Local_X(float global_x, float global_y, float global_z, float det_tilt, float det_rot);
  float Get_Local_Y(float global_x, float global_y, float global_z, float det_tilt, float det_rot);
  float Get_Local_X(float global_x, float global_y, float global_z, PdbParameterMapContainer *det_nodeparams);
//  float Get_Local_X(float global_x, float global_y, float global_z) {return 1;};
  float Get_Local_X(float global_x, float global_y, float global_z, PHParameters Det_params);
  //---------------------
    //---------------------
  // From ejana

  double true_q2;
  double true_x;
  double true_s_e;
  double true_xpi;
  double true_ypi;
  double true_tpi;

  double have_true_dis_info = false;
  
  bool  HIT_IN_ZDC; 
  bool  HIT_IN_HEC;	


  int b0DetNr;

  TLorentzVector r_lelectron;
//  TLorentzVector r_lproton;
//  TLorentzVector r_lproton;

  TLorentzVector r_lscatelec;
  TLorentzVector r_l_scat_nucleon;

  TLorentzVector lproton;

  Int_t ZDC_hit;

  TH2F* h2_ZDC_XY_g; 
  TH2F* h2_ZDC_XY_g_double; 

  TH2F* h2_ZDC_XY_l; 
  TH2F* h2_ZDC_XY_l_double; 

  TH1F* h1_E_dep_smeared;
  TH1F* h1_E_dep;

  // Roman pots
  TH2F* h2_RP_XY_g; 
  TH2F* h2_RP_XY_l; 
  TH2F* h2_RP_XY_signal; 

  // B0
  TH2F* h2_B0_XY; 
  TH2F* h2_B0_XY_l; 

  // Low Q2 Tagger
  TH2F* h2_lowQ2_XY; 
  TH1F* h_Q2_truth; 
  TH1F* h_Q2_truth_LowQ2tag; 
  TH2F* h2_Q2_pos; 
  TH2F* h2_Q2_mom; 
  TH2F* h2_pos_mom; 
  TH2F* h2_Q2_theta; 

  TH2F* h2_Q2_truth_E; 

  TH1F* h_Q2; 
  TH1F* h_log_Q2; 
  TH1F* h_log_Q2_LowQ2tag;

  TH1F* h_E;
  TH1F* h_E_LowQ2tag;

  TH1F* h_eta;
  TH1F* h_eta_LowQ2tag;

  TH1F* h_polar;
  TH1F* h_polar_LowQ2tag; 

  TH2F* h2_E_Q2;
  TH2F* h2_E_Q2_LowQ2tag;


  // Beam parameter

  Float_t e_beam_energy;
  Float_t e_beam_pmag;

  Float_t ion_beam_energy;
  Float_t ion_beam_pmag;

  Float_t crossing_angle;

  Float_t Q2_truth;
  Float_t e_eta_truth;
  Float_t e_Phi_truth;
  Float_t e_E_truth;
  Float_t e_Polar_truth;


  Float_t mProt;
  Float_t mElec;

  TLorentzVector eBeam4Vect;
  TLorentzVector pBeam4Vect;
  TLorentzVector virtphoton4VectTruth;
  TLorentzVector e4VectTruth;

  //-------------------------------
  int m_mpi;
  int m_process_id;
  double m_truthenergy;
  double m_trutheta;
  double m_truthphi;
  double m_truthpx;
  double m_truthpy;
  double m_truthpz;
  double m_truthpt;
  double m_truthp;
  int m_numparticlesinevent;
  int m_truthpid;
  
  PHParameters Enclosure_params{"PHGEnclosure"};
  PHParameters ZDC_params{"PHG4RP"};
  PHParameters RP_1_params{"PHG4RP"};
  PHParameters RP2_params{"PHG4RP2"};
  PHParameters B0_params{"PHG4B0"};
  PHParameters BeamLineMagnet_params{"PHG4BeamLinMagnet"};
  
  PdbParameterMapContainer *encloseure_nodeparams;
  PdbParameterMapContainer *zdc_nodeparams;
  PdbParameterMapContainer *rp_nodeparams;
  PdbParameterMapContainer *rp2_nodeparams;
  PdbParameterMapContainer *b0_nodeparams;
  PdbParameterMapContainer *beamlinemagnet_nodeparams;
  // From ejana


  TString IP_design;

};

#endif // DIFF_TAGG_ANA_H
