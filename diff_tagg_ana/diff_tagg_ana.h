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
class TNtuple;


class diff_tagg_ana : public SubsysReco
{
 public:

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

  int process_g4hits_ZDC(PHCompositeNode *);
  int process_g4hits_RomanPots(PHCompositeNode *);
  int process_g4hits_B0(PHCompositeNode *);

  int process_g4hits_LowQ2Tagger(PHCompositeNode *);


  int process_g4hits(PHCompositeNode *, const std::string&);
  int process_g4clusters(PHCompositeNode *, const std::string&);

  int process_PHG4Truth(PHCompositeNode* topNode);
  int process_PHG4Truth_Primary_Particles(PHCompositeNode* topNode);





// private:

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

//  double e_beam_energy;
//  double ion_beam_energy;
//  double crossing_angle;

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
  TH2F* h2_B0_XY_g; 
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

  TH1F* h_log_Q2; 
  TH1F* h_log_Q2_LowQ2tag;

  // Beam parameter

  Float_t e_beam_energy;
  Float_t e_beam_pmag;

  Float_t ion_beam_energy;
  Float_t ion_beam_pmag;

  Float_t crossing_angle;

  Float_t Q2_truth;
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

  TString IP_design;

};

#endif // DIFF_TAGG_ANA_H
