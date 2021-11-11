#ifndef MACRO_G4TRACKINGEIC_C
#define MACRO_G4TRACKINGEIC_C

#include <GlobalVariables.C>

#include <g4trackfastsim/PHG4TrackFastSim.h>

#include <trackreco/PHRaveVertexing.h>

#include <g4trackfastsim/PHG4TrackFastSimEval.h>

#include <fun4all/Fun4AllServer.h>

#include <vector>

R__LOAD_LIBRARY(libtrack_reco.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

namespace Enable
{
  bool TRACKING = false;
  bool TRACKING_EVAL = false;
  int TRACKING_VERBOSITY = 0;
}  // namespace Enable

namespace G4TRACKING
{
  bool DISPLACED_VERTEX = false;
  bool PROJECTION_EEMC = false;
  bool PROJECTION_EHCAL = false;
  bool PROJECTION_CEMC = false;
  bool PROJECTION_BECAL = false;
  bool PROJECTION_HCALIN = false;
  bool PROJECTION_HCALOUT = false;
  bool PROJECTION_FEMC = false;
  bool PROJECTION_FHCAL = false;
  bool PROJECTION_LFHCAL = false;
}  // namespace G4TRACKING

//-----------------------------------------------------------------------------//
void TrackingInit()
{
  TRACKING::FastKalmanFilter = new PHG4TrackFastSim("PHG4TrackFastSim");
  TRACKING::FastKalmanFilterSiliconTrack = new PHG4TrackFastSim("FastKalmanFilterSiliconTrack");
  TRACKING::FastKalmanFilterInnerTrack = new PHG4TrackFastSim("FastKalmanFilterInnerTrack");
}

void InitFastKalmanFilter(PHG4TrackFastSim *kalman_filter)
{
  //  kalman_filter->Smearing(false);
  if (G4TRACKING::DISPLACED_VERTEX)
  {
    // do not use truth vertex in the track fitting,
    // which would lead to worse momentum resolution for prompt tracks
    // but this allows displaced track analysis including DCA and vertex finding
    kalman_filter->set_use_vertex_in_fitting(false);
    kalman_filter->set_vertex_xy_resolution(0);  // do not smear the vertex used in the built-in DCA calculation
    kalman_filter->set_vertex_z_resolution(0);   // do not smear the vertex used in the built-in DCA calculation
    kalman_filter->enable_vertexing(true);       // enable vertex finding and fitting
  }
  else
  {
    // constraint to a primary vertex and use it as part of the fitting level arm
    kalman_filter->set_use_vertex_in_fitting(true);
    kalman_filter->set_vertex_xy_resolution(50e-4);
    kalman_filter->set_vertex_z_resolution(50e-4);
  }

  kalman_filter->set_sub_top_node_name("TRACKS");
}

//-----------------------------------------------------------------------------//
void Tracking_Reco()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::TRACKING_VERBOSITY);
  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();

  if (TRACKING::FastKalmanFilter == nullptr)
  {
    cout << __PRETTY_FUNCTION__ << " : missing the expected initialization for TRACKING::FastKalmanFilter." << endl;
    exit(1);
  }

  InitFastKalmanFilter(TRACKING::FastKalmanFilter);
  TRACKING::FastKalmanFilter->Verbosity(verbosity);
  TRACKING::FastKalmanFilter->set_trackmap_out_name(TRACKING::TrackNodeName);

  //-------------------------
  // FEMC
  //-------------------------
  // Saved track states (projections)
  if (Enable::FEMC && G4TRACKING::PROJECTION_FEMC)
  {
    TRACKING::FastKalmanFilter->add_state_name("FEMC");
    TRACKING::ProjectionNames.insert("FEMC");
  }

  //-------------------------
  // FHCAL
  //-------------------------
  if (Enable::FHCAL && G4TRACKING::PROJECTION_FHCAL)
  {
    TRACKING::FastKalmanFilter->add_state_name("FHCAL");
    TRACKING::ProjectionNames.insert("FHCAL");
  }
  //-------------------------
  // LFHCAL
  //-------------------------
  if (Enable::LFHCAL && G4TRACKING::PROJECTION_LFHCAL)
  {
    TRACKING::FastKalmanFilter->add_state_name("LFHCAL");
    TRACKING::ProjectionNames.insert("LFHCAL");
  }
  //-------------------------
  // CEMC
  //-------------------------
  if (Enable::CEMC && G4TRACKING::PROJECTION_CEMC)
  {
    TRACKING::FastKalmanFilter->add_state_name("CEMC");
    TRACKING::ProjectionNames.insert("CEMC");
  }
  //-------------------------
  // BECAL
  //-------------------------
  if (Enable::BECAL && G4TRACKING::PROJECTION_BECAL)
  {
    TRACKING::FastKalmanFilter->add_state_name("BECAL");
    TRACKING::ProjectionNames.insert("BECAL");
  }
  //-------------------------
  // HCALIN
  //-------------------------
  if (Enable::HCALIN && G4TRACKING::PROJECTION_HCALIN)
  {
    TRACKING::FastKalmanFilter->add_state_name("HCALIN");
    TRACKING::ProjectionNames.insert("HCALIN");
  }
  //-------------------------
  // HCALOUT
  //-------------------------
  if (Enable::HCALOUT && G4TRACKING::PROJECTION_HCALOUT)
  {
    TRACKING::FastKalmanFilter->add_state_name("HCALOUT");
    TRACKING::ProjectionNames.insert("HCALOUT");
  }
  //-------------------------
  // EEMC
  //-------------------------
  if ((Enable::EEMC or Enable::EEMCH)&& G4TRACKING::PROJECTION_EEMC)
  {
    TRACKING::FastKalmanFilter->add_state_name("EEMC");
    TRACKING::ProjectionNames.insert("EEMC");
  }
  //-------------------------
  // EHCAL
  //-------------------------
  if (Enable::EHCAL && G4TRACKING::PROJECTION_EHCAL)
  {
    TRACKING::FastKalmanFilter->add_state_name("EHCAL");
    TRACKING::ProjectionNames.insert("EHCAL");
  }

  se->registerSubsystem(TRACKING::FastKalmanFilter);

  // next, tracks with partial usage of the tracker stack
  if (TRACKING::FastKalmanFilterInnerTrack == nullptr)
  {
    cout << __PRETTY_FUNCTION__ << " : missing the expected initialization for TRACKING::FastKalmanFilterInnerTrack." << endl;
    exit(1);
  }
  InitFastKalmanFilter(TRACKING::FastKalmanFilterInnerTrack);
  TRACKING::FastKalmanFilterInnerTrack->Verbosity(verbosity);
  TRACKING::FastKalmanFilterInnerTrack->set_trackmap_out_name("InnerTrackMap");
  TRACKING::FastKalmanFilterInnerTrack->enable_vertexing(false);
  se->registerSubsystem(TRACKING::FastKalmanFilterInnerTrack);

  if (TRACKING::FastKalmanFilterSiliconTrack == nullptr)
  {
    cout << __PRETTY_FUNCTION__ << " : missing the expected initialization for TRACKING::FastKalmanFilterSiliconTrack." << endl;
    exit(1);
  }
  InitFastKalmanFilter(TRACKING::FastKalmanFilterSiliconTrack);
  TRACKING::FastKalmanFilterSiliconTrack->Verbosity(verbosity);
  TRACKING::FastKalmanFilterSiliconTrack->set_trackmap_out_name("SiliconTrackMap");
  TRACKING::FastKalmanFilterSiliconTrack->enable_vertexing(false);
  se->registerSubsystem(TRACKING::FastKalmanFilterSiliconTrack);
  return;
}

//-----------------------------------------------------------------------------//

void Tracking_Eval(const std::string &outputfile)
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::TRACKING_VERBOSITY);
  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();

  //----------------
  // Fast Tracking evaluation
  //----------------

  PHG4TrackFastSimEval *fast_sim_eval = new PHG4TrackFastSimEval("FastTrackingEval");
  fast_sim_eval->set_trackmapname(TRACKING::TrackNodeName);
  fast_sim_eval->set_filename(outputfile);
  fast_sim_eval->Verbosity(verbosity);
  //-------------------------
  // FEMC
  //-------------------------

  cout << "Tracking_Eval(): configuration of track projections in PHG4TrackFastSimEval" << endl;
  cout << "/*std::set<std::string>*/ TRACKING::ProjectionNames = {";
  bool first = true;
  for (const string &proj : TRACKING::ProjectionNames)
  {
    if (first)
      first = false;
    else
      cout << ", ";
    cout << "\"" << proj << "\"";

    fast_sim_eval->AddProjection(proj);
  }
  cout << "};" << endl;  // override the TRACKING::ProjectionNames in eval macros

  se->registerSubsystem(fast_sim_eval);

  // now partial track fits
  fast_sim_eval = new PHG4TrackFastSimEval("FastTrackingEval_InnerTrackMap");
  fast_sim_eval->set_trackmapname("InnerTrackMap");
  fast_sim_eval->set_filename(outputfile + ".InnerTrackMap.root");
  fast_sim_eval->Verbosity(verbosity);
  se->registerSubsystem(fast_sim_eval);

  fast_sim_eval = new PHG4TrackFastSimEval("FastTrackingEval_SiliconTrackMap");
  fast_sim_eval->set_trackmapname("SiliconTrackMap");
  fast_sim_eval->set_filename(outputfile + ".SiliconTrackMap.root");
  fast_sim_eval->Verbosity(verbosity);
  se->registerSubsystem(fast_sim_eval);
}
#endif
