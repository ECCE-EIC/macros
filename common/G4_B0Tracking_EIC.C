#ifndef MACRO_G4B0TRACKINGEIC_C
#define MACRO_G4B0TRACKINGEIC_C

#include <GlobalVariables.C>

#include <trackreco/PHRaveVertexing.h>

#include <eicg4b0/B0TrackFastSimEval.h>
#include <eicg4b0/B0TrackFastSim.h>

#include <fun4all/Fun4AllServer.h>

#include <vector>

R__LOAD_LIBRARY(libtrack_reco.so)
R__LOAD_LIBRARY(libEICG4B0.so)//

namespace Enable
{
  bool B0TRACKING = false;
  bool B0TRACKING_EVAL = false;
  int B0TRACKING_VERBOSITY = 0;
}  // namespace Enable

namespace G4B0TRACKING
{
  bool DISPLACED_VERTEX = false;
  bool PROJECTION_B0 = false;
  const double PositionResolution(30e-4);
}  // namespace G4TRACKING

//-----------------------------------------------------------------------------//
void B0TrackingInit()
{
  B0TRACKING::FastKalmanFilter = new B0TrackFastSim("B0TrackFastSim");
  B0TRACKING::FastKalmanFilterB0Track = new B0TrackFastSim("FastKalmanFilterB0Track"); //b0Truth
}

void InitFastKalmanFilter(B0TrackFastSim *kalman_filter)
{
  //  kalman_filter->Smearing(false);
  if (G4B0TRACKING::DISPLACED_VERTEX)
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
    kalman_filter->set_vertex_xy_resolution(G4B0TRACKING::PositionResolution);
    kalman_filter->set_vertex_z_resolution(G4B0TRACKING::PositionResolution);
  }

  kalman_filter->set_sub_top_node_name("TRACKS");
}

//-----------------------------------------------------------------------------//
void B0Tracking_Reco()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::B0TRACKING_VERBOSITY);
  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();

  if (B0TRACKING::FastKalmanFilter == nullptr)
  {
    cout << __PRETTY_FUNCTION__ << " : missing the expected initialization for B0TRACKING::FastKalmanFilter." << endl;
    exit(1);
  }

  InitFastKalmanFilter(B0TRACKING::FastKalmanFilter);
  B0TRACKING::FastKalmanFilter->Verbosity(verbosity);
  B0TRACKING::FastKalmanFilter->set_trackmap_out_name(B0TRACKING::TrackNodeName);

  //-------------------------
  // B0
  //-------------------------
  if (Enable::B0TRACKING && G4B0TRACKING::PROJECTION_B0)
  {
//    TRACKING::FastKalmanFilter->add_state_name("b0Truth");
//    TRACKING::ProjectionNames.insert("b0Truth");
  }

  se->registerSubsystem(B0TRACKING::FastKalmanFilter);

  // next, tracks with partial usage of the tracker stack
  if (B0TRACKING::FastKalmanFilterB0Track == nullptr)
  {
    cout << __PRETTY_FUNCTION__ << " : missing the expected initialization for TRACKING::FastKalmanFilterB0Track." << endl;
    exit(1);
  }
  InitFastKalmanFilter(B0TRACKING::FastKalmanFilterB0Track);
  B0TRACKING::FastKalmanFilterB0Track->Verbosity(verbosity);
  B0TRACKING::FastKalmanFilterB0Track->set_trackmap_out_name("B0TrackMap");
  B0TRACKING::FastKalmanFilterB0Track->enable_vertexing(false);
  se->registerSubsystem(B0TRACKING::FastKalmanFilterB0Track);
  return;
}

//-----------------------------------------------------------------------------//

void B0Tracking_Eval(const std::string &outputfile)
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::B0TRACKING_VERBOSITY);
  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();

  //----------------
  // Fast Tracking evaluation
  //----------------

  B0TrackFastSimEval *b0_fast_sim_eval = new B0TrackFastSimEval("B0FastTrackingEval");
  b0_fast_sim_eval->set_trackmapname(B0TRACKING::TrackNodeName);
  b0_fast_sim_eval->set_filename(outputfile);
  b0_fast_sim_eval->Verbosity(verbosity);
  //-------------------------
  // FEMC
  //-------------------------

  cout << "Tracking_Eval(): configuration of track projections in B0TrackFastSimEval" << endl;
  cout << "std::set<std::string> B0TRACKING::ProjectionNames = {";
  bool first = true;
  for (const string &proj : B0TRACKING::B0ProjectionNames)
  {
    if (first)
      first = false;
    else
      cout << ", ";
    cout << "\"" << proj << "\"";

    b0_fast_sim_eval->AddProjection(proj);
  }
  cout << "};" << endl;  // override the TRACKING::ProjectionNames in eval macros

  se->registerSubsystem(b0_fast_sim_eval);
  // now partial track fits
  //B0TrackFastSimEval *b0_fast_sim_eval = new B0TrackFastSimEval("FastTrackingEval_B0TrackMap");
  b0_fast_sim_eval = new B0TrackFastSimEval("FastTrackingEval_B0TrackMap");
  b0_fast_sim_eval->set_trackmapname("B0TrackMap");
  b0_fast_sim_eval->set_filename(outputfile + ".B0TrackMap_debug.root");
  b0_fast_sim_eval->Verbosity(verbosity);
  se->registerSubsystem(b0_fast_sim_eval);
}
#endif
