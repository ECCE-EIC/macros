#ifndef MACRO_G4TRD_C
#define MACRO_G4TRD_C

#include <GlobalVariables.C>

#include <fun4all/Fun4AllServer.h>
#include <g4trd/PHG4TRDSubsystem.h>
#include <g4main/PHG4Reco.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>

R__LOAD_LIBRARY(libg4trd.so)
R__LOAD_LIBRARY(libg4detectors.so)

namespace Enable
{
  bool TRD = false;
  bool TRD_OVERLAPCHECK = false;
  int TRD_VERBOSITY = 0;
}

namespace TRD
{
  double R_max =  20.; 
  double R_min = 5.; 
  double z_mid = 200.;
  double half_length = 12.5;
  double z_min = z_mid - half_length;
  double z_max = z_mid + half_length;  
}
void TRDInit()
{

  BlackHoleGeometry::max_radius =  std::max(BlackHoleGeometry::max_radius, TRD::R_max);
  BlackHoleGeometry::max_z =  std::max(BlackHoleGeometry::max_z, TRD::z_max);
  BlackHoleGeometry::min_z =  std::min(BlackHoleGeometry::min_z, TRD::z_min);

}

void TRDSetup(PHG4Reco *g4Reco)
{

  Fun4AllServer* se = Fun4AllServer::instance();
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::TRD_OVERLAPCHECK;
  PHG4TRDSubsystem* trd_hcap = new PHG4TRDSubsystem("TRD_hcap", 1);
  //// Mother volume dimensions, 
  double ThicknessZ = 25. ; //  thickness along Z in cm
  double RIn  = 5. ; // Inner radius in cm
  double ROut = 140; // Outer radius in cm
  double PosZ  = 200.;  // Center of TRD @ Z in cm 
  //Daughter volume (radiator and absorber) inner and outer radii, rest are in the detector construction class.
  double det_RIn = 5;
  double det_ROut = 140;


  trd_hcap->set_double_param("ThicknessZ", ThicknessZ);
  trd_hcap->set_double_param("RIn", RIn);
  trd_hcap->set_double_param("ROut", ROut);
  trd_hcap->set_double_param("PosZ", PosZ);
  trd_hcap->set_double_param("det_RIn", det_RIn);
  trd_hcap->set_double_param("det_ROut", det_ROut);
  trd_hcap->SetActive(1);
  trd_hcap->SuperDetector("TRD");
  trd_hcap->OverlapCheck(OverlapCheck);
  g4Reco->registerSubsystem(trd_hcap);

 
  if (TRACKING::FastKalmanFilter)
    {
      TRACKING::FastKalmanFilter-> add_zplane_state("TRD", 212.0 );
      TRACKING::ProjectionNames.insert("TRD");

    }


  return;

}

void TRD_Reco()
{

  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4detectors.so");


  int verbosity = std::max(Enable::VERBOSITY, Enable::TRD_VERBOSITY);
  Fun4AllServer* se = Fun4AllServer::instance();

  return;

}
#endif




