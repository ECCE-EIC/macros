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
  //bool TRD_ABSORBER = false;
  bool TRD_GAS = false;
  bool TRD_OVERLAPCHECK = true;
  int TRD_VERBOSITY = 1;
}

namespace TRD
{
  double R_max =  120.; 
  double R_min = 5.; 
  double z_mid = 200.;
  double half_length = 6.5; // Mother volume, radiator = 10 cm
  double z_min = z_mid - half_length;
  double z_max = z_mid + half_length + 10.;  
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
  se->Verbosity(INT_MAX-10);
  
  //bool AbsorberActive = Enable::ABSORBER || Enable::TRD_ABSORBER;
  bool GasActive = Enable::ABSORBER || Enable::TRD_GAS;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::TRD_OVERLAPCHECK;
  PHG4TRDSubsystem* trd_hcap = new PHG4TRDSubsystem("TRD_hcap", 1);
  //// Mother volume dimensions, 
  double ThicknessZ = 13. ; // Mother  thickness along Z in cm, not to go below this
  double RIn  = 5. ; // Inner radius in cm
  double ROut = 120; // Outer radius in cm
  double PosZ  = 200.;  // Center of TRD mother volume  @ Z in cm 
  //Daughter volume (radiator and absorber) inner and outer radii, rest are in the detector construction class.
  double det_RIn = 5;
  double det_ROut = 120;


  trd_hcap->set_double_param("ThicknessZ", ThicknessZ);
  trd_hcap->set_double_param("RIn", RIn);
  trd_hcap->set_double_param("ROut", ROut);
  trd_hcap->set_double_param("PosZ", PosZ);
  trd_hcap->set_double_param("det_RIn", det_RIn);
  trd_hcap->set_double_param("det_ROut", det_ROut);
  trd_hcap->SetActive(1);
  trd_hcap->SuperDetector("TRD");
  //if(AbsorberActive)
  if(GasActive)
    {
      trd_hcap->SetAbsorberActive(1);
    }
  trd_hcap->OverlapCheck(OverlapCheck);
  //trd_hcap->OverlapCheck(1);
  g4Reco->registerSubsystem(trd_hcap);
  
  cout << " Min Z  :" << " " << TRD::z_min  << " Max Z  :" << TRD::z_max << endl; 
  cout << "=======End setting parameters to geometry : ============" << endl; 
  if (TRACKING::FastKalmanFilter)
    {
      TRACKING::FastKalmanFilter-> add_zplane_state("TRD", 208. );
      TRACKING::ProjectionNames.insert("TRD");
     
      //Use hits from  absorber (Active Gas volume)
      TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + string(Form("ACTIVEGAS_TRD")),  //      const std::string& phg4hitsNames,
					       PHG4TrackFastSim::Vertical_Plane,                       //      const DETECTOR_TYPE phg4dettype,
					       1. / sqrt(12.),                                     //      const float radres,
					       70e-4,                                              //      const float phires,
					       100e-4,                                              //      const float lonres,
					       1,                                                  //      const float eff,
                                               0);      
      
    }
 
  cout << " =====End Fast Kalman Filter ======  " << endl; 
  return;
 
}

void TRD_Reco()
{
  
  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4detectors.so");
  
  
  int verbosity = std::max(Enable::VERBOSITY, Enable::TRD_VERBOSITY);
  Fun4AllServer* se = Fun4AllServer::instance();
  //se->Verbosity(INT_MAX-10);
  
  return;

}
#endif




