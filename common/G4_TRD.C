#ifndef MACRO_G4TRD_C
#define MACRO_G4TRD_C

#include <GlobalVariables.C>

#include <fun4all/Fun4AllServer.h>
#include <g4main/PHG4Reco.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4trd/PHG4TRDSubsystem.h>

R__LOAD_LIBRARY(libg4trd.so)
R__LOAD_LIBRARY(libg4detectors.so)

namespace Enable
{
  bool TRD = false;
  //bool TRD_ABSORBER = false;
  bool TRD_GAS = false;
  bool TRD_OVERLAPCHECK = true;
  int TRD_VERBOSITY = 0;
}  // namespace Enable

namespace TRD
{
  double R_max = 180.;
  double R_min = 15.;
  double z_mid = 200. + 20;
  double half_length = 6.5;  // Mother volume, radiator = 10 cm
  double z_min = z_mid - half_length;
  double z_max = z_mid + half_length + 10.;
}  // namespace TRD
void TRDInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, TRD::R_max);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, TRD::z_max);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, TRD::z_min);
}

void TRDSetup(PHG4Reco* g4Reco)
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::TRD_VERBOSITY);
  Fun4AllServer* se = Fun4AllServer::instance();
  se->Verbosity(verbosity);

  //bool AbsorberActive = Enable::ABSORBER || Enable::TRD_ABSORBER;
  bool GasActive = Enable::ABSORBER || Enable::TRD_GAS;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::TRD_OVERLAPCHECK;
  PHG4TRDSubsystem* trd_hcap = new PHG4TRDSubsystem("TRD_hcap", 1);

  trd_hcap->set_double_param("ThicknessZ", 2. * TRD::half_length);
  trd_hcap->set_double_param("RIn", TRD::R_min);
  trd_hcap->set_double_param("ROut", TRD::R_max);
  trd_hcap->set_double_param("PosZ", TRD::z_mid);
  trd_hcap->set_double_param("det_RIn", TRD::R_min);
  trd_hcap->set_double_param("det_ROut", TRD::R_max);
  trd_hcap->SetActive(1);
  trd_hcap->SuperDetector("TRD");
  //if(AbsorberActive)
  if (GasActive)
  {
    trd_hcap->SetAbsorberActive(1);
  }
  trd_hcap->OverlapCheck(OverlapCheck);
  g4Reco->registerSubsystem(trd_hcap);

  if (verbosity > 0)
  {
    cout << " Min Z  :"
         << " " << TRD::z_min << " Max Z  :" << TRD::z_max << endl;
    cout << "=======End setting parameters to geometry : ============" << endl;
  }
  if (TRACKING::FastKalmanFilter)
  {
    TRACKING::FastKalmanFilter->add_zplane_state("TRD", 284.);
    TRACKING::ProjectionNames.insert("TRD");

    //Use hits from  absorber (Active Gas volume)
    TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + string(Form("ACTIVEGAS_TRD")),  //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Vertical_Plane,                  //      const DETECTOR_TYPE phg4dettype,
                                             1. / sqrt(12.),                                    //      const float radres,
                                             70e-4,                                             //      const float phires,
                                             100e-4,                                            //      const float lonres,
                                             1,                                                 //      const float eff,
                                             0);
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
