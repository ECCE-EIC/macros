#ifndef MACRO_G4BMMG_C
#define MACRO_G4BMMG_C

#include <GlobalVariables.C>

#include <fun4all/Fun4AllServer.h>
#include <g4barrelmmg/CreateCZHitContainer.h>
#include <g4barrelmmg/PHG4CylinderStripSubsystem.h>
#include <g4main/PHG4Reco.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>

R__LOAD_LIBRARY(libg4barrelmmg.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4testbench.so)
R__LOAD_LIBRARY(libg4detectors.so)

namespace Enable
{
  bool BMMG = false;
  bool use_2Dreadout = true;
  bool BMMG_OVERLAPCHECK = true;
  int BMMG_VERBOSITY = 0;
}  // namespace Enable

namespace BMMG
{
  const int n_layer = 3;

  const double rad[BMMG::n_layer] = {45., 47.4, 67.4};  // approximate radial location
  const double len[BMMG::n_layer] = {140, 150, 280.0};
}  // namespace BMMG

void BMMGInit(int verbosity = 1)
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, BMMG::rad[BMMG::n_layer - 1] / 10. + 0.7);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, BMMG::len[BMMG::n_layer - 1] / 2.0 + 1.0);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -BMMG::len[BMMG::n_layer - 1] / 2.0 - 10.);
}

void BMMGSetup(PHG4Reco* g4Reco)
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::BMMG_VERBOSITY);
  Fun4AllServer* se = Fun4AllServer::instance();
  se->Verbosity(verbosity);

  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::BMMG_OVERLAPCHECK;

  gSystem->Load("libfun4all");
  gSystem->Load("libg4detectors.so");
  gSystem->Load("libg4testbench.so");
  gSystem->Load("libg4trackfastsim.so");

  PHG4CylinderStripSubsystem* barrel_mmg;
  double gap_betweenCZ = 1.5;
  double Gap_betweenlayer = 1.5;
  double thickness = 0.36499;
  int nCZlayer = 2;
  if (Enable::use_2Dreadout)
  {
    gap_betweenCZ = 0;
    nCZlayer = 1;
  }

  const double prapidity = 1;

  for (int ilayer = 0; ilayer < BMMG::n_layer; ilayer++)
  {
    barrel_mmg = new PHG4CylinderStripSubsystem(Form("BMT_%d", ilayer), ilayer);
    barrel_mmg->set_double_param("radius", BMMG::rad[ilayer]);
    barrel_mmg->set_string_param("gas", "myMMGas");
    //barrel_mmg->set_double_param("steplimits", 300e-4);
    barrel_mmg->set_double_param("phi0", 15 * ilayer);
    barrel_mmg->set_double_param("gap", gap_betweenCZ);
    barrel_mmg->SetActive();
    barrel_mmg->SuperDetector("BMT");
    barrel_mmg->set_int_param("lengthviarapidity", 0);
    barrel_mmg->set_double_param("gas1thickness", 0.15);
    barrel_mmg->set_double_param("length", BMMG::len[ilayer]);
    barrel_mmg->set_double_param("deadzone", 0.2);
    barrel_mmg->set_int_param("nhit", 1);
    barrel_mmg->OverlapCheck(true);
    barrel_mmg->set_int_param("use_2Dreadout", Enable::use_2Dreadout);
    g4Reco->registerSubsystem(barrel_mmg);
  }  // ilayer loop

  if (TRACKING::FastKalmanFilter)
  {
    if (Enable::use_2Dreadout)
    {
      TRACKING::FastKalmanFilter->add_phg4hits(
          "G4HIT_BMT",                 //      const std::string& phg4hitsNames,
          PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
          2.5 / 2 / sqrt(12),          //1./sqrt(12),                      //       radial-resolution , only used for Vertical Plane Detector Type
          75e-4,                       //150e-4,                       //        azimuthal-resolution [cm]
          75e-4,                       //150e-4,                           //      z-resolution [cm]
          1,                           //      efficiency,
          0);                          //      noise hits
    }
    else
    {
      TRACKING::FastKalmanFilter->add_phg4hits(
          "G4HIT_CZBMT",               //      const std::string& phg4hitsNames,
          PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
          2.5 / 2 / sqrt(12),          //1/sqrt(12),                      //       radial-resolution [cm], only used for Vertical Plane Detector Type
          75e-4,                       //150e-4,                       //        azimuthal-resolution [cm]
          75e-4,                       //150e-4,                           //      z-resolution [cm]
          1,                           //      efficiency,
          0);                          //      noise hits
    }
  }

  return;
}

void BMMGT_Reco()
{
  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4detectors.so");

  int verbosity = std::max(Enable::VERBOSITY, Enable::BMMG_VERBOSITY);
  Fun4AllServer* se = Fun4AllServer::instance();

  return;
}
#endif
