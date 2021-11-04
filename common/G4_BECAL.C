#ifndef MACRO_G4BECAL_C
#define MACRO_G4BECAL_C

#include <GlobalVariables.C>

#include <g4calo/RawTowerDigitizer.h>

#include <g4eiccalos/PHG4BarrelEcalSubsystem.h>
#include <g4eiccalos/RawTowerBuilderByHitIndexBECAL.h>

#include <g4eval/CaloEvaluator.h>

#include <g4main/PHG4Reco.h>

#include <eiccaloreco/RawClusterBuilderkV3.h>
#include <eiccaloreco/RawClusterBuilderHelper.h>

#include <caloreco/RawClusterBuilderFwd.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawTowerCalibration.h>

#include <fun4all/Fun4AllServer.h>

R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libg4calo.so)
R__LOAD_LIBRARY(libg4eiccalos.so)
R__LOAD_LIBRARY(libg4eval.so)

namespace Enable
{
  bool BECAL = false;
  bool BECAL_ABSORBER = false;
  bool BECAL_CELL = false;
  bool BECAL_TOWER = false;
  bool BECAL_CLUSTER = false;
  bool BECAL_EVAL = false;
  bool BECAL_OVERLAPCHECK = false;
  int BECAL_VERBOSITY = 0;
}  // namespace Enable



namespace G4BECAL
{

  //double minz =  -273.6*cm;
  //double maxz =  142.4*cm;
  double minz = -453;
  double maxz = 371;
  double topradius =  138;
  double radius =  80;

  // this is default set to -1.5<eta<1.24 for 2018 Letter of Intent
  // if the user changes these, the z position of the

   // Digitization (default photon digi):
  RawTowerDigitizer::enu_digi_algorithm TowerDigi = RawTowerDigitizer::kSiPM_photon_digitization;
  // directly pass the energy of sim tower to digitized tower
  // kNo_digitization
  // simple digitization with photon statistics, single amplitude ADC conversion and pedestal
  // kSimple_photon_digitization
  // digitization with photon statistics on SiPM with an effective pixel N, ADC conversion and pedestal
  // kSiPM_photon_digitization

}  // namespace G4BECAL

void BECALInit()
{

  // update black hole settings
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4BECAL::topradius);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, G4BECAL::maxz);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, G4BECAL::minz);

}

double BECALSetup(PHG4Reco *g4Reco)
{

  bool AbsorberActive = Enable::ABSORBER || Enable::BECAL_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::BECAL_OVERLAPCHECK;

  // Update IR of BCAL to R~80cm
  // From Nathaly Santiesteban:
  // https://raw.githubusercontent.com/eic/fun4all_eiccalibrations/main/BarrelEcal/mapping/towerMap_BEMC_v002.txt
  // It uses IR of 80.3.
  ostringstream mapping_becal;
  mapping_becal << getenv("CALIBRATIONROOT") << "/BarrelEcal/mapping/towerMap_BEMC_v002.txt";

  PHG4BarrelEcalSubsystem *becal = new PHG4BarrelEcalSubsystem("BECAL");
  becal->set_string_param("mapping_file", mapping_becal.str());
  becal->OverlapCheck(OverlapCheck);
  becal->SetActive();
  becal->SuperDetector("BECAL");
  if (AbsorberActive) becal->SetAbsorberActive();

  g4Reco->registerSubsystem(becal);

  return  G4BECAL::topradius;

}

void BECAL_Cells(int verbosity = 0)
{
  return;
}

void BECAL_Towers()
{
  
  int verbosity = std::max(Enable::VERBOSITY, Enable::BECAL_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();

  ostringstream mapping_BECAL;
  mapping_BECAL << getenv("CALIBRATIONROOT") << "/BarrelEcal/mapping/towerMap_BEMC_v001.txt";
  
  const double photoelectron_per_GeV = 5000;

  RawTowerBuilderByHitIndexBECAL *tower_BECAL = new RawTowerBuilderByHitIndexBECAL("TowerBuilder_BECAL");
  tower_BECAL->Detector("BECAL");
  tower_BECAL->set_sim_tower_node_prefix("SIM");
  tower_BECAL->EminCut(1e-7);
  tower_BECAL->GeometryTableFile(mapping_BECAL.str());
  tower_BECAL->Verbosity(verbosity);
  se->registerSubsystem(tower_BECAL);

  RawTowerDigitizer *TowerDigitizer_BECAL = new RawTowerDigitizer("BECALRawTowerDigitizer");
  TowerDigitizer_BECAL->Detector("BECAL");
  TowerDigitizer_BECAL->Verbosity(verbosity);
//   TowerDigitizer_BECAL->Verbosity(verbosity);
  TowerDigitizer_BECAL->set_digi_algorithm(G4BECAL::TowerDigi);
  TowerDigitizer_BECAL->set_raw_tower_node_prefix("RAW");
  TowerDigitizer_BECAL->set_pedstal_central_ADC(0);
  TowerDigitizer_BECAL->set_pedstal_width_ADC(0);  // eRD1 test beam setting
  TowerDigitizer_BECAL->set_photonelec_ADC(1);     // not simulating ADC discretization error
  TowerDigitizer_BECAL->set_photonelec_yield_visible_GeV(photoelectron_per_GeV);
  TowerDigitizer_BECAL->set_zero_suppression_ADC(0);  // eRD1 test beam setting
  se->registerSubsystem(TowerDigitizer_BECAL); 

  RawTowerCalibration *TowerCalibration_BECAL = new RawTowerCalibration("BECALRawTowerCalibration");
  TowerCalibration_BECAL->Detector("BECAL");
  TowerCalibration_BECAL->Verbosity(verbosity);
  TowerCalibration_BECAL->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
  TowerCalibration_BECAL->set_calib_const_GeV_ADC(1. / photoelectron_per_GeV);
  TowerCalibration_BECAL->set_pedstal_ADC(0);
  se->registerSubsystem(TowerCalibration_BECAL);

}

void BECAL_Clusters()
{
  Fun4AllServer *se = Fun4AllServer::instance();

  RawClusterBuilderHelper *ClusterBuilder = new RawClusterBuilderkV3("BECALRawClusterBuilderkV3");
  ClusterBuilder->Detector("BECAL");
  ClusterBuilder->set_seed_e(0.5);
  ClusterBuilder->set_agg_e(0.1);
  se->registerSubsystem(ClusterBuilder);

  return;
}

void BECAL_Eval(const std::string &outputfile)
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::BECAL_VERBOSITY);
  Fun4AllServer *se = Fun4AllServer::instance();

  CaloEvaluator *eval = new CaloEvaluator("BECALEVALUATOR", "BECAL", outputfile.c_str());
  eval->set_do_cluster_eval(false);
  eval->Verbosity(1);
  se->registerSubsystem(eval);

  return;
}
#endif
