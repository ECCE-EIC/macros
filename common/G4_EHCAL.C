#ifndef MACRO_G4EHCAL_C
#define MACRO_G4EHCAL_C

#include <GlobalVariables.C>

#include <g4calo/RawTowerBuilderByHitIndex.h>
#include <g4calo/RawTowerDigitizer.h>

#include <g4eiccalos/PHG4ForwardCalCellReco.h>
#include <g4eiccalos/PHG4ForwardHcalSubsystem.h>

#include <g4eval/CaloEvaluator.h>

#include <g4main/PHG4Reco.h>

#include <eiccaloreco/RawClusterBuilderkMA.h>
#include <eiccaloreco/RawClusterBuilderHelper.h>

#include <caloreco/RawClusterBuilderFwd.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawTowerCalibration.h>

#include <fun4all/Fun4AllServer.h>

R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libg4calo.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4eval.so)

namespace Enable
{
  bool EHCAL = false;
  bool EHCAL_ABSORBER = false;
  bool EHCAL_CELL = false;
  bool EHCAL_TOWER = false;
  bool EHCAL_CLUSTER = false;
  bool EHCAL_EVAL = false;
  bool EHCAL_OVERLAPCHECK = false;
  int EHCAL_VERBOSITY = 0;
}  // namespace Enable

namespace G4EHCAL
{
  // from ForwardHcal/mapping/towerMap_EHCAL_v005.txt
  double Gz0 = 400.;
  double Gdz = 100.;
  double outer_radius = 262.;
  enum enu_EHCAL_clusterizer
  {
    kEHCALGraphClusterizer,
    kEHCALTemplateClusterizer
  };
  //template clusterizer, as developed by Sasha Bazilevsky
  enu_EHCAL_clusterizer EHCAL_clusterizer = kEHCALTemplateClusterizer;
  // graph clusterizer
  //enu_FHcal_clusterizer FHcal_clusterizer = kFHcalGraphClusterizer;
  namespace SETTING
  {
    bool FullEtaAcc = false;
    bool HC2x = false;
    bool HC4x = false;
    bool towercalib1 = false;
    bool towercalibSiPM = false;
    bool towercalibHCALIN = false;
    bool towercalib3 = false;
  }  // namespace SETTING
}  // namespace G4EHCAL

void EHCALInit()
{
  // simple way to check if only 1 of the settings is true
  if ((G4EHCAL::SETTING::FullEtaAcc ? 1 : 0) + (G4EHCAL::SETTING::HC4x ? 1 : 0) + (G4EHCAL::SETTING::HC2x ? 1 : 0) > 1)
  {
    cout << "use only  G4EHCAL::SETTING::FullEtaAcc=true or G4EHCAL::SETTING::HC2x=true or G4EHCAL::SETTING::HC4x=true" << endl;
    gSystem->Exit(1);
  }
  if ((G4EHCAL::SETTING::towercalib1 ? 1 : 0) + (G4EHCAL::SETTING::towercalibSiPM ? 1 : 0) +
          (G4EHCAL::SETTING::towercalibHCALIN ? 1 : 0) + (G4EHCAL::SETTING::towercalib3 ? 1 : 0) >
      1)
  {
    cout << "use only G4EHCAL::SETTING::towercalib1 = true or G4EHCAL::SETTING::towercalibSiPM = true"
         << " or G4EHCAL::SETTING::towercalibHCALIN = true or G4EHCAL::SETTING::towercalib3 = true" << endl;
    gSystem->Exit(1);
  }

  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4EHCAL::outer_radius);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, G4EHCAL::Gz0 + G4EHCAL::Gdz / 2.);
}

void EHCALSetup(PHG4Reco *g4Reco)
{
  const bool AbsorberActive = Enable::ABSORBER || Enable::EHCAL_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::EHCAL_OVERLAPCHECK;
  Fun4AllServer *se = Fun4AllServer::instance();

  /** Use dedicated EHCAL module */
  PHG4ForwardHcalSubsystem *ehcal = new PHG4ForwardHcalSubsystem("EHCAL");

  ostringstream mapping_EHCAL;

  // Switch to desired calo setup
  // HCal Fe-Scint with doubled granularity
  if (G4EHCAL::SETTING::HC2x )
  {
    mapping_EHCAL << getenv("CALIBRATIONROOT") << "/BackwardHcal/mapping/towerMap_EHCAL_2x.txt";
  }
  // full HCal Fe-Scint with nominal acceptance doubled granularity
  else if (G4EHCAL::SETTING::HC2x && G4EHCAL::SETTING::FullEtaAcc)
  {
    mapping_EHCAL << getenv("CALIBRATIONROOT") << "/BackwardHcal/mapping/towerMap_EHCAL_2x_fullEtaCov.txt";
  }
  // HCal Fe-Scint with four times granularity
  else if (G4EHCAL::SETTING::HC4x )
  {
    mapping_EHCAL << getenv("CALIBRATIONROOT") << "/BackwardHcal/mapping/towerMap_EHCAL_4x.txt";
  }
  // full HCal Fe-Scint with nominal acceptance four times granularity
  else if (G4EHCAL::SETTING::HC4x && G4EHCAL::SETTING::FullEtaAcc)
  {
    mapping_EHCAL << getenv("CALIBRATIONROOT") << "/BackwardHcal/mapping/towerMap_EHCAL_4x_fullEtaCov.txt";
  }
  // full HCal Fe-Scint with nominal acceptance
  else if (G4EHCAL::SETTING::FullEtaAcc)
  {
    mapping_EHCAL << getenv("CALIBRATIONROOT") << "/BackwardHcal/mapping/towerMap_EHCAL_default_fullEtaCov.txt";
  }
  // full HCal Fe-Scint with enlarged beam pipe opening for Mar 2020 beam pipe
  else
  {
    mapping_EHCAL << getenv("CALIBRATIONROOT")
                  << "/BackwardHcal/mapping/towerMap_EHCAL_default.txt";
  }

  ehcal->SetTowerMappingFile(mapping_EHCAL.str());
  ehcal->OverlapCheck(OverlapCheck);
  ehcal->SetActive();
  ehcal->SuperDetector("EHCAL");
  if (AbsorberActive) ehcal->SetAbsorberActive();

  g4Reco->registerSubsystem(ehcal);
}

void EHCAL_Cells(int verbosity = 0)
{
  return;
}

void EHCAL_Towers()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::EHCAL_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();

  ostringstream mapping_EHCAL;

  // Switch to desired calo setup
  // HCal Fe-Scint with doubled granularity
  if (G4EHCAL::SETTING::HC2x )
  {
    mapping_EHCAL << getenv("CALIBRATIONROOT") << "/BackwardHcal/mapping/towerMap_EHCAL_2x.txt";
  }
  // full HCal Fe-Scint with nominal acceptance doubled granularity
  else if (G4EHCAL::SETTING::HC2x && G4EHCAL::SETTING::FullEtaAcc)
  {
    mapping_EHCAL << getenv("CALIBRATIONROOT") << "/BackwardHcal/mapping/towerMap_EHCAL_2x_fullEtaCov.txt";
  }
  // HCal Fe-Scint with four times granularity
  else if (G4EHCAL::SETTING::HC4x )
  {
    mapping_EHCAL << getenv("CALIBRATIONROOT") << "/BackwardHcal/mapping/towerMap_EHCAL_4x.txt";
  }
  // full HCal Fe-Scint with nominal acceptance four times granularity
  else if (G4EHCAL::SETTING::HC4x && G4EHCAL::SETTING::FullEtaAcc)
  {
    mapping_EHCAL << getenv("CALIBRATIONROOT") << "/BackwardHcal/mapping/towerMap_EHCAL_4x_fullEtaCov.txt";
  }
  // full HCal Fe-Scint with nominal acceptance
  else if (G4EHCAL::SETTING::FullEtaAcc)
  {
    mapping_EHCAL << getenv("CALIBRATIONROOT") << "/BackwardHcal/mapping/towerMap_EHCAL_default_fullEtaCov.txt";
  }
  // full HCal Fe-Scint with enlarged beam pipe opening for Mar 2020 beam pipe
  else
  {
    mapping_EHCAL << getenv("CALIBRATIONROOT") << "/BackwardHcal/mapping/towerMap_EHCAL_default.txt";
  }

  RawTowerBuilderByHitIndex *tower_EHCAL = new RawTowerBuilderByHitIndex("TowerBuilder_EHCAL");
  tower_EHCAL->Detector("EHCAL");
  tower_EHCAL->set_sim_tower_node_prefix("SIM");
  tower_EHCAL->GeometryTableFile(mapping_EHCAL.str());

  se->registerSubsystem(tower_EHCAL);

  // enable usage of different tower calibrations for systematic studies
  if (G4EHCAL::SETTING::towercalib1)
  {
    cout << "1: using towercalib1 for EHCAL towers" << endl;
    const double EHCAL_photoelectron_per_GeV = 500;
    RawTowerDigitizer *TowerDigitizer_EHCAL = new RawTowerDigitizer("EHCALRawTowerDigitizer");

    TowerDigitizer_EHCAL->Detector("EHCAL");
    TowerDigitizer_EHCAL->Verbosity(verbosity);
    TowerDigitizer_EHCAL->set_raw_tower_node_prefix("RAW");
    TowerDigitizer_EHCAL->set_digi_algorithm(RawTowerDigitizer::kSiPM_photon_digitization);
    TowerDigitizer_EHCAL->set_pedstal_central_ADC(0);
    TowerDigitizer_EHCAL->set_pedstal_width_ADC(8);  // eRD1 test beam setting
    TowerDigitizer_EHCAL->set_photonelec_ADC(1);     //not simulating ADC discretization error
    TowerDigitizer_EHCAL->set_photonelec_yield_visible_GeV(EHCAL_photoelectron_per_GeV);
    TowerDigitizer_EHCAL->set_zero_suppression_ADC(16);  // eRD1 test beam setting

    se->registerSubsystem(TowerDigitizer_EHCAL);

    RawTowerCalibration *TowerCalibration_EHCAL = new RawTowerCalibration("EHCALRawTowerCalibration");
    TowerCalibration_EHCAL->Detector("EHCAL");
    TowerCalibration_EHCAL->Verbosity(verbosity);
    TowerCalibration_EHCAL->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
    TowerCalibration_EHCAL->set_calib_const_GeV_ADC(1. / EHCAL_photoelectron_per_GeV);
    TowerCalibration_EHCAL->set_pedstal_ADC(0);

    se->registerSubsystem(TowerCalibration_EHCAL);
  }
  else if (G4EHCAL::SETTING::towercalibSiPM)
  {
    //from https://sphenix-collaboration.github.io/doxygen/d4/d58/Fun4All__G4__Prototype4_8C_source.html
    const double sampling_fraction = 0.019441;     //  +/-  0.019441 from 0 Degree indenting 12 GeV electron showers
    const double photoelectron_per_GeV = 500;      //500 photon per total GeV deposition
    const double ADC_per_photoelectron_HG = 3.8;   // From Sean Stoll, Mar 29
    const double ADC_per_photoelectron_LG = 0.24;  // From Sean Stoll, Mar 29

    cout << "2: using towercalibSiPM for EHCAL towers" << endl;
    RawTowerDigitizer *TowerDigitizer = new RawTowerDigitizer("EHCALRawTowerDigitizer");
    TowerDigitizer->Detector("EHCAL");
    TowerDigitizer->set_raw_tower_node_prefix("RAW");
    TowerDigitizer->set_digi_algorithm(RawTowerDigitizer::kSiPM_photon_digitization);
    TowerDigitizer->set_pedstal_central_ADC(0);
    TowerDigitizer->set_pedstal_width_ADC(1);
    TowerDigitizer->set_photonelec_ADC(1. / ADC_per_photoelectron_LG);
    TowerDigitizer->set_photonelec_yield_visible_GeV(photoelectron_per_GeV / sampling_fraction);
    TowerDigitizer->set_zero_suppression_ADC(-1000);  // no-zero suppression
    se->registerSubsystem(TowerDigitizer);

    RawTowerCalibration *TowerCalibration = new RawTowerCalibration("EHCALRawTowerCalibration");
    TowerCalibration->Detector("EHCAL");
    TowerCalibration->set_raw_tower_node_prefix("RAW");
    TowerCalibration->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
    TowerCalibration->set_calib_const_GeV_ADC(1. / ADC_per_photoelectron_LG / photoelectron_per_GeV);
    TowerCalibration->set_pedstal_ADC(0);
    TowerCalibration->set_zero_suppression_GeV(-1);  // no-zero suppression
    se->registerSubsystem(TowerCalibration);
  }
  else if (G4EHCAL::SETTING::towercalibHCALIN)
  {
    const double visible_sample_fraction_HCALIN = 7.19505e-02;  // 1.34152e-02
    RawTowerDigitizer *TowerDigitizer = new RawTowerDigitizer("EHCALRawTowerDigitizer");
    TowerDigitizer->Detector("EHCAL");
    TowerDigitizer->set_raw_tower_node_prefix("RAW");
    TowerDigitizer->set_digi_algorithm(RawTowerDigitizer::kSimple_photon_digitalization);
    TowerDigitizer->set_pedstal_central_ADC(0);
    TowerDigitizer->set_pedstal_width_ADC(1);
    TowerDigitizer->set_photonelec_ADC(32. / 5.);
    TowerDigitizer->set_photonelec_yield_visible_GeV(32. / 5 / (0.4e-3));
    TowerDigitizer->set_zero_suppression_ADC(-1000);  // no-zero suppression
    se->registerSubsystem(TowerDigitizer);

    RawTowerCalibration *TowerCalibration = new RawTowerCalibration("EHCALRawTowerCalibration");
    TowerCalibration->Detector("EHCAL");
    TowerCalibration->set_raw_tower_node_prefix("RAW");
    TowerCalibration->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
    TowerCalibration->set_calib_const_GeV_ADC(0.4e-3 / visible_sample_fraction_HCALIN);
    TowerCalibration->set_pedstal_ADC(0);
    TowerCalibration->set_zero_suppression_GeV(-1);  // no-zero suppression
    se->registerSubsystem(TowerCalibration);
  }
  else if (G4EHCAL::SETTING::towercalib3)
  {
    cout << "3: using towercalib3 for EHCAL towers" << endl;
    RawTowerDigitizer *TowerDigitizer = new RawTowerDigitizer("EHCALRawTowerDigitizer");
    TowerDigitizer->Detector("EHCAL");
    TowerDigitizer->set_pedstal_central_ADC(0);
    TowerDigitizer->set_pedstal_width_ADC(8);  // eRD1 test beam setting
    TowerDigitizer->Verbosity(verbosity);
    TowerDigitizer->set_digi_algorithm(RawTowerDigitizer::kNo_digitization);
    se->registerSubsystem(TowerDigitizer);

    RawTowerCalibration *TowerCalibration = new RawTowerCalibration("EHCALRawTowerCalibration");
    TowerCalibration->Detector("EHCAL");
    TowerCalibration->Verbosity(verbosity);
    TowerCalibration->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
    TowerCalibration->set_calib_const_GeV_ADC(1. / 0.03898);  // calibrated with muons
    TowerCalibration->set_pedstal_ADC(0);
    se->registerSubsystem(TowerCalibration);
  }
  else
  {
    cout << "def: using default for EHCAL towers" << endl;
    RawTowerDigitizer *TowerDigitizer = new RawTowerDigitizer("EHCALRawTowerDigitizer");
    TowerDigitizer->Detector("EHCAL");
    TowerDigitizer->Verbosity(verbosity);
    TowerDigitizer->set_digi_algorithm(RawTowerDigitizer::kNo_digitization);
    se->registerSubsystem(TowerDigitizer);

    RawTowerCalibration *TowerCalibration = new RawTowerCalibration("EHCALRawTowerCalibration");
    TowerCalibration->Detector("EHCAL");
    TowerCalibration->Verbosity(verbosity);
    TowerCalibration->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
    TowerCalibration->set_calib_const_GeV_ADC(1. / 0.03898);  // calibrated with muons
    TowerCalibration->set_pedstal_ADC(0);
    se->registerSubsystem(TowerCalibration);
  }
}

void EHCAL_Clusters()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::EHCAL_VERBOSITY);
  Fun4AllServer *se = Fun4AllServer::instance();

  if (G4EHCAL::EHCAL_clusterizer == G4EHCAL::kEHCALTemplateClusterizer)
  {
    RawClusterBuilderHelper *ehcal_clusterbuilder = new RawClusterBuilderkMA("EhcalRawClusterBuilderkMA");
    ehcal_clusterbuilder->Detector("EHCAL");
    ehcal_clusterbuilder->set_seed_e(0.1);
    ehcal_clusterbuilder->set_agg_e(0.005);
    se->registerSubsystem(ehcal_clusterbuilder);
    /*
    RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EHCALRawClusterBuilderTemplate");
    ClusterBuilder->Detector("EHCAL");
    ClusterBuilder->SetPlanarGeometry();  // has to be called after Detector()
    ClusterBuilder->Verbosity(verbosity);
    ClusterBuilder->set_threshold_energy(0.100);
    se->registerSubsystem(ClusterBuilder);
    */
  }
  else if (G4EHCAL::EHCAL_clusterizer == G4EHCAL::kEHCALTemplateClusterizer)
  {
    RawClusterBuilderFwd *ClusterBuilder = new RawClusterBuilderFwd("EHCALRawClusterBuilderFwd");
    ClusterBuilder->Detector("EHCAL");
    ClusterBuilder->Verbosity(verbosity);
    ClusterBuilder->set_threshold_energy(0.100);
    se->registerSubsystem(ClusterBuilder);
  }
  else
  {
    cout << "EHCAL_Clusters - unknown clusterizer setting " << G4EHCAL::EHCAL_clusterizer << endl;
    gSystem->Exit(1);
  }

  return;
}

void EHCAL_Eval(const std::string &outputfile)
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::EHCAL_VERBOSITY);
  Fun4AllServer *se = Fun4AllServer::instance();

  CaloEvaluator *eval = new CaloEvaluator("EHCALEVALUATOR", "EHCAL", outputfile.c_str());
  eval->Verbosity(verbosity);
  se->registerSubsystem(eval);

  return;
}
#endif
