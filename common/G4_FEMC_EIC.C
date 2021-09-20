#ifndef MACRO_G4FEMCEIC_C
#define MACRO_G4FEMCEIC_C

#include <GlobalVariables.C>

#include <g4calo/RawTowerBuilderByHitIndex.h>
#include <g4calo/RawTowerDigitizer.h>

#include <g4eiccalos/PHG4ForwardCalCellReco.h>
#include <g4eiccalos/PHG4ForwardEcalSubsystem.h>

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
R__LOAD_LIBRARY(libg4eiccalos.so)
R__LOAD_LIBRARY(libg4eval.so)

namespace Enable
{
  bool FEMC = false;
  bool FEMC_ABSORBER = false;
  bool FEMC_CELL = false;
  bool FEMC_TOWER = false;
  bool FEMC_CLUSTER = false;
  bool FEMC_EVAL = false;
  bool FEMC_OVERLAPCHECK = false;
  int FEMC_VERBOSITY = 0;
}  // namespace Enable

namespace G4FEMC
{
  // from ForwardEcal/mapping/towerMap_FEMC_v007.txt
  const double Gz0 = 310.;
  const double Gdz = 36.5;
  const double outer_radius = 182.655;
  enum enu_Femc_clusterizer
  {
    kFemcGraphClusterizer,
    kFemcTemplateClusterizer
  };
  //template clusterizer, as developed by Sasha Bazilevsky
  enu_Femc_clusterizer Femc_clusterizer = kFemcTemplateClusterizer;
  // graph clusterizer
  //enu_Femc_clusterizer Femc_clusterizer = kFemcGraphClusterizer;
  namespace SETTING
  {
    bool FullEtaAcc = false;
    bool fsPHENIX = false;
    bool EC2x = false;
    bool readoutsplit = true;
    bool asymmetric = true;
    bool wDR = false;
    bool FwdSquare = false;
  }  // namespace SETTING
}  // namespace G4FEMC

void FEMCInit()
{
  // simple way to check if only 1 of the settings is true
  if ((G4FEMC::SETTING::FullEtaAcc ? 1 : 0) + (G4FEMC::SETTING::fsPHENIX ? 1 : 0) + (G4FEMC::SETTING::wDR ? 1 : 0) + (G4FEMC::SETTING::FwdSquare ? 1 : 0) + (G4FEMC::SETTING::asymmetric ? 1 : 0) > 1)
  {
    cout << "use only  G4FHCAL::SETTING::FullEtaAcc=true or G4FHCAL::SETTING::fsPHENIX=true or G4FHCAL::SETTING::wDR=true or G4FHCAL::SETTING::asymmetric=true" << endl;
    gSystem->Exit(1);
  }

  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4FEMC::outer_radius);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, G4FEMC::Gz0 + G4FEMC::Gdz / 2.);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -10.);
}

void FEMCSetup(PHG4Reco *g4Reco)
{
  bool AbsorberActive = Enable::ABSORBER || Enable::FEMC_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::FEMC_OVERLAPCHECK;

  Fun4AllServer *se = Fun4AllServer::instance();

  /** Use dedicated FEMC module */
  PHG4ForwardEcalSubsystem *femc = new PHG4ForwardEcalSubsystem("FEMC");

  ostringstream mapping_femc;

  // PbScint ECAL with nominal eta coverage
  if (G4FEMC::SETTING::FullEtaAcc)
  {
    mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_fullEtaCov.txt";
  }
  // doubled granularity ECAL
  else if (G4FEMC::SETTING::EC2x)
  {
    mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_2x.txt";
  }
  // fsPHENIX ECAL
  else if (G4FEMC::SETTING::fsPHENIX)
  {
    mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_fsPHENIX_v004.txt";
  }
  // asymmetric ECAL around beampipe
  else if (G4FEMC::SETTING::asymmetric)
  {
    if (Enable::IP6){
      if (G4FEMC::SETTING::readoutsplit)
        mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_IP6-asymmetric_ROS.txt";
      else
        mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_IP6-asymmetric.txt";
    } else {
      if (G4FEMC::SETTING::readoutsplit)
        mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_asymmetric_ROS.txt";
      else
        mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_asymmetric.txt";
    }
  }
  // ECAL surrounding dual readout calorimeter
  else if (G4FEMC::SETTING::FwdSquare)
  {
    if (G4FEMC::SETTING::readoutsplit)
      mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_FwdSquare_ROS.txt";
    else
      mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_FwdSquare.txt";
  }
  // ECAL surrounding dual readout calorimeter
  else if (G4FEMC::SETTING::wDR)
  {
    mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_wDR.txt";
  }
  // PbScint ECAL with enlarged beam pipe opening for Mar 2020 beam pipe
  else
  {
    mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_v007.txt";
  }
  cout << mapping_femc.str() << endl;
  femc->SetTowerMappingFile(mapping_femc.str());
  femc->OverlapCheck(OverlapCheck);
  femc->SetActive();
  femc->SetDetailed(false);
  femc->SuperDetector("FEMC");
  if (AbsorberActive) femc->SetAbsorberActive();

  g4Reco->registerSubsystem(femc);
}

void FEMC_Cells()
{
  return;
}

void FEMC_Towers()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::FEMC_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();

  ostringstream mapping_femc;

  //  // fsPHENIX ECAL
  //  mapping_femc << getenv("CALIBRATIONROOT") <<
  //   	"/ForwardEcal/mapping/towerMap_FEMC_fsPHENIX_v004.txt";
  // PbScint ECAL with enlarged beam pipe opening for Mar 2020 beam pipe
  // PbScint ECAL with nominal eta coverage
  if (G4FEMC::SETTING::FullEtaAcc)
  {
    mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_fullEtaCov.txt";
  }
  // doubled granularity ECAL
  else if (G4FEMC::SETTING::EC2x)
  {
    mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_2x.txt";
  }
  // fsPHENIX ECAL
  else if (G4FEMC::SETTING::fsPHENIX)
  {
    mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_fsPHENIX_v004.txt";
  }
  // ECAL surrounding dual readout calorimeter
  else if (G4FEMC::SETTING::FwdSquare)
  {
    mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_FwdSquare.txt";
  }
  // ECAL surrounding dual readout calorimeter
  else if (G4FEMC::SETTING::wDR)
  {
    mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_wDR.txt";
  }
  // asymmetric ECAL around beampipe
  else if (G4FEMC::SETTING::asymmetric)
  {
    if (Enable::IP6){
      if (G4FEMC::SETTING::readoutsplit)
        mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_IP6-asymmetric_ROS.txt";
      else
        mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_IP6-asymmetric.txt";
    } else {
      if (G4FEMC::SETTING::readoutsplit)
        mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_asymmetric_ROS.txt";
      else
        mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_asymmetric.txt";
    }
  }
  // ECAL surrounding dual readout calorimeter
  else if (G4FEMC::SETTING::FwdSquare)
  {
    if (G4FEMC::SETTING::readoutsplit)
      mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_FwdSquare_ROS.txt";
    else
      mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_FwdSquare.txt";
  }
  // PbScint ECAL with enlarged beam pipe opening for Mar 2020 beam pipe
  else
  {
    mapping_femc << getenv("CALIBRATIONROOT") << "/ForwardEcal/mapping/towerMap_FEMC_v007.txt";
  }

  RawTowerBuilderByHitIndex *tower_FEMC = new RawTowerBuilderByHitIndex("TowerBuilder_FEMC");
  tower_FEMC->Detector("FEMC");
  tower_FEMC->set_sim_tower_node_prefix("SIM");
  tower_FEMC->GeometryTableFile(mapping_femc.str());

  se->registerSubsystem(tower_FEMC);

  // PbW crystals
  //RawTowerDigitizer *TowerDigitizer1 = new RawTowerDigitizer("FEMCRawTowerDigitizer1");
  //TowerDigitizer1->Detector("FEMC");
  //TowerDigitizer1->TowerType(1);
  //TowerDigitizer1->Verbosity(verbosity);
  //TowerDigitizer1->set_digi_algorithm(RawTowerDigitizer::kNo_digitization);
  //se->registerSubsystem( TowerDigitizer1 );

  // PbSc towers
  RawTowerDigitizer *TowerDigitizer2 = new RawTowerDigitizer("FEMCRawTowerDigitizer2");
  TowerDigitizer2->Detector("FEMC");
  TowerDigitizer2->TowerType(2);
  TowerDigitizer2->Verbosity(verbosity);
  TowerDigitizer2->set_digi_algorithm(RawTowerDigitizer::kNo_digitization);
  se->registerSubsystem(TowerDigitizer2);

  //  // E864 towers (three types for three sizes)
  //  RawTowerDigitizer *TowerDigitizer3 = new RawTowerDigitizer("FEMCRawTowerDigitizer3");
  //  TowerDigitizer3->Detector("FEMC");
  //  TowerDigitizer3->TowerType(3);
  //  TowerDigitizer3->Verbosity(verbosity);
  //  TowerDigitizer3->set_digi_algorithm(RawTowerDigitizer::kNo_digitization);
  //  se->registerSubsystem( TowerDigitizer3 );
  //
  //  RawTowerDigitizer *TowerDigitizer4 = new RawTowerDigitizer("FEMCRawTowerDigitizer4");
  //  TowerDigitizer4->Detector("FEMC");
  //  TowerDigitizer4->TowerType(4);
  //  TowerDigitizer4->Verbosity(verbosity);
  //  TowerDigitizer4->set_digi_algorithm(RawTowerDigitizer::kNo_digitization);
  //  se->registerSubsystem( TowerDigitizer4 );
  //
  //  RawTowerDigitizer *TowerDigitizer5 = new RawTowerDigitizer("FEMCRawTowerDigitizer5");
  //  TowerDigitizer5->Detector("FEMC");
  //  TowerDigitizer5->TowerType(5);
  //  TowerDigitizer5->Verbosity(verbosity);
  //  TowerDigitizer5->set_digi_algorithm(RawTowerDigitizer::kNo_digitization);
  //  se->registerSubsystem( TowerDigitizer5 );
  //
  //  RawTowerDigitizer *TowerDigitizer6 = new RawTowerDigitizer("FEMCRawTowerDigitizer6");
  //  TowerDigitizer6->Detector("FEMC");
  //  TowerDigitizer6->TowerType(6);
  //  TowerDigitizer6->Verbosity(verbosity);
  //  TowerDigitizer6->set_digi_algorithm(RawTowerDigitizer::kNo_digitization);
  //  se->registerSubsystem( TowerDigitizer6 );

  // PbW crystals
  //RawTowerCalibration *TowerCalibration1 = new RawTowerCalibration("FEMCRawTowerCalibration1");
  //TowerCalibration1->Detector("FEMC");
  //TowerCalibration1->TowerType(1);
  //TowerCalibration1->Verbosity(verbosity);
  //TowerCalibration1->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
  //TowerCalibration1->set_calib_const_GeV_ADC(1.0);  // sampling fraction = 1.0
  //TowerCalibration1->set_pedstal_ADC(0);
  //se->registerSubsystem( TowerCalibration1 );

  // PbSc towers
  RawTowerCalibration *TowerCalibration2 = new RawTowerCalibration("FEMCRawTowerCalibration2");
  TowerCalibration2->Detector("FEMC");
  TowerCalibration2->TowerType(2);
  TowerCalibration2->Verbosity(verbosity);
  TowerCalibration2->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
  if (G4FEMC::SETTING::readoutsplit)
    TowerCalibration2->set_calib_const_GeV_ADC(1.0 / (0.249*0.84));  // sampling fraction = 0.249 for e-
  else 
    TowerCalibration2->set_calib_const_GeV_ADC(1.0 / 0.249);  // sampling fraction = 0.249 for e-
  TowerCalibration2->set_pedstal_ADC(0);
  se->registerSubsystem(TowerCalibration2);

  //  // E864 towers (three types for three sizes)
  //  RawTowerCalibration *TowerCalibration3 = new RawTowerCalibration("FEMCRawTowerCalibration3");
  //  TowerCalibration3->Detector("FEMC");
  //  TowerCalibration3->TowerType(3);
  //  TowerCalibration3->Verbosity(verbosity);
  //  TowerCalibration3->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
  //  TowerCalibration3->set_calib_const_GeV_ADC(1.0/0.030);  // sampling fraction = 0.030
  //  TowerCalibration3->set_pedstal_ADC(0);
  //  se->registerSubsystem( TowerCalibration3 );
  //
  //  RawTowerCalibration *TowerCalibration4 = new RawTowerCalibration("FEMCRawTowerCalibration4");
  //  TowerCalibration4->Detector("FEMC");
  //  TowerCalibration4->TowerType(4);
  //  TowerCalibration4->Verbosity(verbosity);
  //  TowerCalibration4->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
  //  TowerCalibration4->set_calib_const_GeV_ADC(1.0/0.030);  // sampling fraction = 0.030
  //  TowerCalibration4->set_pedstal_ADC(0);
  //  se->registerSubsystem( TowerCalibration4 );
  //
  //  RawTowerCalibration *TowerCalibration5 = new RawTowerCalibration("FEMCRawTowerCalibration5");
  //  TowerCalibration5->Detector("FEMC");
  //  TowerCalibration5->TowerType(5);
  //  TowerCalibration5->Verbosity(verbosity);
  //  TowerCalibration5->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
  //  TowerCalibration5->set_calib_const_GeV_ADC(1.0/0.030);  // sampling fraction = 0.030
  //  TowerCalibration5->set_pedstal_ADC(0);
  //  se->registerSubsystem( TowerCalibration5 );
  //
  //  RawTowerCalibration *TowerCalibration6 = new RawTowerCalibration("FEMCRawTowerCalibration6");
  //  TowerCalibration6->Detector("FEMC");
  //  TowerCalibration6->TowerType(6);
  //  TowerCalibration6->Verbosity(verbosity);
  //  TowerCalibration6->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
  //  TowerCalibration6->set_calib_const_GeV_ADC(1.0/0.030);  // sampling fraction = 0.030
  //  TowerCalibration6->set_pedstal_ADC(0);
  //  se->registerSubsystem( TowerCalibration6 );
}

void FEMC_Clusters()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::FEMC_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();

  if (G4FEMC::Femc_clusterizer == G4FEMC::kFemcTemplateClusterizer)
  {
    RawClusterBuilderHelper *ClusterBuilder = new RawClusterBuilderkMA("FEMCRawClusterBuilderkMA");
    ClusterBuilder->Detector("FEMC");
    ClusterBuilder->set_seed_e(0.1);
    ClusterBuilder->set_agg_e(0.005);
    se->registerSubsystem(ClusterBuilder);
    /*
    RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplateFEMC");
    ClusterBuilder->Detector("FEMC");
    ClusterBuilder->Verbosity(verbosity);
    ClusterBuilder->set_threshold_energy(0.020);  // This threshold should be the same as in FEMCprof_Thresh**.root file below
    std::string femc_prof = getenv("CALIBRATIONROOT");
    femc_prof += "/EmcProfile/FEMCprof_Thresh20MeV.root";
    ClusterBuilder->LoadProfile(femc_prof.c_str());
    se->registerSubsystem(ClusterBuilder);
    */
  }
  else if (G4FEMC::Femc_clusterizer == G4FEMC::kFemcGraphClusterizer)
  {
    RawClusterBuilderFwd *ClusterBuilder = new RawClusterBuilderFwd("FEMCRawClusterBuilderFwd");

    ClusterBuilder->Detector("FEMC");
    ClusterBuilder->Verbosity(verbosity);
    ClusterBuilder->set_threshold_energy(0.010);
    se->registerSubsystem(ClusterBuilder);
  }
  else
  {
    cout << "FEMC_Clusters - unknown clusterizer setting!" << endl;
    exit(1);
  }

  return;
}

void FEMC_Eval(const std::string &outputfile)
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::FEMC_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();

  CaloEvaluator *eval = new CaloEvaluator("FEMCEVALUATOR", "FEMC", outputfile);
  eval->Verbosity(verbosity);
  se->registerSubsystem(eval);

  return;
}
#endif
