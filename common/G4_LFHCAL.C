#ifndef MACRO_G4LFHCAL_C
#define MACRO_G4LFHCAL_C

#include <GlobalVariables.C>

#include <g4calo/RawTowerDigitizer.h>

#include <g4eiccalos/PHG4LFHcalSubsystem.h>
#include <g4eiccalos/RawTowerBuilderByHitIndexLHCal.h>

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
  bool LFHCAL = false;
  bool LFHCAL_ABSORBER = false;
  bool LFHCAL_CELL = false;
  bool LFHCAL_TOWER = false;
  bool LFHCAL_CLUSTER = false;
  bool LFHCAL_EVAL = false;
  bool LFHCAL_OVERLAPCHECK = false;
  int LFHCAL_VERBOSITY = 0;
}  // namespace Enable

namespace G4LFHCAL
{
  // from LFHcal/mapping/towerMap_LFHCAL_v005.txt
  double Gz0 = 400.;
  double Gdz = 100.;
  double outer_radius = 265.;
  enum enu_FHcal_clusterizer
  {
    kFHcalGraphClusterizer,
    kFHcalTemplateClusterizer
  };
  //template clusterizer, as developed by Sasha Bazilevsky
  enu_FHcal_clusterizer FHcal_clusterizer = kFHcalTemplateClusterizer;
  // graph clusterizer
  //enu_FHcal_clusterizer FHcal_clusterizer = kFHcalGraphClusterizer;
  namespace SETTING
  {
    bool FullEtaAcc   = false;
    bool HC2x         = false;
    bool asymmetric   = true;
    bool wDR          = false;
    bool FwdSquare    = false;
    bool FwdConfig    = false;
    bool longer       = true;
    bool tailcatcher  = true;
  }  // namespace SETTING
}  // namespace G4LFHCAL

TString GetMappingFile(){
  TString mappinFileName = getenv("CALIBRATIONROOT");
  if (G4LFHCAL::SETTING::HC2x )
  {
    if (G4LFHCAL::SETTING::longer)
      mappinFileName += "/LFHcal/mapping/towerMap_LFHCAL_2x-long.txt";
    else 
      mappinFileName +=  "/LFHcal/mapping/towerMap_LFHCAL_2x.txt";
  }
  // HCal Fe-Scint surrounding dual readout calorimeter R>50cm
  else if (G4LFHCAL::SETTING::wDR)
  {
    if (G4LFHCAL::SETTING::longer)
      mappinFileName +=  "/LFHcal/mapping/towerMap_LFHCAL_wDR-long.txt";
    else 
      mappinFileName +=  "/LFHcal/mapping/towerMap_LFHCAL_wDR.txt";
  }
  // HCal Fe-Scint surrounding dual readout calorimeter R>50cm
  else if (G4LFHCAL::SETTING::FwdConfig)
  {
    if (G4LFHCAL::SETTING::longer)
      mappinFileName += "/LFHcal/mapping/towerMap_LFHCAL_FwdConfig-long.txt";
    else 
      mappinFileName += "/LFHcal/mapping/towerMap_LFHCAL_FwdConfig.txt";
  }
  // HCal Fe-Scint surrounding dual readout calorimeter R>50cm
  else if (G4LFHCAL::SETTING::FwdSquare)
  {
    if (G4LFHCAL::SETTING::longer)
      mappinFileName += "/LFHcal/mapping/towerMap_LFHCAL_FwdSquare-long.txt";
    else 
      mappinFileName += "/LFHcal/mapping/towerMap_LFHCAL_FwdSquare.txt";
  }
  // full HCal Fe-Scint with asymmetric centering around beampipe
  else if (G4LFHCAL::SETTING::asymmetric)
  {
    if (Enable::IP6){
      if (G4LFHCAL::SETTING::longer){
        if (G4LFHCAL::SETTING::tailcatcher)
          mappinFileName += "/LFHcal/mapping/towerMap_LFHCAL_IP6-asymmetric-long-tailcatcher.txt";
        else
          mappinFileName += "/LFHcal/mapping/towerMap_LFHCAL_IP6-asymmetric-long.txt";
      } else 
        mappinFileName += "/LFHcal/mapping/towerMap_LFHCAL_IP6-asymmetric.txt";
    } else {
      if (G4LFHCAL::SETTING::longer)
        mappinFileName += "/LFHcal/mapping/towerMap_LFHCAL_asymmetric-long.txt";
      else 
        mappinFileName += "/LFHcal/mapping/towerMap_LFHCAL_asymmetric.txt";      
    }
  }
  //  PSD like HCal Fe-Scint with enlarged beam pipe opening for Mar 2020 beam pipe
  else
  {
    if (G4LFHCAL::SETTING::longer)
      mappinFileName +=  "/LFHcal/mapping/towerMap_LFHCAL_default-long.txt";
    else
      mappinFileName +=  "/LFHcal/mapping/towerMap_LFHCAL_default.txt";
  }

  return mappinFileName;
  
}


void LFHCALInit()
{
  if (G4LFHCAL::SETTING::longer){
    G4LFHCAL::Gz0 = 420;
    G4LFHCAL::Gdz = 140;
  }
  // simple way to check if only 1 of the settings is true
  if ((G4LFHCAL::SETTING::FullEtaAcc ? 1 : 0) + (G4LFHCAL::SETTING::HC2x ? 1 : 0) > 1)
  {
    cout << "use only  G4LFHCAL::SETTING::FullEtaAcc=true or G4LFHCAL::SETTING::HC2x=true or G4LFHCAL::SETTING::HC4x=true" << endl;
    gSystem->Exit(1);
  }
 
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4LFHCAL::outer_radius);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, G4LFHCAL::Gz0 + G4LFHCAL::Gdz / 2.);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -10*cm);
}

void LFHCALSetup(PHG4Reco *g4Reco)
{
  const bool AbsorberActive = Enable::ABSORBER || Enable::LFHCAL_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::LFHCAL_OVERLAPCHECK;
  Fun4AllServer *se = Fun4AllServer::instance();

  /** Use dedicated LFHCAL module */
  PHG4LFHcalSubsystem *fhcal = new PHG4LFHcalSubsystem("LFHCAL");

  TString mapping_fhcal = GetMappingFile();
  cout << "LFHCAL: "<< mapping_fhcal.Data() << endl;
  ostringstream mapping_fhcal_s;
  mapping_fhcal_s << mapping_fhcal.Data();
  
  fhcal->SetTowerMappingFile(mapping_fhcal_s.str());
  fhcal->OverlapCheck(OverlapCheck);
  fhcal->SetActive();
  //fhcal->SetDetailed(true);
  fhcal->SuperDetector("LFHCAL");
  if (AbsorberActive) fhcal->SetAbsorberActive();

  g4Reco->registerSubsystem(fhcal);
}

void LFHCAL_Cells(int verbosity = 0)
{
  return;
}

void LFHCAL_Towers()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::LFHCAL_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();

  // Switch to desired calo setup;
  // PSD like HCal Fe-Scint with doubled granularity  
  TString mapping_fhcal = GetMappingFile();
  ostringstream mapping_fhcal_s;
  mapping_fhcal_s << mapping_fhcal.Data();

  RawTowerBuilderByHitIndexLHCal *tower_LFHCAL = new RawTowerBuilderByHitIndexLHCal("TowerBuilder_LFHCAL");
  tower_LFHCAL->Detector("LFHCAL");
  tower_LFHCAL->set_sim_tower_node_prefix("SIM");
  tower_LFHCAL->GeometryTableFile(mapping_fhcal_s.str());

  se->registerSubsystem(tower_LFHCAL);

  cout << "def: using default for LFHCAL towers" << endl;
  RawTowerDigitizer *TowerDigitizer = new RawTowerDigitizer("LFHCALRawTowerDigitizer");
  TowerDigitizer->Detector("LFHCAL");
  TowerDigitizer->Verbosity(verbosity);
  TowerDigitizer->set_digi_algorithm(RawTowerDigitizer::kNo_digitization);
  se->registerSubsystem(TowerDigitizer);

  
  RawTowerCalibration *TowerCalibration = new RawTowerCalibration("LFHCALRawTowerCalibration");
  TowerCalibration->Detector("LFHCAL");
  TowerCalibration->Verbosity(verbosity);
  TowerCalibration->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
  TowerCalibration->set_calib_const_GeV_ADC(1. / (0.03898*0.93));  // temporary factor 0.93 to fix calibration for new tower design
  TowerCalibration->set_pedstal_ADC(0);
  se->registerSubsystem(TowerCalibration);
}

void LFHCAL_Clusters()
{
  Fun4AllServer *se = Fun4AllServer::instance();

  RawClusterBuilderHelper *ClusterBuilder = new RawClusterBuilderkMA("LFHCALRawClusterBuilderkMA");
  ClusterBuilder->Detector("LFHCAL");
  ClusterBuilder->set_seed_e(0.1);
  ClusterBuilder->set_agg_e(0.001);
  se->registerSubsystem(ClusterBuilder);

  return;
}

void LFHCAL_Eval(const std::string &outputfile)
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::LFHCAL_VERBOSITY);
  Fun4AllServer *se = Fun4AllServer::instance();

  CaloEvaluator *eval = new CaloEvaluator("LFHCALEVALUATOR", "LFHCAL", outputfile.c_str());
  eval->Verbosity(verbosity);
  se->registerSubsystem(eval);

  return;
}
#endif
