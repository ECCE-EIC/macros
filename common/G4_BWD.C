#ifndef MACRO_G4BWD_C
#define MACRO_G4BWD_C

#include <GlobalVariables.C>

//include our own Bwd Raw Tower Builder
#include <eicg4bwd/BwdRawTowerBuilderByHitIndex.h>

#include <g4calo/RawTowerDigitizer.h>

#include <g4eiccalos/PHG4ForwardCalCellReco.h>
// Use Forward Cal Cell Reco .

#include <eicg4bwd/EICG4BwdSubsystem.h>
// Include our Subsystem

//Standard RawTowerDefs.h is modified

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
  bool BWD = false;
  bool BWDN[5]={true,false,false,false,false};
  bool BWD_ABSORBER = false;
  bool BWD_CELL = false;
  bool BWD_TOWER = false;
  bool BWD_CLUSTER = false;
  bool BWD_EVAL = false;
  bool BWD_OVERLAPCHECK = false;
  int BWD_VERBOSITY = 0;
}  // namespace Enable


namespace G4BWD
{
  
  string mapname[5]={"BWD_mapping_v1.txt","BWD_mapping_v2.txt","BWD_mapping_v3.txt","BWD_mapping_v4.txt","BWD_mapping_v5.txt"};
  double minz = -5000;
  double maxz = -300;
  double radius =  100;

  // Default set to B0 Ecal position at IP6

  // Digitization (default photon digi):
  RawTowerDigitizer::enu_digi_algorithm TowerDigi = RawTowerDigitizer::kNo_digitization;
  // directly pass the energy of sim tower to digitized tower
  // kNo_digitization
  // simple digitization with photon statistics, single amplitude ADC conversion and pedestal
  // kSimple_photon_digitization
  // digitization with photon statistics on SiPM with an effective pixel N, ADC conversion and pedestal
  // kSiPM_photon_digitization

}  // namespace G4B0ECAL

void BWDInit()
{
}

void BWDSetup(PHG4Reco *g4Reco)
{
//Done in G4_hFarBwdBeamLine.C
}

void BWD_Cells(int verbosity = 0)
{
  return;
}

void BWD_Towers()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::BWD_VERBOSITY);


for (int i = 0; i < 5; i++){
	if (!Enable::BWDN[i])continue;
  Fun4AllServer *se = Fun4AllServer::instance();
  ostringstream mapping_bwd;
  mapping_bwd << getenv("CALIBRATIONROOT") << "/BWD/mapping/"<<G4BWD::mapname[i];
  //mapping_bwd << G4BWD::mapname[i];
  BwdRawTowerBuilderByHitIndex *tower_BWD = new BwdRawTowerBuilderByHitIndex(Form("TowerBuilder_BWD_%d", i));
  tower_BWD->Detector(Form("BWD_%d", i));
  tower_BWD->set_sim_tower_node_prefix("SIM");
  tower_BWD->GeometryTableFile(mapping_bwd.str());

  se->registerSubsystem(tower_BWD);

  
  RawTowerDigitizer *TowerDigitizer = new RawTowerDigitizer(Form("BWDRawTowerDigitizer_%d",i));
  TowerDigitizer->Detector(Form("BWD_%d", i));
  TowerDigitizer->Verbosity(verbosity);
  TowerDigitizer->set_digi_algorithm(RawTowerDigitizer::kNo_digitization);
  se->registerSubsystem(TowerDigitizer);

  RawTowerCalibration *TowerCalibration = new RawTowerCalibration(Form("BWDRawTowerCalibration_%d",i));
  TowerCalibration->Detector(Form("BWD_%d", i));
  TowerCalibration->Verbosity(verbosity);
  TowerCalibration->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
  TowerCalibration->set_calib_const_GeV_ADC(1. ); 
  TowerCalibration->set_pedstal_ADC(0);
  se->registerSubsystem(TowerCalibration);
}
}

void BWD_Clusters()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::BWD_VERBOSITY);

for (int i = 0; i < 5; i++){
	if (!Enable::BWDN[i])continue;
  Fun4AllServer *se = Fun4AllServer::instance();
  RawClusterBuilderFwd *ClusterBuilder = new RawClusterBuilderFwd(Form("BWDRawClusterBuilderFwd_%d",i));
  ClusterBuilder->Detector(Form("BWD_%d", i));
  ClusterBuilder->Verbosity(verbosity);
  ClusterBuilder->set_threshold_energy(0.100);
  se->registerSubsystem(ClusterBuilder);
}
  return;
}

void BWD_Eval(const std::string &outputfile)
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::BWD_VERBOSITY);
  string filename=outputfile.c_str();
for (int i = 0; i < 5; i++){
	if (!Enable::BWDN[i])continue;
  Fun4AllServer *se = Fun4AllServer::instance();
  CaloEvaluator *eval = new CaloEvaluator(Form("BWDEVALUATOR_%d",i), Form("BWD_%d", i), (filename+Form("_%d.root",i)) );
  eval->Verbosity(verbosity);
  se->registerSubsystem(eval);
}
  return;
}
#endif
