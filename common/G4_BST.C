#ifndef MACRO_G4BST_C
#define MACRO_G4BST_C

#include <GlobalVariables.C>

// #include <g4calo/RawTowerDigitizer.h>

// #include <g4eiccalos/PHG4ForwardCalCellReco.h>
#include <bst/PHG4BSTSubsystem.h>

// #include <g4eval/CaloEvaluator.h>

#include <g4main/PHG4Reco.h>


#include <fun4all/Fun4AllServer.h>

R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libg4calo.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4eval.so)

namespace Enable
{
  bool BST = false;
  bool BST_ABSORBER = false;
  bool BST_CELL = false;
  bool BST_TOWER = false;
  bool BST_CLUSTER = false;
  bool BST_EVAL = false;
  bool BST_OVERLAPCHECK = false;
  int BST_VERBOSITY = 0;
}  // namespace Enable

namespace G4BST
{
  double Gz0 = 0.;
  double Gdz = 100.;
  double outer_radius = 60.;
  namespace SETTING
  {
    bool Tungsten = false;
  }  // namespace SETTING
}  // namespace G4BST

void BSTInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4BST::outer_radius);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, G4BST::Gz0 + G4BST::Gdz / 2.);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, G4BST::Gz0 - G4BST::Gdz / 2.);
}

void BSTSetup(PHG4Reco *g4Reco)
{
  const bool AbsorberActive = Enable::ABSORBER || Enable::BST_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::BST_OVERLAPCHECK;
  Fun4AllServer *se = Fun4AllServer::instance();

  /** Use dedicated BST module */
  PHG4BSTSubsystem *hhcal = new PHG4BSTSubsystem("BST");
  // hhcal->OverlapCheck(OverlapCheck);
  hhcal->OverlapCheck(true);
  hhcal->SuperDetector("BST");
  hhcal->SetActive();
  if (AbsorberActive) hhcal->SetAbsorberActive();
  g4Reco->registerSubsystem(hhcal);
}

void BST_Cells(int verbosity = 0)
{
  return;
}

void BST_Towers()
{
  return;
}

void BST_Clusters()
{
  return;
}

void BST_Eval(const std::string &outputfile)
{
  // int verbosity = std::max(Enable::VERBOSITY, Enable::BST_VERBOSITY);
  // Fun4AllServer *se = Fun4AllServer::instance();

  // CaloEvaluator *eval = new CaloEvaluator("BSTEVALUATOR", "BST", outputfile.c_str());
  // eval->Verbosity(verbosity);
  // se->registerSubsystem(eval);

  return;
}
#endif
