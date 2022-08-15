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
    int optionMat = 0;
    int BSTnoSupport = 0;
    int bent_sagitta_default = 0;
    int bent_sagitta_mod = 0;
    int ECCE_with_OuterStave = 0;
    int use_OuterStave_lowerR = 0;
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
  PHG4BSTSubsystem *hbst = new PHG4BSTSubsystem("BST");
  hbst->OverlapCheck(OverlapCheck);
  // hbst->OverlapCheck(true);
  cout << "BST material special setting " << G4BST::SETTING::optionMat << endl;
  if(G4BST::SETTING::BSTnoSupport==1){
    hbst->set_int_param("do_internal_supports", 0);
    hbst->set_int_param("do_external_supports", 0);
  }
  if(Enable::EPIC_TRACKINGGEO){
    hbst->set_int_param("use_EPIC_setup", 1);
    if(G4BST::SETTING::bent_sagitta_default){
      hbst->set_int_param("use_bent_wafer_sagittas_default", 1);
    } else if(G4BST::SETTING::bent_sagitta_mod){
      hbst->set_int_param("use_bent_wafer_sagittas_mod", 1);
    } else if(G4BST::SETTING::ECCE_with_OuterStave){
      hbst->set_int_param("use_ECCE_with_OuterStave", 1);
      if(G4BST::SETTING::use_OuterStave_lowerR){
        hbst->set_double_param("radius_outer_stave", 33.);
      }
    }
  }
  if(G4BST::SETTING::optionMat==1){
    hbst->set_double_param("layer_backing_thickness", 0.035 / 100 * 28.57 * cm); // 100 microns of Kapton
  }
  else if(G4BST::SETTING::optionMat==2){
    hbst->set_double_param("layer_backing_thickness", 0.07 / 100 * 28.57 * cm); // 200 microns of Kapton
  }
  else if(G4BST::SETTING::optionMat==3){
    hbst->set_double_param("layer_backing_thickness", 0.14 / 100 * 28.57 * cm); // 400 microns of Kapton
  }
  hbst->SuperDetector("BST");
  hbst->SetActive();
  if (AbsorberActive) hbst->SetAbsorberActive();
  g4Reco->registerSubsystem(hbst);

  int ilayersBST = 5;
  // TODO FIX currently hardcoded values!
  G4double layerradii[] = {
    36.16,
    48.2239,
    60.1924,
    198.307,
    210.276
  };
  G4double layerradii_ECCE_wStave[] = {
    36.16,
    48.2239,
    60.1924,
    198.307,
    210.276,
    420.0
  };
  G4double layerradii_EPIC[] = {
    36.16,
    48.2239,
    120.070,
    270.0,
    420.0,
    0.0
  };
  G4double layerradii_EPIC_sagdef[] = {
    36.16,
    48.2239,
    120.070,
    198.307,
    210.276,
    420.0
  };
  G4double layerradii_EPIC_sagmod[] = {
    36.16,
    48.2239,
    120.070,
    240.3,
    270.0,
    420.0
  };

  if(Enable::EPIC_TRACKINGGEO){
    if(G4BST::SETTING::bent_sagitta_default){
      ilayersBST = 6;
    } else if(G4BST::SETTING::bent_sagitta_mod){
      ilayersBST = 6;
    } else if(G4BST::SETTING::ECCE_with_OuterStave){
      ilayersBST = 6;
    }
  }

  for(int ilay=0; ilay<ilayersBST; ilay++){
    if(Enable::EPIC_TRACKINGGEO){
      if(G4BST::SETTING::bent_sagitta_default){
        layerradii_EPIC[ilay] = layerradii_EPIC_sagdef[ilay];
      } else if(G4BST::SETTING::bent_sagitta_mod){
        layerradii_EPIC[ilay] = layerradii_EPIC_sagmod[ilay];
      } else if(G4BST::SETTING::ECCE_with_OuterStave){
        layerradii_EPIC[ilay] = layerradii_ECCE_wStave[ilay];
      }
    }
    if (TRACKING::FastKalmanFilter)
    {
      TRACKING::FastKalmanFilter->add_phg4hits(string(Form("G4HIT_BST_%d",ilay)),     //      const std::string& phg4hitsNames,
                                              PHG4TrackFastSim::Cylinder,
                                              999.,                      // radial-resolution [cm]
                                              10. / 10000. / sqrt(12.),  // azimuthal-resolution [cm]
                                              10. / 10000. / sqrt(12.),  // z-resolution [cm]
                                              0.95,                         // efficiency,
                                              0                          // noise hits
      );
      TRACKING::FastKalmanFilter->add_cylinder_state(Form("BST_%d",ilay), Enable::EPIC_TRACKINGGEO ? layerradii_EPIC[ilay] : layerradii[ilay]);
      TRACKING::ProjectionNames.insert(Form("BST_%d",ilay));
    }
    if (TRACKING::FastKalmanFilterInnerTrack)
    {
      TRACKING::FastKalmanFilterInnerTrack->add_phg4hits(string(Form("G4HIT_BST_%d",ilay)),     //      const std::string& phg4hitsNames,
                                              PHG4TrackFastSim::Cylinder,
                                              999.,                      // radial-resolution [cm]
                                              10. / 10000. / sqrt(12.),  // azimuthal-resolution [cm]
                                              10. / 10000. / sqrt(12.),  // z-resolution [cm]
                                              0.95,                         // efficiency,
                                              0                          // noise hits
      );
      TRACKING::FastKalmanFilterInnerTrack->add_cylinder_state(Form("BST_%d",ilay), Enable::EPIC_TRACKINGGEO ? layerradii_EPIC[ilay] : layerradii[ilay]);
    }
  }
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
