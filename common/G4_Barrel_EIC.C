#ifndef MACRO_G4BARREL_C
#define MACRO_G4BARREL_C

#include <g4lblvtx/AllSiliconTrackerSubsystem.h>
#include <g4lblvtx/SimpleNtuple.h>
#include <g4main/PHG4Reco.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include "GlobalVariables.C"

#include <string>

R__LOAD_LIBRARY(libg4detectors.so)

using namespace std;

namespace Enable
{
  bool BARREL = false;
  bool BARREL_OVERLAPCHECK = false;
}  // namespace Enable

namespace G4BARRELEIC
{
  namespace SETTING
  {
    float SAGITTAX0 = 0.05;
    float TRACKING_EFFICIENCY = 1.00;
  }  // namespace SETTING
}  // namespace G4BARRELEIC

//-----------------------------------------------------------------------------------//
void BarrelInit()
{
}

void BarrelFastKalmanFilterConfigSVTX(PHG4TrackFastSim * kalman_filter, int ilay, double radius, double silicon_thickness, bool addproj)
{

  // import Kalman filter config (lines 226 to 246 here: https://github.com/eic/g4lblvtx/blob/master/macros/auxiliary_studies/simplified_geometry/Fun4All_G4_simplified_v2.C):

  // add Vertexing Layers
  kalman_filter->add_phg4hits(
      Form("G4HIT_SVTX_%d",ilay),  // const std::string& phg4hitsNames,
      PHG4TrackFastSim::Cylinder,
      silicon_thickness / sqrt(12.),                      // radial-resolution [cm]
      10. / 10000. / sqrt(12.),  // azimuthal-resolution [cm]
      10. / 10000. / sqrt(12.),  // z-resolution [cm]
      G4BARRELEIC::SETTING::TRACKING_EFFICIENCY,                         // efficiency,
      0                          // noise hits
  );
  kalman_filter->add_cylinder_state(Form("SVTX_%d",ilay), radius);
  if(addproj)TRACKING::ProjectionNames.insert(Form("SVTX_%d",ilay));
}

void BarrelFastKalmanFilterConfigBARR(PHG4TrackFastSim * kalman_filter, int ilay, double radius, double silicon_thickness, bool addproj)
{
  // add Barrel Layers
  kalman_filter->add_phg4hits(
      Form("G4HIT_BARR_%d",ilay),  // const std::string& phg4hitsNames,
      PHG4TrackFastSim::Cylinder,
      silicon_thickness / sqrt(12.),                      // radial-resolution [cm]
      10. / 10000. / sqrt(12.),  // azimuthal-resolution [cm]
      10. / 10000. / sqrt(12.),  // z-resolution [cm]
      G4BARRELEIC::SETTING::TRACKING_EFFICIENCY,                         // efficiency,
      0                          // noise hits
  );
  kalman_filter->add_cylinder_state(Form("BARR_%d",ilay), radius);
  if(addproj)TRACKING::ProjectionNames.insert(Form("BARR_%d",ilay));

}


//-----------------------------------------------------------------------------------//
void Barrel(PHG4Reco *g4Reco, int det_ver = 3)
{
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::BARREL_OVERLAPCHECK;
  // import Geometry (lines 111 to 148 in https://github.com/eic/g4lblvtx/blob/master/macros/auxiliary_studies/simplified_geometry/Fun4All_G4_simplified_v2.C):

  cout << "BARREL: Using tracking efficiency of: " << G4BARRELEIC::SETTING::TRACKING_EFFICIENCY << endl;

  PHG4CylinderSubsystem *cyl(nullptr);

  if(Enable::AI_TRACKINGGEO){
    //---------------------------
    // Vertexing
    double si_vtx_r_pos[] = {3.40, 5.67, 7.93};
    const int nVtxLayers = sizeof(si_vtx_r_pos) / sizeof(*si_vtx_r_pos);
    double silicon_thickness = 0.05 / 100. * 9.37;
    for (int ilayer = 0; ilayer < nVtxLayers; ilayer++)
    {
      cyl = new PHG4CylinderSubsystem("SVTX", ilayer);
      cyl->set_string_param("material", "G4_Si");
      cyl->set_double_param("radius", si_vtx_r_pos[ilayer]);
      cyl->set_double_param("thickness", silicon_thickness);
      cyl->set_double_param("place_z", 0.);
      cyl->set_double_param("length", 30.);
      cyl->SetActive();
  //    cyl->SuperDetector("SVTX");  breakout SVTX into individual layers
      cyl->OverlapCheck(OverlapCheck);
      g4Reco->registerSubsystem(cyl);
      BarrelFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilter, ilayer, si_vtx_r_pos[ilayer],silicon_thickness,true);
      BarrelFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilterInnerTrack, ilayer, si_vtx_r_pos[ilayer],silicon_thickness,false);
      BarrelFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilterSiliconTrack, ilayer, si_vtx_r_pos[ilayer],silicon_thickness,false);
    }
  } else if(Enable::EPIC_TRACKINGGEO){
    //---------------------------
    // Vertexing
    // 270 mm vertex layer with X/X0 ~ 0.05% at r = 36mm
    //   270 mm vertex layer with X/X0 ~ 0.05% at r = 48mm
    //   270 mm vertex layer with X/X0 ~ 0.05% at r = 120mm
    double si_vtx_r_pos[] = {3.60, 4.8, 12.00};
    const int nVtxLayers = sizeof(si_vtx_r_pos) / sizeof(*si_vtx_r_pos);
    double silicon_thickness = 0.05 / 100. * 9.37;
    for (int ilayer = 0; ilayer < nVtxLayers; ilayer++)
    {
      cyl = new PHG4CylinderSubsystem("SVTX", ilayer);
      cyl->set_string_param("material", "G4_Si");
      cyl->set_double_param("radius", si_vtx_r_pos[ilayer]);
      cyl->set_double_param("thickness", silicon_thickness);
      cyl->set_double_param("place_z", 0.);
      cyl->set_double_param("length", 27.);
      cyl->SetActive();
  //    cyl->SuperDetector("SVTX");  breakout SVTX into individual layers
      cyl->OverlapCheck(OverlapCheck);
      g4Reco->registerSubsystem(cyl);
      BarrelFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilter, ilayer, si_vtx_r_pos[ilayer],silicon_thickness,true);
      BarrelFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilterInnerTrack, ilayer, si_vtx_r_pos[ilayer],silicon_thickness,false);
      BarrelFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilterSiliconTrack, ilayer, si_vtx_r_pos[ilayer],silicon_thickness,false);
    }
  } else {
    //---------------------------
    // Vertexing
    double si_vtx_r_pos[] = {3.3, 4.35, 5.4};
    const int nVtxLayers = sizeof(si_vtx_r_pos) / sizeof(*si_vtx_r_pos);
    double silicon_thickness = 0.05 / 100. * 9.37;
    for (int ilayer = 0; ilayer < nVtxLayers; ilayer++)
    {
      cyl = new PHG4CylinderSubsystem("SVTX", ilayer);
      cyl->set_string_param("material", "G4_Si");
      cyl->set_double_param("radius", si_vtx_r_pos[ilayer]);
      cyl->set_double_param("thickness", silicon_thickness);
      cyl->set_double_param("place_z", 0.);
      cyl->set_double_param("length", 27);
      cyl->SetActive();
  //    cyl->SuperDetector("SVTX");  breakout SVTX into individual layers
      cyl->OverlapCheck(OverlapCheck);
      g4Reco->registerSubsystem(cyl);

      BarrelFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilter, ilayer, si_vtx_r_pos[ilayer],silicon_thickness,true);
      BarrelFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilterInnerTrack, ilayer, si_vtx_r_pos[ilayer],silicon_thickness,false);
      BarrelFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilterSiliconTrack, ilayer, si_vtx_r_pos[ilayer],silicon_thickness,false);
    }
  }


  //---------------------------
  // Barrel
  if(Enable::AI_TRACKINGGEO){
    double si_len = 60.0;
    double z_e_length[] = {30., 30.};
    double z_h_length[] = {30., 30.};
    double si_z_length[] = {si_len, si_len};
    double si_r_pos[] = {15.3, 17.0}; // Modified on 29th Oct to account for new struture design
    const int nTrckLayers = sizeof(si_r_pos) / sizeof(*si_r_pos);
    double silicon_thickness = G4BARRELEIC::SETTING::SAGITTAX0 / 100. * 9.37;
    
    for (int ilayer = 0; ilayer < nTrckLayers; ilayer++)
    {
      //cout << "Radius " << ilayer + 1 << ": " << si_r_pos[ilayer] << "cms \t e-length : " << z_e_length[ilayer] << "cms \t h-length : " << z_h_length[ilayer] << "cms"<< endl;
      //cout << "eslope : " << e_slope1 << " \n hslope : " << h_slope1 << " \n zleft : " << z_e_length[ilayer] << "\n zright : " << z_h_length[ilayer] << " \n radius : " << si_r_pos[ilayer] << endl ;
      cyl = new PHG4CylinderSubsystem("BARR", ilayer);
      cyl->set_string_param("material", "G4_Si");
      cyl->set_double_param("radius", si_r_pos[ilayer]);
      cyl->set_double_param("thickness", silicon_thickness);
      cyl->set_double_param("place_z", (z_h_length[ilayer] - z_e_length[ilayer])/2);
      cyl->set_double_param("length", si_z_length[ilayer]);
      cyl->SetActive();
      cyl->OverlapCheck(OverlapCheck);
  //    cyl->SuperDetector("BARR");   breakout BARR into individual layers
      g4Reco->registerSubsystem(cyl);

      BarrelFastKalmanFilterConfigBARR(TRACKING::FastKalmanFilter, ilayer, si_r_pos[ilayer],silicon_thickness,true);
      BarrelFastKalmanFilterConfigBARR(TRACKING::FastKalmanFilterInnerTrack, ilayer, si_r_pos[ilayer],silicon_thickness,false);
      BarrelFastKalmanFilterConfigBARR(TRACKING::FastKalmanFilterSiliconTrack, ilayer, si_r_pos[ilayer],silicon_thickness,false);
    }
  } else if(Enable::EPIC_TRACKINGGEO){
    double si_z_length[] = {54.0, 84.0};
    double si_r_pos[] = {27.0, 42.0}; // Modified on 29th Oct to account for new struture design
    const int nTrckLayers = sizeof(si_r_pos) / sizeof(*si_r_pos);
    float sagittaX0layer[] = {0.25, 0.55};
    
    for (int ilayer = 0; ilayer < nTrckLayers; ilayer++)
    {
      double silicon_thickness = sagittaX0layer[ilayer] / 100. * 9.37;
      //cout << "Radius " << ilayer + 1 << ": " << si_r_pos[ilayer] << "cms \t e-length : " << z_e_length[ilayer] << "cms \t h-length : " << z_h_length[ilayer] << "cms"<< endl;
      //cout << "eslope : " << e_slope1 << " \n hslope : " << h_slope1 << " \n zleft : " << z_e_length[ilayer] << "\n zright : " << z_h_length[ilayer] << " \n radius : " << si_r_pos[ilayer] << endl ;
      cyl = new PHG4CylinderSubsystem("BARR", ilayer);
      cyl->set_string_param("material", "G4_Si");
      cyl->set_double_param("radius", si_r_pos[ilayer]);
      cyl->set_double_param("thickness", silicon_thickness);
      cyl->set_double_param("place_z", 0.0);
      cyl->set_double_param("length", si_z_length[ilayer]);
      cyl->SetActive();
      cyl->OverlapCheck(OverlapCheck);
  //    cyl->SuperDetector("BARR");   breakout BARR into individual layers
      g4Reco->registerSubsystem(cyl);

      BarrelFastKalmanFilterConfigBARR(TRACKING::FastKalmanFilter, ilayer, si_r_pos[ilayer],silicon_thickness,true);
      BarrelFastKalmanFilterConfigBARR(TRACKING::FastKalmanFilterInnerTrack, ilayer, si_r_pos[ilayer],silicon_thickness,false);
      BarrelFastKalmanFilterConfigBARR(TRACKING::FastKalmanFilterSiliconTrack, ilayer, si_r_pos[ilayer],silicon_thickness,false);
    }
  } else {
    double z_e_length[] = {-27, -29.0};
    double z_h_length[] = {27., 29.0};
    double si_r_pos[] = {21, 22.68}; // Modified on 29th Oct to account for new struture design
    const int nTrckLayers = sizeof(si_r_pos) / sizeof(*si_r_pos);
    double silicon_thickness = G4BARRELEIC::SETTING::SAGITTAX0 / 100. * 9.37;
    
    for (int ilayer = 0; ilayer < nTrckLayers; ilayer++)
    {
      cyl = new PHG4CylinderSubsystem("BARR", ilayer);
      cyl->set_string_param("material", "G4_Si");
      cyl->set_double_param("radius", si_r_pos[ilayer]);
      cyl->set_double_param("thickness", silicon_thickness);
      cyl->set_double_param("place_z", (z_h_length[ilayer] + z_e_length[ilayer])/2);
      cyl->set_double_param("length", (z_h_length[ilayer] - z_e_length[ilayer]));
      cyl->SetActive();
      cyl->OverlapCheck(OverlapCheck);
  //    cyl->SuperDetector("BARR");   breakout BARR into individual layers
      g4Reco->registerSubsystem(cyl);

      BarrelFastKalmanFilterConfigBARR(TRACKING::FastKalmanFilter, ilayer, si_r_pos[ilayer],silicon_thickness,true);
      BarrelFastKalmanFilterConfigBARR(TRACKING::FastKalmanFilterInnerTrack, ilayer, si_r_pos[ilayer],silicon_thickness,false);
      BarrelFastKalmanFilterConfigBARR(TRACKING::FastKalmanFilterSiliconTrack, ilayer, si_r_pos[ilayer],silicon_thickness,false);
    }
  }

  return;
}

#endif

