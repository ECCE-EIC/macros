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
}  // namespace Enable
//-----------------------------------------------------------------------------------//
void BarrelInit()
{
}


void BarrelFastKalmanFilterConfigSVTX(PHG4TrackFastSim * kalman_filter, int ilay, double radius, bool addproj)
{

  // import Kalman filter config (lines 226 to 246 here: https://github.com/eic/g4lblvtx/blob/master/macros/auxiliary_studies/simplified_geometry/Fun4All_G4_simplified_v2.C):

  // add Vertexing Layers
  kalman_filter->add_phg4hits(
      Form("G4HIT_SVTX_%d",ilay),  // const std::string& phg4hitsNames,
      PHG4TrackFastSim::Cylinder,
      999.,                      // radial-resolution [cm]
      10. / 10000. / sqrt(12.),  // azimuthal-resolution [cm]
      10. / 10000. / sqrt(12.),  // z-resolution [cm]
      1,                         // efficiency,
      0                          // noise hits
  );
  kalman_filter->add_cylinder_state(Form("SVTX_%d",ilay), radius);
  if(addproj)TRACKING::ProjectionNames.insert(Form("SVTX_%d",ilay));
}

void BarrelFastKalmanFilterConfigBARR(PHG4TrackFastSim * kalman_filter, int ilay, double radius, bool addproj)
{
  // add Barrel Layers
  kalman_filter->add_phg4hits(
      Form("G4HIT_BARR_%d",ilay),  // const std::string& phg4hitsNames,
      PHG4TrackFastSim::Cylinder,
      999.,                      // radial-resolution [cm]
      10. / 10000. / sqrt(12.),  // azimuthal-resolution [cm]
      10. / 10000. / sqrt(12.),  // z-resolution [cm]
      1,                         // efficiency,
      0                          // noise hits
  );
  kalman_filter->add_cylinder_state(Form("BARR_%d",ilay), radius);
  if(addproj)TRACKING::ProjectionNames.insert(Form("BARR_%d",ilay));

}


//-----------------------------------------------------------------------------------//
void Barrel(PHG4Reco *g4Reco, int det_ver = 3)
{
  // import Geometry (lines 111 to 148 in https://github.com/eic/g4lblvtx/blob/master/macros/auxiliary_studies/simplified_geometry/Fun4All_G4_simplified_v2.C):
  PHG4CylinderSubsystem *cyl(nullptr);

  //---------------------------
  // Vertexing
  double si_vtx_r_pos[] = {3.30, 5.70, 8.10};
  const int nVtxLayers = sizeof(si_vtx_r_pos) / sizeof(*si_vtx_r_pos);
  for (int ilayer = 0; ilayer < nVtxLayers; ilayer++)
  {
    cyl = new PHG4CylinderSubsystem(Form("SVTX_%d",ilayer), ilayer);
    cyl->set_string_param("material", "G4_Si");
    cyl->set_double_param("radius", si_vtx_r_pos[ilayer]);
    cyl->set_double_param("thickness", 0.05 / 100. * 9.37);
    cyl->set_double_param("place_z", 0);
    cyl->set_double_param("length", 30.);
    cyl->SetActive();
    cyl->SuperDetector(Form("SVTX_%d",ilayer));
    g4Reco->registerSubsystem(cyl);

    BarrelFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilter, ilayer, si_vtx_r_pos[ilayer],true);
    BarrelFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilterInnerTrack, ilayer, si_vtx_r_pos[ilayer],false);
    BarrelFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilterSiliconTrack, ilayer, si_vtx_r_pos[ilayer],false);
  }
  //---------------------------
  // Barrel
  double si_r_pos[] = {21., 22.68};
  const int nTrckLayers = sizeof(si_r_pos) / sizeof(*si_r_pos);
  double si_z_length[] = {60., 60.};
  for (int ilayer = 0; ilayer < nTrckLayers; ilayer++)
  {
    cyl = new PHG4CylinderSubsystem(Form("BARR_%d",ilayer), ilayer);
    cyl->set_string_param("material", "G4_Si");
    cyl->set_double_param("radius", si_r_pos[ilayer]);
    cyl->set_double_param("thickness", 0.05 / 100. * 9.37);
    cyl->set_double_param("place_z", 0);
    cyl->set_double_param("length", si_z_length[ilayer]);
    cyl->SetActive();
    cyl->SuperDetector(Form("BARR_%d",ilayer));
    g4Reco->registerSubsystem(cyl);

    BarrelFastKalmanFilterConfigBARR(TRACKING::FastKalmanFilter, ilayer, si_r_pos[ilayer],true);
    BarrelFastKalmanFilterConfigBARR(TRACKING::FastKalmanFilterInnerTrack, ilayer, si_r_pos[ilayer],false);
    BarrelFastKalmanFilterConfigBARR(TRACKING::FastKalmanFilterSiliconTrack, ilayer, si_r_pos[ilayer],false);
  }

}

#endif
