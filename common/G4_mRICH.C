#pragma once

#include "GlobalVariables.C"

#include <g4main/PHG4Reco.h>
#include <g4detectors/PHG4SectorSubsystem.h>
#include <g4mrich/PHG4mRICHSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <eccefastpidreco/ECCEFastPIDReco.h>
#include <eccefastpidreco/ECCEmRICHFastPIDMap.h>
/*!
 * \file G4_mRICH.C
 * \brief Aerogel mRICH for EIC detector
 * \author Murad Sarsour
 * \date $Date: 2020/7/2  $
 */

R__LOAD_LIBRARY(libECCEFastPIDReco.so)

namespace Enable
{
  bool mRICH = false;
  bool mRICH_RECO = false;
  bool mRICH_OVERLAPCHECK = false;
  int mRICH_VERBOSITY = 0;
}  // namespace Enable

void mRICHInit()
{
}

//-1: single module
// 0: build h-side sectors and e-side wall
// 1: build h-side sectors
// 2: build e-side wall
// 3: build h-side wall
// 4: build h-side wall and e-side wall
// 5: build b-side sectors
// 6: build projected e-side wall

void mRICHSetup(PHG4Reco *g4Reco, const int detectorSetup = 1,  //1: full setup; 0:skeleton
                const int mRICHsystemSetup = 6                  //5//2//-1//
)
{
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::mRICH_OVERLAPCHECK;

  PHG4mRICHSubsystem *mRICH = new PHG4mRICHSubsystem("mRICH", 1);
  mRICH->set_int_param("detectorSetup", detectorSetup);
  mRICH->set_int_param("subsystemSetup", mRICHsystemSetup);
  mRICH->UseCalibFiles(PHG4DetectorSubsystem::xml);
  mRICH->SetCalibrationFileDir(string(getenv("CALIBRATIONROOT")) + string("/mRICH/GeometryV2/"));
  mRICH->OverlapCheck(OverlapCheck);
  //  mRICH->Verbosity(5);

  g4Reco->registerSubsystem(mRICH);


  double Rmin = 6;
  double Rmax = 57;
  double zposOrg = -133.5 ;
  double zpos = zposOrg;
  
  double min_polar_angle = atan2(Rmin, zpos);
  double max_polar_angle = atan2(Rmax, zpos);

  // always facing the interaction point
  double polar_angle = 0;
  if (zpos < 0)
  {
    zpos = -zpos;
    polar_angle = M_PI;
  }
  if (max_polar_angle < min_polar_angle)
  {
    double t = max_polar_angle;
    max_polar_angle = min_polar_angle;
    min_polar_angle = t;
  }
  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = .1 * cm;
  const double um = 1e-3 * mm;

  PHG4SectorSubsystem *mRichDummy;
  mRichDummy = new PHG4SectorSubsystem("mRICH_Plane_0");

  mRichDummy->SuperDetector("mRICH_Plane_0");

  mRichDummy->get_geometry().set_normal_polar_angle(polar_angle);
  mRichDummy->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm());
  mRichDummy->get_geometry().set_min_polar_angle(min_polar_angle);
  mRichDummy->get_geometry().set_max_polar_angle(max_polar_angle);
  mRichDummy->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::ConeEdge());
  mRichDummy->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::ConeEdge());
  mRichDummy->get_geometry().set_N_Sector(1);
  mRichDummy->get_geometry().set_material("G4_AIR");
  mRichDummy->OverlapCheck(true);  //true);//overlapcheck);
  mRichDummy->get_geometry().AddLayer("mRICH_Plane_0_Air", "G4_AIR", 100* um, true, 100);

  g4Reco->registerSubsystem(mRichDummy);

  Rmin = 7;
  Rmax = 57;
  double zposOrg2 = -162. ;
  zpos = zposOrg2;
  
  min_polar_angle = atan2(Rmin, zpos);
  max_polar_angle = atan2(Rmax, zpos);

  // always facing the interaction point
  polar_angle = 0;
  if (zpos < 0)
  {
    zpos = -zpos;
    polar_angle = M_PI;
  }
  if (max_polar_angle < min_polar_angle)
  {
    double t = max_polar_angle;
    max_polar_angle = min_polar_angle;
    min_polar_angle = t;
  }

  PHG4SectorSubsystem *mRichDummy2;
  mRichDummy2 = new PHG4SectorSubsystem("mRICH_Plane_1");

  mRichDummy2->SuperDetector("mRICH_Plane_1");

  mRichDummy2->get_geometry().set_normal_polar_angle(polar_angle);
  mRichDummy2->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm());
  mRichDummy2->get_geometry().set_min_polar_angle(min_polar_angle);
  mRichDummy2->get_geometry().set_max_polar_angle(max_polar_angle);
  mRichDummy2->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::ConeEdge());
  mRichDummy2->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::ConeEdge());
  mRichDummy2->get_geometry().set_N_Sector(1);
  mRichDummy2->get_geometry().set_material("G4_AIR");
  mRichDummy2->OverlapCheck(true);  //true);//overlapcheck);
  mRichDummy2->get_geometry().AddLayer("mRICH_Plane_1_Air", "G4_AIR", 100* um, true, 100);

  g4Reco->registerSubsystem(mRichDummy2);

  if (TRACKING::FastKalmanFilter)
  {
    TRACKING::FastKalmanFilter->add_phg4hits(string(Form("G4HIT_mRICH_Plane_0")),     //      const std::string& phg4hitsNames,
                                                  PHG4TrackFastSim::Vertical_Plane,
                                                  999.,                      // radial-resolution [cm]
                                                  999.,  // azimuthal-resolution [cm]
                                                  999.,  // z-resolution [cm]
                                                  0.0001,                         // efficiency,
                                                  0);                          // noise hits

    TRACKING::FastKalmanFilter-> add_zplane_state("mRICH_Plane_0", zposOrg);
    TRACKING::ProjectionNames.insert("mRICH_Plane_0");

    TRACKING::FastKalmanFilter->add_phg4hits(string(Form("G4HIT_mRICH_Plane_1")),     //      const std::string& phg4hitsNames,
                                                  PHG4TrackFastSim::Vertical_Plane,
                                                  999.,                      // radial-resolution [cm]
                                                  999.,  // azimuthal-resolution [cm]
                                                  999.,  // z-resolution [cm]
                                                  0.0001,                         // efficiency,
                                                  0);                          // noise hits

    TRACKING::FastKalmanFilter-> add_zplane_state("mRICH_Plane_1", zposOrg2);
    TRACKING::ProjectionNames.insert("mRICH_Plane_1");
  }
}



void mRICHReco(int verbosity = 0)
{
  verbosity = std::max(Enable::VERBOSITY, verbosity);
  verbosity = std::max(Enable::mRICH_VERBOSITY, verbosity);

  Fun4AllServer *se = Fun4AllServer::instance();

  ECCEmRICHFastPIDMap * pidmap = new ECCEmRICHFastPIDMap();
  pidmap->Verbosity(verbosity);

  ECCEFastPIDReco * reco = new ECCEFastPIDReco(pidmap, EICPIDDefs::mRICH, "ECCEFastPIDReco-mRICH");
  reco->Verbosity(verbosity);

  se->registerSubsystem(reco);
}

void mRICH_Eval(std::string outputfile, int verbosity = 0)
{
  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4eval.so");
  Fun4AllServer *se = Fun4AllServer::instance();

  return;
}
