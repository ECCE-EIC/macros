/*!
 * \file G4_RICH.C
 * \brief Setup the gas RICH detector as designed in ePHENIX LOI
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.2 $
 * \date $Date: 2013/10/09 01:08:17 $
 */
#ifndef MACRO_G4dRICH_C
#define MACRO_G4dRICH_C

#include <GlobalVariables.C>

#include <eccefastpidreco/ECCEFastPIDReco.h>
#include <eccefastpidreco/ECCEdRICHFastPIDMap.h>
#include <g4drich/EICG4dRICHSubsystem.h>
#include <g4detectors/PHG4SectorSubsystem.h>
#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>

#include <g4main/PHG4Reco.h>

R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libEICG4dRICH.so)

R__LOAD_LIBRARY(libECCEFastPIDReco.so)

namespace Enable
{
  bool RICH = false;
  bool RICH_RECO = false;
  bool RICH_OVERLAPCHECK = false;
  int RICH_VERBOSITY = 0;
}  // namespace Enable

void RICHInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 210.);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 280.);
}

//! Fast geometry model for dRICH
//Grzegorz Kalicy <gkalicy@jlab.org>
//- it starts 180 cm from IP
//- radius of aerogel part starts at 110 cm at rises up to 120cm over 20 cm of length.
//- at 175 cm from IP it rapidly grows radially to radius 210cm at stays constant like a cylinder until end at 280cm from the IP
void RICHSetup(PHG4Reco *g4Reco)
{
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::RICH_OVERLAPCHECK;
  int verbosity = std::max(Enable::VERBOSITY, Enable::RICH_VERBOSITY);

  double z = 185;   //Start of dRICH
  double dz = 100;  //Length of dRICH

  EICG4dRICHSubsystem *drichSubsys = new EICG4dRICHSubsystem("dRICh");
  drichSubsys->SetGeometryFile(string(getenv("CALIBRATIONROOT")) + "/dRICH/mapping/drich-g4model_v5.txt");
  drichSubsys->set_double_param("place_z", z + dz * 0.5);  // relative position to mother vol.
  drichSubsys->OverlapCheck(OverlapCheck);
  drichSubsys->Verbosity(verbosity);
  drichSubsys->SetActive();

  g4Reco->registerSubsystem(drichSubsys);



  double Rmin = 10;
  double Rmax = 88;
  double zposOrg = 184.7 ;
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

  PHG4SectorSubsystem *dRICHDummy;
  dRICHDummy = new PHG4SectorSubsystem("dRICH_Plane_0");

  dRICHDummy->SuperDetector("dRICH_Plane_0");

  dRICHDummy->get_geometry().set_normal_polar_angle(polar_angle);
  dRICHDummy->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm());
  dRICHDummy->get_geometry().set_min_polar_angle(min_polar_angle);
  dRICHDummy->get_geometry().set_max_polar_angle(max_polar_angle);
  dRICHDummy->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::ConeEdge());
  dRICHDummy->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::ConeEdge());
  dRICHDummy->get_geometry().set_N_Sector(1);
  dRICHDummy->get_geometry().set_material("G4_AIR");
  dRICHDummy->OverlapCheck(true);  //true);//overlapcheck);
  dRICHDummy->get_geometry().AddLayer("dRICH_Plane_0_Air", "G4_AIR", 100* um, true, 100);
  dRICHDummy->SaveAllHits(true);

  g4Reco->registerSubsystem(dRICHDummy);


  Rmin = 13;
  Rmax = 170;
  double zposOrg2 = 185.3+100;
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

  PHG4SectorSubsystem *dRICHDummy2;
  dRICHDummy2 = new PHG4SectorSubsystem("dRICH_Plane_1");

  dRICHDummy2->SuperDetector("dRICH_Plane_1");

  dRICHDummy2->get_geometry().set_normal_polar_angle(polar_angle);
  dRICHDummy2->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm());
  dRICHDummy2->get_geometry().set_min_polar_angle(min_polar_angle);
  dRICHDummy2->get_geometry().set_max_polar_angle(max_polar_angle);
  dRICHDummy2->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::ConeEdge());
  dRICHDummy2->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::ConeEdge());
  dRICHDummy2->get_geometry().set_N_Sector(1);
  dRICHDummy2->get_geometry().set_material("G4_AIR");
  dRICHDummy2->OverlapCheck(true);  //true);//overlapcheck);
  dRICHDummy2->get_geometry().AddLayer("dRICH_Plane_1_Air", "G4_AIR", 100* um, true, 100);
  dRICHDummy2->SaveAllHits(true);

  g4Reco->registerSubsystem(dRICHDummy2);

  
   if (TRACKING::FastKalmanFilter)
  {
    TRACKING::FastKalmanFilter->add_phg4hits(string(Form("G4HIT_dRICH_Plane_0")),     //      const std::string& phg4hitsNames,
                                                  PHG4TrackFastSim::Vertical_Plane,
                                                  999.,                      // radial-resolution [cm]
                                                  999.,  // azimuthal-resolution [cm]
                                                  999.,  // z-resolution [cm]
                                                  0.0001,                         // efficiency,
                                                  0);
    // project to an reference plane at z=170 cm
    TRACKING::FastKalmanFilter->add_zplane_state("dRICH_Plane_0", zposOrg);
    TRACKING::ProjectionNames.insert("dRICH_Plane_0");

    TRACKING::FastKalmanFilter->add_phg4hits(string(Form("G4HIT_dRICH_Plane_1")),     //      const std::string& phg4hitsNames,
                                                  PHG4TrackFastSim::Vertical_Plane,
                                                  999.,                      // radial-resolution [cm]
                                                  999.,  // azimuthal-resolution [cm]
                                                  999.,  // z-resolution [cm]
                                                  0.0001,                         // efficiency,
                                                  0);
    // project to an reference plane at z=170 cm
    TRACKING::FastKalmanFilter->add_zplane_state("dRICH_Plane_1", zposOrg2);
    TRACKING::ProjectionNames.insert("dRICH_Plane_1");
  }
}

void RICHReco()
{
  const int verbosity = std::max(Enable::VERBOSITY, Enable::RICH_VERBOSITY);
  Fun4AllServer *se = Fun4AllServer::instance();

  ECCEdRICHFastPIDMap *pidmap = new ECCEdRICHFastPIDMap();
  pidmap->Verbosity(verbosity);
  pidmap->dualRICH_aerogel();

  ECCEFastPIDReco *reco = new ECCEFastPIDReco(pidmap, EICPIDDefs::dRICH_AeroGel, "ECCEFastPIDReco-dRICH_AeroGel");
  reco->Verbosity(verbosity);

  se->registerSubsystem(reco);

  pidmap = new ECCEdRICHFastPIDMap();
  pidmap->Verbosity(verbosity);
  pidmap->dualRICH_C2F6();

  reco = new ECCEFastPIDReco(pidmap, EICPIDDefs::dRICH_Gas, "ECCEFastPIDReco-dRICH_Gas");
  reco->Verbosity(verbosity);

  se->registerSubsystem(reco);
}

#endif
