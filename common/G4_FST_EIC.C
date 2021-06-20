/*---------------------------------------------------------------------*
 * Barrel tracker designed by LANL EIC team                            *
 * See technical notes for details: arXiv:2009.02888                   *
 * Contact Ping and Xuan @LANL for questions:                          *
 *   Xuan: xuanli@lanl.gov                                             *
 *   Ping: cpwong@lanl.gov                                             *
 *---------------------------------------------------------------------*/

#ifndef MACRO_G4FSTEIC_C
#define MACRO_G4FSTEIC_C

#include "GlobalVariables.C"

#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4detectors/PHG4SectorSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>

#include <g4main/PHG4Reco.h>

#include <string>

R__LOAD_LIBRARY(libg4detectors.so)

int make_LANL_FST_station(string name, PHG4Reco *g4Reco, double zpos, double Rmin,
                          double Rmax, double tSilicon, double pitch);
int make_supportCyl(string name, PHG4Reco *g4Reco,
                    double r, double t, double length);
//-----------------------------------------------------------------------------------//
namespace Enable
{
  static bool FST = false;
  bool FST_OVERLAPCHECK = false;
}  // namespace Enable

namespace G4FST
{
  namespace SETTING
  {
    bool FST_TPC = false;
    bool SUPPORTCYL = true;
  }  // namespace SETTING
}  // namespace G4FST

//-----------------------------------------------------------------------------------//
void FST_Init()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 48.);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 127.);
  if (G4FST::SETTING::SUPPORTCYL)
  {
    BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -127.);
  }
}
//-----------------------------------------------------------------------------------//
void FSTSetup(PHG4Reco *g4Reco, const double min_eta = 1.245)
{
  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = .1 * cm;
  const double um = 1e-3 * mm;

  //Design from Xuan Li @LANL
  make_LANL_FST_station("FST_0", g4Reco, 35, 4, 22, 35 * um, 20e-4);  //cm
  make_LANL_FST_station("FST_1", g4Reco, 57.5, 4.5, 42, 35 * um, 20e-4);
  make_LANL_FST_station("FST_2", g4Reco, 80, 6, 43.5, 35 * um, 20e-4);
  make_LANL_FST_station("FST_3", g4Reco, 115, 9.3, 46.8, 85 * um, 36.4e-4);
  make_LANL_FST_station("FST_4", g4Reco, 125, 9.6, 47.1, 85 * um, 36.4e-4);

  //mirror for e-going FST
  make_LANL_FST_station("EFST_0", g4Reco, -35, 4, 22, 35 * um, 20e-4);  //cm
  make_LANL_FST_station("EFST_1", g4Reco, -57.5, 4.5, 42, 35 * um, 20e-4);
  make_LANL_FST_station("EFST_2", g4Reco, -80, 6, 43.5, 35 * um, 20e-4);
  make_LANL_FST_station("EFST_3", g4Reco, -115, 9.3, 46.8, 85 * um, 36.4e-4);
  make_LANL_FST_station("EFST_4", g4Reco, -125, 9.6, 47.1, 85 * um, 36.4e-4);

  if (G4FST::SETTING::SUPPORTCYL)
  {
    double gap = 8;                                                               //cm
    double tSupport = 0.2;                                                        //cm
    make_supportCyl("FSTSupportCyl", g4Reco, 50.1 + gap, tSupport, 125.0 * 2.0);  //cm
  }
}
//-----------------------------------------------------------------------------------//
int make_LANL_FST_station(string name, PHG4Reco *g4Reco,
                          double zpos, double Rmin, double Rmax, double tSilicon, double pitch)  //silicon thickness
{
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::FST_OVERLAPCHECK;

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
  PHG4SectorSubsystem *fst;
  fst = new PHG4SectorSubsystem(name);

  fst->SuperDetector(name);

  fst->get_geometry().set_normal_polar_angle(polar_angle);
  fst->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm());
  fst->get_geometry().set_min_polar_angle(min_polar_angle);
  fst->get_geometry().set_max_polar_angle(max_polar_angle);
  fst->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::ConeEdge());
  fst->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::ConeEdge());
  fst->get_geometry().set_N_Sector(1);
  fst->get_geometry().set_material("G4_AIR");
  fst->OverlapCheck(OverlapCheck);  //true);//overlapcheck);

  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = .1 * cm;
  const double um = 1e-3 * mm;
  // build up layers

  fst->get_geometry().AddLayer("SiliconSensor", "G4_Si", tSilicon, true, 100);
  fst->get_geometry().AddLayer("Metalconnection", "G4_Al", 15 * um, false, 100);
  fst->get_geometry().AddLayer("HDI", "G4_KAPTON", 20 * um, false, 100);
  fst->get_geometry().AddLayer("Cooling", "G4_WATER", 100 * um, false, 100);
  fst->get_geometry().AddLayer("Support", "G4_GRAPHITE", 50 * um, false, 100);
  fst->get_geometry().AddLayer("Support_Gap", "G4_AIR", 1 * cm, false, 100);
  fst->get_geometry().AddLayer("Support2", "G4_GRAPHITE", 50 * um, false, 100);

  g4Reco->registerSubsystem(fst);

  if (TRACKING::FastKalmanFilter)
  {
    TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             pitch / sqrt(12.),                 //      const float radres,
                                             pitch / sqrt(12.),                 //      const float phires,
                                             50e-4 / sqrt(12.),                 //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
  }
  return 0;
}
//-----------------------------------------------------------------------------------//
int make_supportCyl(string name, PHG4Reco *g4Reco, double r, double t, double length)
{
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::FST_OVERLAPCHECK;

  PHG4CylinderSubsystem *cyl = new PHG4CylinderSubsystem(name, 5);
  cyl->set_double_param("radius", r);
  cyl->set_double_param("length", length);
  cyl->set_string_param("material", "CFRP_INTT");  // borrow carbon fiber reinforced polymer used in sPHENIX silicon tracker support
  cyl->set_double_param("thickness", t);
  cyl->set_double_param("place_x", 0.);
  cyl->set_double_param("place_y", 0.);
  cyl->set_double_param("place_z", 0);
  cyl->OverlapCheck(OverlapCheck);  //true);//overlapcheck);
  cyl->SetActive(0);
  //cyl->SuperDetector("");
  cyl->OverlapCheck(Enable::FST_OVERLAPCHECK);  //OverlapCheck);

  g4Reco->registerSubsystem(cyl);
  return 0;
}
#endif
