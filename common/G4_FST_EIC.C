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

#include <g4detectors/PHG4SectorSubsystem.h>
#include <g4detectors/PHG4CylinderSubsystem.h>

#include <g4main/PHG4Reco.h>

#include <string>

R__LOAD_LIBRARY(libg4detectors.so)

int make_LANL_FST_station(string name, PHG4Reco *g4Reco, double zpos, double Rmin,
                          double Rmax,double tSilicon);
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
    bool SUPPORTCYL = false;
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

  if (!Enable::FST_OVERLAPCHECK) 
  {
    Enable::FST_OVERLAPCHECK=Enable::OVERLAPCHECK;
  }
  //bool OverlapCheck = Enable::OVERLAPCHECK || Enable::FST_OVERLAPCHECK;

  //Design from Xuan Li @LANL
  if (G4FST::SETTING::FST_TPC) // this ver. fits in the TPC
    {
      make_LANL_FST_station("FST_0", g4Reco, 35, 4, 17, 35 * um);  //cm
      make_LANL_FST_station("FST_1", g4Reco, 53, 4.5, 17, 35 * um);
      make_LANL_FST_station("FST_2", g4Reco, 77, 5, 17, 35 * um);
      make_LANL_FST_station("FST_3", g4Reco, 101, 7.5, 17, 85 * um);
      make_LANL_FST_station("FST_4", g4Reco, 125, 9.5, 45, 85 * um);
    }
  else
    {
      make_LANL_FST_station("FST_0", g4Reco, 35, 4, 22, 35 * um);  //cm
      make_LANL_FST_station("FST_1", g4Reco, 57.5, 4.5, 42, 35 * um);
      make_LANL_FST_station("FST_2", g4Reco, 80, 6, 43.5, 35 * um);
      make_LANL_FST_station("FST_3", g4Reco, 115, 9.3, 46.8, 85 * um);
      make_LANL_FST_station("FST_4", g4Reco, 125, 9.6, 47.1, 85 * um);

      if (G4FST::SETTING::SUPPORTCYL)
	{
	  double gap=8;  //cm
	  double tSupport=0.5;  //cm
	  make_supportCyl("FSTSupportCyl",g4Reco, 50.1+gap, tSupport, 125.0*2.0);  //cm
	}
    }
}
//-----------------------------------------------------------------------------------//
int make_LANL_FST_station(string name, PHG4Reco *g4Reco,
                          double zpos, double Rmin, double Rmax,double tSilicon) //silicon thickness
{
  // always facing the interaction point
  double polar_angle = 0;
  if (zpos < 0)
    {
      zpos = -zpos;
      polar_angle = M_PI;
    }

  double min_polar_angle = atan2(Rmin, zpos);
  double max_polar_angle = atan2(Rmax, zpos);

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
  fst->OverlapCheck(Enable::FST_OVERLAPCHECK);//true);//overlapcheck);

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
 return 0;
}
//-----------------------------------------------------------------------------------//
int make_supportCyl(string name, PHG4Reco *g4Reco,double r, double t, double length)
{

  PHG4CylinderSubsystem* cyl = new PHG4CylinderSubsystem(name, 5);
  cyl->set_double_param("radius", r);
  cyl->set_double_param("length", length);
  cyl->set_string_param("material", "G4_GRAPHITE");
  cyl->set_double_param("thickness", t);
  cyl->set_double_param("place_x", 0.);
  cyl->set_double_param("place_y", 0.);
  cyl->set_double_param("place_z", 0);
  cyl->SetActive(0);
  //cyl->SuperDetector("");
  cyl->OverlapCheck(Enable::FST_OVERLAPCHECK);//OverlapCheck);

  g4Reco->registerSubsystem(cyl);
  return 0;
}
#endif
