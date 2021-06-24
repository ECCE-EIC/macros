#ifndef MACRO_G4TTLEIC_C
#define MACRO_G4TTLEIC_C

#include "GlobalVariables.C"

#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4detectors/PHG4SectorSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>

#include <g4main/PHG4Reco.h>

#include <string>

R__LOAD_LIBRARY(libg4detectors.so)

int make_forward_station(string name, PHG4Reco *g4Reco, double zpos, double Rmin,
                         double Rmax, double tSilicon);
int make_barrel_layer(string name, PHG4Reco *g4Reco,
                      double radius, double halflength, double tSilicon);

//-----------------------------------------------------------------------------------//
namespace Enable
{
  bool FTTL = false;
  bool ETTL = false;
  bool CTTL = false;
}  // namespace Enable

namespace TTL
{
  // 2, LGAD based ToF by Wei Li:
  // Present the recent simulation studies and detector layout for the LGAD based ToF,
  // which can be used as an outer tracker. The current detector can reach around 30 micron spatial
  // resolution with 500 micron AC-LGAD sensor.
  // This detector can provide better coverage and provide further constraints in track
  // reconstruction within high pseudorapidity regions.

  const double PositionResolution(30e-4);
}  // namespace TTL

//-----------------------------------------------------------------------------------//
void TTL_Init()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 200.);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 350.);
}
//-----------------------------------------------------------------------------------//
void FTTLSetup(PHG4Reco *g4Reco)
{
  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = .1 * cm;
  const double um = 1e-3 * mm;

  make_forward_station("FTTL_0", g4Reco, 287, 3.5, 1.3, 85 * um);
  make_forward_station("FTTL_1", g4Reco, 289, 3.5, 1.3, 85 * um);
}

//-----------------------------------------------------------------------------------//
void ETTLSetup(PHG4Reco *g4Reco)
{
  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = .1 * cm;
  const double um = 1e-3 * mm;

  make_forward_station("ETTL_0", g4Reco, -190, -1.85, -3.5, 85 * um);  // define wit eta
  make_forward_station("ETTL_1", g4Reco, -192, -1.85, -3.5, 85 * um);  // define wit eta
}

//-----------------------------------------------------------------------------------//
void CTTLSetup(PHG4Reco *g4Reco)
{
  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = .1 * cm;
  const double um = 1e-3 * mm;

  make_barrel_layer("CTTL_0", g4Reco, 83, 180, 85 * um);
}

//-----------------------------------------------------------------------------------//
int make_forward_station(string name, PHG4Reco *g4Reco,
                         double zpos, double etamin, double etamax,
                         double tSilicon)  //silicon thickness
{
  //  cout
  //      << "make_GEM_station - GEM construction with PHG4SectorSubsystem - make_GEM_station_EdgeReadout  of "
  //      << name << endl;

  // always facing the interaction point
  double polar_angle = 0;
  if (zpos < 0)
  {
    zpos = -zpos;
    polar_angle = M_PI;
  }
  if (etamax < etamin)
  {
    double t = etamax;
    etamax = etamin;
    etamin = t;
  }

  PHG4SectorSubsystem *ttl;
  ttl = new PHG4SectorSubsystem(name);

  ttl->SuperDetector(name);

  ttl->get_geometry().set_normal_polar_angle(polar_angle);
  ttl->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm());
  ttl->get_geometry().set_min_polar_angle(PHG4Sector::Sector_Geometry::eta_to_polar_angle(etamax));
  ttl->get_geometry().set_max_polar_angle(PHG4Sector::Sector_Geometry::eta_to_polar_angle(etamin));
  ttl->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::ConeEdge());
  ttl->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::ConeEdge());
  ttl->get_geometry().set_N_Sector(1);
  ttl->get_geometry().set_material("G4_AIR");
  ttl->OverlapCheck(true);

  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = .1 * cm;
  const double um = 1e-3 * mm;
  // build up layers

  ttl->get_geometry().AddLayer("SiliconSensor", "G4_Si", tSilicon, true, 100);
  ttl->get_geometry().AddLayer("Metalconnection", "G4_Al", 100 * um, false, 100);
  ttl->get_geometry().AddLayer("HDI", "G4_KAPTON", 20 * um, false, 100);
  ttl->get_geometry().AddLayer("Cooling", "G4_WATER", 100 * um, false, 100);
  ttl->get_geometry().AddLayer("Support", "G4_GRAPHITE", 50 * um, false, 100);
  ttl->get_geometry().AddLayer("Support_Gap", "G4_AIR", 1 * cm, false, 100);
  ttl->get_geometry().AddLayer("Support2", "G4_GRAPHITE", 50 * um, false, 100);

  g4Reco->registerSubsystem(ttl);

  if (TRACKING::FastKalmanFilter)
  {
    TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             TTL::PositionResolution,           //      const float radres,
                                             TTL::PositionResolution,           //      const float phires,
                                             tSilicon / sqrt(12.),              //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
    TRACKING::FastKalmanFilter-> add_zplane_state(name, zpos);
    TRACKING::ProjectionNames.insert(name);
  }

  return 0;
}

//-----------------------------------------------------------------------------------//
int make_barrel_layer(string name, PHG4Reco *g4Reco,
                      double radius, double halflength, double tSilicon)
{
  //---------------------------------
  //build barrel layer
  //---------------------------------
  const int nSubLayer = 7;
  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = .1 * cm;
  const double um = 1e-3 * mm;

  string layerName[nSubLayer] = {"SiliconSensor", "Metalconnection", "HDI", "Cooling",
                                 "Support1", "Support_Gap", "Support2"};
  string material[nSubLayer] = {"G4_Si", "G4_Al", "G4_KAPTON", "G4_WATER",
                                "G4_GRAPHITE", "G4_AIR", "G4_GRAPHITE"};
  double thickness[nSubLayer] = {tSilicon, 15 * um, 20 * um, 100 * um,
                                 50 * um, 1, 50 * um};

  double max_bh_radius = 0.;
  PHG4CylinderSubsystem *cyl;
  //   cout << "started to create cylinder layer: " << name << endl;

  double currRadius = radius;
  //   cout << currRadius << endl;
  for (int l = 0; l < nSubLayer; l++)
  {
    //     cout << name <<"_"<< layerName[l] << endl;
    cyl = new PHG4CylinderSubsystem(name + "_" + layerName[l], l);
    cyl->SuperDetector(name);
    cyl->set_double_param("radius", currRadius);
    cyl->set_double_param("length", 2.0 * halflength);
    cyl->set_string_param("material", material[l]);
    cyl->set_double_param("thickness", thickness[l]);
    if (l == 0) cyl->SetActive();  //only the Silicon Sensor is active
    cyl->OverlapCheck(true);
    g4Reco->registerSubsystem(cyl);
    currRadius = currRadius + thickness[l];
    //     cout << currRadius << endl;
  }

  if (TRACKING::FastKalmanFilter)
  {
    TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + name,     //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
                                             tSilicon / sqrt(12.),        //      const float radres,
                                             TTL::PositionResolution,     //      const float phires,
                                             TTL::PositionResolution,     //      const float lonres,
                                             1,                           //      const float eff,
                                             0);                          //      const float noise
    TRACKING::FastKalmanFilter -> add_cylinder_state(name, radius);

    TRACKING::ProjectionNames.insert(name);
  }
  return 0;
}

#endif

//-----------------------------------------------------------------------------------//
