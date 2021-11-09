#ifndef MACRO_G4GEMEIC_C
#define MACRO_G4GEMEIC_C

#include <GlobalVariables.C>

#include <g4trackfastsim/PHG4TrackFastSim.h>

#include <g4detectors/PHG4SectorSubsystem.h>

#include <g4main/PHG4Reco.h>

#include <string>

R__LOAD_LIBRARY(libg4detectors.so)

int make_GEM_station(string name, PHG4Reco *g4Reco, double zpos, double etamin, double etamax, const int N_Sector = 16, double tilt = 0, bool doTilt = false);
void AddLayers_MiniTPCDrift(PHG4SectorSubsystem *gem);
void AddLayers_GEMDrift(PHG4SectorSubsystem *gem);

namespace Enable
{
  bool EGEM = false;
  bool FGEM = false;
  bool BGEM = false;
}  // namespace Enable

void EGEM_Init()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 80.);
  // extends only to -z
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -160.);
}

void FGEM_Init()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 150.);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 282.);
}

void BGEM_Init()
{
  cout << __PRETTY_FUNCTION__ << ": BGEM not yet defined" << endl;
  exit(1);

  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 150.);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 150.);
}

void EGEMSetup(PHG4Reco *g4Reco)
{
  make_GEM_station("EGEM_0", g4Reco, -121.0, -1.668, -3.7);
}

void FGEMSetup(PHG4Reco *g4Reco, const int N_Sector = 16)
{
  make_GEM_station("FGEM_0", g4Reco, 287.0, 1.3, 3.6, N_Sector);
}

//! Add drift layers to mini TPC
void AddLayers_MiniTPCDrift(PHG4SectorSubsystem *gem)
{
  assert(gem);

  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = 0.1 * cm;
  const double um = 1e-3 * mm;

  //  const int N_Layers = 70; // used for mini-drift TPC timing digitalization
  const int N_Layers = 1;  // simplified setup
  const double thickness = 2 * cm;

  gem->get_geometry().AddLayer("EntranceWindow", "G4_MYLAR", 25 * um, false, 100);
  gem->get_geometry().AddLayer("Cathode", "G4_GRAPHITE", 10 * um, false, 100);

  for (int d = 1; d <= N_Layers; d++)
  {
    ostringstream s;
    s << "DriftLayer_";
    s << d;

    gem->get_geometry().AddLayer(s.str(), "G4_METHANE", thickness / N_Layers, true);
  }
}

//! Add drift layers to mini TPC
void AddLayers_GEMDrift(PHG4SectorSubsystem *gem)
{
  assert(gem);

  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = 0.1 * cm;
  const double um = 1e-3 * mm;

  const int N_Layers = 1;  // simplified setup
  const double thickness = 3 * mm;

  gem->get_geometry().AddLayer("EntranceWindow", "G4_MYLAR", 25 * um, false, 100);
  gem->get_geometry().AddLayer("Cathode", "G4_GRAPHITE", 10 * um, false, 100);

  for (int d = 1; d <= N_Layers; d++)
  {
    ostringstream s;
    s << "DriftLayer_";
    s << d;

    gem->get_geometry().AddLayer(s.str(), "G4_METHANE", thickness / N_Layers, true);
  }
}

int make_GEM_station(string name, PHG4Reco *g4Reco, double zpos, double etamin,
                     double etamax, const int N_Sector = 16, double tilt = 0, bool doTilt = false)
{
  //  cout
  //      << "make_GEM_station - GEM construction with PHG4SectorSubsystem - make_GEM_station_EdgeReadout  of "
  //      << name << endl;

  double zpos_lab(zpos);
  double polar_angle = 0;

  if (doTilt)
  {
    zpos = zpos - (zpos * sin(tilt) + zpos * cos(tilt) * tan(PHG4Sector::Sector_Geometry::eta_to_polar_angle(2) - tilt)) * sin(tilt);
  }
  else
  {
    if (zpos < 0)
    {
      zpos = -zpos;
      polar_angle = M_PI;
    }
  }
  if (etamax < etamin)
  {
    double t = etamax;
    etamax = etamin;
    etamin = t;
  }

  PHG4SectorSubsystem *gem;
  gem = new PHG4SectorSubsystem(name);

  gem->SuperDetector(name);

  if (doTilt)
  {
    gem->get_geometry().set_normal_polar_angle((PHG4Sector::Sector_Geometry::eta_to_polar_angle(etamin) + PHG4Sector::Sector_Geometry::eta_to_polar_angle(etamax)) / 2);
    gem->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm(), PHG4Sector::Sector_Geometry::eta_to_polar_angle(etamax));
  }
  else
  {
    gem->get_geometry().set_normal_polar_angle(polar_angle);
    gem->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm());
  }
  gem->get_geometry().set_min_polar_angle(PHG4Sector::Sector_Geometry::eta_to_polar_angle(etamax));
  gem->get_geometry().set_max_polar_angle(PHG4Sector::Sector_Geometry::eta_to_polar_angle(etamin));
  gem->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
  gem->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
  gem->get_geometry().set_N_Sector(N_Sector);
  gem->get_geometry().set_material("G4_METHANE");
  gem->OverlapCheck(Enable::OVERLAPCHECK);

  AddLayers_GEMDrift(gem);
  gem->get_geometry().AddLayers_HBD_GEM();
  gem->OverlapCheck(Enable::OVERLAPCHECK);
  g4Reco->registerSubsystem(gem);

  // Following Nov-5 tracking meeting, update to muRwell performance of 55um resolution 2D readout
  if (TRACKING::FastKalmanFilter)
  {
    TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             55e-4,                    //      const float radres,
                                             55e-4,                             //      const float phires,
                                             100e-4,                            //      const float lonres,
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise

    TRACKING::FastKalmanFilter->add_zplane_state(name, zpos_lab);
    TRACKING::ProjectionNames.insert(name);
  }

  if (TRACKING::FastKalmanFilterInnerTrack and zpos_lab<0)
    TRACKING::FastKalmanFilterInnerTrack->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                                       PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                                       55e-4,                    //      const float radres,
                                                       55e-4,                             //      const float phires,
                                                       100e-4,                            //      const float lonres,
                                                       1,                                 //      const float eff,
                                                       0);                                //      const float noise

  return 0;
}
#endif

