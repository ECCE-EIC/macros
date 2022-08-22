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
#include <fst/PHG4FSTSubsystem.h>

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
  bool FST_ABSORBER = false;
  bool FST_CELL = false;
  bool FST_TOWER = false;
  bool FST_CLUSTER = false;
  bool FST_EVAL = false;
  int FST_VERBOSITY = 0;
}  // namespace Enable

namespace G4FST
{
  namespace SETTING
  {
    bool FST_TPC = false;
    bool SUPPORTCYL = false;
    bool EPIC_TRACKINGGEO_VARIANT = false;
    double TRACKING_EFFICIENCY = 1.00;
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

/* 

The idea is to calculate all the rmax rmin based on z and theta1 and theta2 in a python script (obj_fun.py)
Once calculated just pass the disk dimensions into this script

*/

void FSTSetup(PHG4Reco *g4Reco)
{
  // if(Enable::EPIC_TRACKINGGEO){

  cout << "FST/EST: Using tracking efficiency of: " << G4FST::SETTING::TRACKING_EFFICIENCY << endl;


  const bool AbsorberActive = Enable::ABSORBER || Enable::FST_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::FST_OVERLAPCHECK;
  Fun4AllServer *se = Fun4AllServer::instance();

  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = .1 * cm;
  const double um = 1e-3 * mm;
  //const double max_radius = 50.;


  // AI optimized values
  const double bkwd_z[] = {33.2, 58.29, 80.05, 107.4};
  double bkwd_rmin[] = {3.6, 3.6, 5.25, 6.4};
  double bkwd_rmax[] = {15.3, 27.3, 35.25, 48.4};

  // EPIC setup based on Ernst
  const double bkwd_z_EPIC[] = {25.0, 45.0, 70.0, 100.0, 135.0};
  double bkwd_rmin_EPIC[] = {3.6, 3.6, 3.6, 3.9, 4.5};
  double bkwd_offset_EPIC[] = {0.0, 0.0, 0.0, 0.2, 0.7};
  double bkwd_rmax_EPIC[] = {19.0, 43.0, 43.0, 43.0, 59.0};

  // AI optimized z positions in EPIC service cone
  const double bkwd_z_EPIC2[] = {33.2, 58.29, 80.05, 107.4};
  double bkwd_rmin_EPIC2[] = {3.6, 3.6, 5.25, 6.4};
  double bkwd_offset_EPIC2[] = {0.0, 0.0, 0.0, 0.0};
  double bkwd_rmax_EPIC2[] = {15.3, 27.3, 35.25, 46.4};
  // const double bkwd_z_EPIC2[] = {33.2, 58.29, 80.05, 107.4};
  // double bkwd_rmin_EPIC2[] = {3.6, 3.6, 3.6, 3.9};
  // double bkwd_offset_EPIC2[] = {0.0, 0.0, 0.0, 0.3};
  // double bkwd_rmax_EPIC2[] = {28.3, 43.0, 43.0, 45.0};


  // non-projective values
  const double bkwd_z_np[] = {25, 52, 79, 106};
  double bkwd_rmin_np[] = {3.6, 3.6, 4.5, 5.5};
  double bkwd_rmax_np[] = {18.5, 36.5, 40.5, 41.5};

  int n_bkwd_disk = sizeof(bkwd_z) / sizeof(*bkwd_z);
  if( Enable::EPIC_TRACKINGGEO)
  {
    if(G4FST::SETTING::EPIC_TRACKINGGEO_VARIANT)
    {
      n_bkwd_disk = sizeof(bkwd_z) / sizeof(*bkwd_z);
    } else {
      n_bkwd_disk = sizeof(bkwd_z_EPIC) / sizeof(*bkwd_z_EPIC);
    }
  }

  for (unsigned int i = 0; i < n_bkwd_disk; i++)
  {
    if(Enable::EPIC_TRACKINGGEO){
      string name = Form("EST_%i", i);
      double pitch = 10e-4;
      PHG4FSTSubsystem *hfst = new PHG4FSTSubsystem(name);
      hfst->OverlapCheck(OverlapCheck);
      hfst->SuperDetector(name);
      hfst->SetActive();
      if (AbsorberActive) hfst->SetAbsorberActive();
      hfst->set_double_param("z_position", G4FST::SETTING::EPIC_TRACKINGGEO_VARIANT ? -1*bkwd_z_EPIC2[i] : -1*bkwd_z_EPIC[i]);
      hfst->set_double_param("r_min", G4FST::SETTING::EPIC_TRACKINGGEO_VARIANT ? bkwd_rmin_EPIC2[i] : bkwd_rmin_EPIC[i]);
      hfst->set_double_param("r_max", G4FST::SETTING::EPIC_TRACKINGGEO_VARIANT ? bkwd_rmax_EPIC2[i] : bkwd_rmax_EPIC[i]);
      hfst->set_double_param("silicon_thickness", 35 * um);
      hfst->set_double_param("offset_cutout", G4FST::SETTING::EPIC_TRACKINGGEO_VARIANT ? bkwd_offset_EPIC2[i] : bkwd_offset_EPIC[i]);
      g4Reco->registerSubsystem(hfst);

      if (TRACKING::FastKalmanFilter)
      {
        TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                                PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                                pitch / sqrt(12.),                 //      const float radres,
                                                pitch / sqrt(12.),                 //      const float phires,
                                                50e-4 / sqrt(12.),                 //      const float lonres, *ignored in plane detector*
                                                G4FST::SETTING::TRACKING_EFFICIENCY,                                 //      const float eff,
                                                0);                                //      const float noise
        TRACKING::FastKalmanFilterInnerTrack->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                                PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                                pitch / sqrt(12.),                 //      const float radres,
                                                pitch / sqrt(12.),                 //      const float phires,
                                                50e-4 / sqrt(12.),                 //      const float lonres, *ignored in plane detector*
                                                G4FST::SETTING::TRACKING_EFFICIENCY,                                 //      const float eff,
                                                0);                                //      const float noise
        TRACKING::FastKalmanFilterSiliconTrack->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                                PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                                pitch / sqrt(12.),                 //      const float radres,
                                                pitch / sqrt(12.),                 //      const float phires,
                                                50e-4 / sqrt(12.),                 //      const float lonres, *ignored in plane detector*
                                                G4FST::SETTING::TRACKING_EFFICIENCY,                                 //      const float eff,
                                                0);                                //      const float noise
      }
    } else if(Enable::AI_TRACKINGGEO){
      make_LANL_FST_station(Form("EST_%i", i), g4Reco, -1*bkwd_z[i], bkwd_rmin[i], bkwd_rmax[i], 35 * um, 10e-4);  //cm
    } else {
      make_LANL_FST_station(Form("EST_%i", i), g4Reco, -1*bkwd_z_np[i], bkwd_rmin_np[i], bkwd_rmax_np[i], 35 * um, 10e-4);  //cm
    }
  }

  // AI optimized values
  const double fwd_z[] = {33.2, 58.29, 79.85, 115, 125};
  double fwd_rmin[] = {3.6, 3.6, 5.2, 6.5, 7.5};
  double fwd_rmax[] = {15.3, 27.3, 35.2, 49.5, 49.5};

  // Enable::EPIC_TRACKINGGEO
  const double fwd_z_EPIC[] = {25.0, 45.0, 70.0, 100.0, 135.0};
  double fwd_rmin_EPIC[] = {3.6, 3.6, 3.6, 4.5, 5.4};
  double fwd_offset_EPIC[] = {0.0, 0.0, 0.0, -0.8, -1.7};
  double fwd_rmax_EPIC[] = {19.0, 43.0, 43.0, 43.0, 53.0};

  // AI optimized z positions in EPIC service cone
  // const double fwd_z_EPIC2[] = {33.2, 58.29, 79.85, 115, 125};
  // double fwd_rmin_EPIC2[] = {3.6, 3.6, 3.9, 4.85, 5.1};
  // double fwd_offset_EPIC2[] = {0.0, 0.0, -0.2, -1.2, -1.4};
  // double fwd_rmax_EPIC2[] = {15.3, 27.3, 35.2, 49.5, 49.5};
  const double fwd_z_EPIC2[] = {33.2, 58.29, 79.85, 115, 125};
  double fwd_rmin_EPIC2[] = {3.6, 3.6, 5.2, 6.5, 7.5};
  double fwd_offset_EPIC2[] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double fwd_rmax_EPIC2[] = {15.3, 27.3, 35.2, 44.5, 49.5};

  // non-projective values
  const double fwd_z_np[] = {25, 52, 73, 106, 125};
  double fwd_rmin_np[] = {3.6, 3.6, 4.5, 5.5, 7.5};
  double fwd_rmax_np[] = {18.5, 36.5, 40.5, 41.5, 43.4};

  const int n_fwd_disk = sizeof(fwd_z) / sizeof(*fwd_z);
  if(G4FST::SETTING::EPIC_TRACKINGGEO_VARIANT)
  {
    n_bkwd_disk = sizeof(fwd_z) / sizeof(*fwd_z);
  } else {
    n_bkwd_disk = sizeof(fwd_z_EPIC) / sizeof(*fwd_z_EPIC);
  }
  for (unsigned int i = 0; i < n_fwd_disk; i++)
  {
    if(Enable::EPIC_TRACKINGGEO){
      string name = Form("FST_%i", i);
      double pitch = 10e-4;
      PHG4FSTSubsystem *hfst = new PHG4FSTSubsystem(name);
      hfst->OverlapCheck(OverlapCheck);
      hfst->SuperDetector(name);
      hfst->SetActive();
      if (AbsorberActive) hfst->SetAbsorberActive();
      hfst->set_double_param("z_position", G4FST::SETTING::EPIC_TRACKINGGEO_VARIANT ? fwd_z_EPIC2[i] : fwd_z_EPIC[i]);
      hfst->set_double_param("r_min", G4FST::SETTING::EPIC_TRACKINGGEO_VARIANT ? fwd_rmin_EPIC2[i] : fwd_rmin_EPIC[i]);
      hfst->set_double_param("r_max", G4FST::SETTING::EPIC_TRACKINGGEO_VARIANT ? fwd_rmax_EPIC2[i] : fwd_rmax_EPIC[i]);
      hfst->set_double_param("silicon_thickness", 35 * um);
      hfst->set_double_param("offset_cutout", G4FST::SETTING::EPIC_TRACKINGGEO_VARIANT ? fwd_offset_EPIC2[i] : fwd_offset_EPIC[i]);
      g4Reco->registerSubsystem(hfst);

      if (TRACKING::FastKalmanFilter)
      {
        TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                                PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                                pitch / sqrt(12.),                 //      const float radres,
                                                pitch / sqrt(12.),                 //      const float phires,
                                                50e-4 / sqrt(12.),                 //      const float lonres, *ignored in plane detector*
                                                G4FST::SETTING::TRACKING_EFFICIENCY,                                 //      const float eff,
                                                0);                                //      const float noise
        TRACKING::FastKalmanFilterInnerTrack->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                                PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                                pitch / sqrt(12.),                 //      const float radres,
                                                pitch / sqrt(12.),                 //      const float phires,
                                                50e-4 / sqrt(12.),                 //      const float lonres, *ignored in plane detector*
                                                G4FST::SETTING::TRACKING_EFFICIENCY,                                 //      const float eff,
                                                0);                                //      const float noise
        TRACKING::FastKalmanFilterSiliconTrack->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                                PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                                pitch / sqrt(12.),                 //      const float radres,
                                                pitch / sqrt(12.),                 //      const float phires,
                                                50e-4 / sqrt(12.),                 //      const float lonres, *ignored in plane detector*
                                                G4FST::SETTING::TRACKING_EFFICIENCY,                                 //      const float eff,
                                                0);                                //      const float noise
      }
    } else if(Enable::AI_TRACKINGGEO){
      make_LANL_FST_station(Form("FST_%i", i), g4Reco, fwd_z[i], fwd_rmin[i], fwd_rmax[i], 35 * um, 10e-4);  //cm
    } else {
      make_LANL_FST_station(Form("FST_%i", i), g4Reco, fwd_z_np[i], fwd_rmin_np[i], fwd_rmax_np[i], 35 * um, 10e-4);  //cm
    }
  }

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
  // fst->OverlapCheck(true);  //true);//overlapcheck);

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
                                             G4FST::SETTING::TRACKING_EFFICIENCY,                                 //      const float eff,
                                             0);                                //      const float noise
    TRACKING::FastKalmanFilterInnerTrack->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             pitch / sqrt(12.),                 //      const float radres,
                                             pitch / sqrt(12.),                 //      const float phires,
                                             50e-4 / sqrt(12.),                 //      const float lonres, *ignored in plane detector*
                                             G4FST::SETTING::TRACKING_EFFICIENCY,                                 //      const float eff,
                                             0);                                //      const float noise
    TRACKING::FastKalmanFilterSiliconTrack->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             pitch / sqrt(12.),                 //      const float radres,
                                             pitch / sqrt(12.),                 //      const float phires,
                                             50e-4 / sqrt(12.),                 //      const float lonres, *ignored in plane detector*
                                             G4FST::SETTING::TRACKING_EFFICIENCY,                                 //      const float eff,
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

