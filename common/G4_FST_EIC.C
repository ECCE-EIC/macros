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

/* 

The idea is to calculate all the rmax rmin based on z and theta1 and theta2 in a python script (obj_fun.py)
Once calculated just pass the disk dimensions into this script

*/

void FSTSetup(PHG4Reco *g4Reco)
{
  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = .1 * cm;
  const double um = 1e-3 * mm;
  //const double max_radius = 50.;



//  const double bkwd_z[] = {33.2, 58.29, 80.05, 107.4};
//  double bkwd_rmin[] = {3.3, 3.3, 5.25, 6.4};
//  double bkwd_rmax[] = {15.3, 27.3, 35.25, 48.4};
  const double bkwd_z[] = {25, 52, 79, 106};
  double bkwd_rmin[] = {3.5, 3.5, 4.5, 5.5};
  double bkwd_rmax[] = {18.5, 36.5, 40.5, 41.5};
  const int n_bkwd_disk = sizeof(bkwd_z) / sizeof(*bkwd_z);
  for (unsigned int i = 0; i < n_bkwd_disk; i++)
  {
/*
    // Below was made to auto calculate the min and max Radius
    if(bkwd_z[i] < uRwell1_e_length) bkwd_rmax[i] = std::min(max_radius, e_slope1*bkwd_z[i] + e_intercept1) - 0.5;
    else if (bkwd_z[i] >= uRwell1_e_length && bkwd_z[i] <= (uRwell1_e_length + uRwell_plateau_length)){bkwd_rmax[i] = uRwell1_radius - 1.5;}
    else if(bkwd_z[i] > (uRwell1_e_length + uRwell_plateau_length)){bkwd_rmax[i] = std::min(max_radius, e_slope2*bkwd_z[i] + e_intercept2) - 0.5;}
    else {cout << "Cannot calculate the RMax exiting" << endl; gSystem->Exit(0); }

    if(bkwd_z[i]>79.8 && bkwd_z[i]>0) bkwd_rmin[i] = (0.0521*bkwd_z[i] + 1.0);
    else bkwd_rmin[i] = 3.3;
*/  
    make_LANL_FST_station(Form("EST_%i", i), g4Reco, -1*bkwd_z[i], bkwd_rmin[i], bkwd_rmax[i], 35 * um, 10e-4);  //cm
  }


  const double fwd_z[] = {25, 52, 73, 106, 125};
  double fwd_rmin[] = {3.5, 3.5, 4.5, 5.5, 7.5};
  double fwd_rmax[] = {18.5, 36.5, 40.5, 41.5, 43.4};
  const int n_fwd_disk = sizeof(fwd_z) / sizeof(*fwd_z);
  for (unsigned int i = 0; i < n_fwd_disk; i++)
  {

/*
    if(fwd_z[i] < uRwell1_h_length) fwd_rmax[i] = std::min(max_radius, h_slope1*fwd_z[i] + h_intercept1) - 0.5;
    else if (fwd_z[i] >= uRwell1_h_length && fwd_z[i] <= (uRwell1_h_length + uRwell_plateau_length)){fwd_rmax[i] = uRwell1_radius - 1.5;}
    //else if(fwd_z[i] > (uRwell1_h_length + uRwell_plateau_length)){fwd_rmax[i] = std::min(max_radius, h_slope2*fwd_z[i] + h_intercept2) - 0.5;}
    else if(fwd_z[i] > (uRwell1_h_length + uRwell_plateau_length) && fwd_z[i] <= 130.){ fwd_rmax[i] = std::min(max_radius, h_slope2*fwd_z[i] + h_intercept2) - 0.5;}
    //else if(fwd_z[i] > 113.){fwd_rmax[i] = h_slope2*(fwd_z[i] + 5.) + h_intercept2 - 0.25;}
    else {cout << "Cannot calculate the RMax exiting" << endl; gSystem->Exit(0); }

    if(fwd_z[i]>66.8 && fwd_z[i]>0) fwd_rmin[i] = (0.0521*fwd_z[i] + 1.0);
    else fwd_rmin[i] = 3.3;
*/
    make_LANL_FST_station(Form("FST_%i", i), g4Reco, fwd_z[i], fwd_rmin[i], fwd_rmax[i], 35 * um, 10e-4);  //cm
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
                                             0.9,                                 //      const float eff,
                                             0);                                //      const float noise
    TRACKING::FastKalmanFilterInnerTrack->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             pitch / sqrt(12.),                 //      const float radres,
                                             pitch / sqrt(12.),                 //      const float phires,
                                             50e-4 / sqrt(12.),                 //      const float lonres, *ignored in plane detector*
                                             0.9,                                 //      const float eff,
                                             0);                                //      const float noise
    TRACKING::FastKalmanFilterSiliconTrack->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             pitch / sqrt(12.),                 //      const float radres,
                                             pitch / sqrt(12.),                 //      const float phires,
                                             50e-4 / sqrt(12.),                 //      const float lonres, *ignored in plane detector*
                                             0.9,                                 //      const float eff,
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

