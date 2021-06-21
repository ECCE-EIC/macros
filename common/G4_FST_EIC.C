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
//#include <g4mvtx/PHG4EICMvtxSubsystem.h>
#include <g4main/PHG4Reco.h>

#include <string>

R__LOAD_LIBRARY(libg4detectors.so)

int make_LANL_FST_station(string name, PHG4Reco *g4Reco, double zpos, double Rmin,
                          double Rmax, double tSilicon, double pitch);
int make_supportCyl(string name, PHG4Reco *g4Reco,
                    double r, double t, double length);
double barrelService(PHG4Reco* g4Reco, double radius);
void calculateMaterialBoundaries(int& service_layer_ID, double& outer_copper_radius,
                                 double& outer_water_radius, double& outer_plastic_radius);
double calculateOR(double inner_radius, double area);
//-----------------------------------------------------------------------------------//
namespace Enable
{
  static bool FST = false;
  bool FST_OVERLAPCHECK = false;
  bool FST_ABSORBER = false;
  int FST_VERBOSITY = 0;
}  // namespace Enable

namespace G4FST
{
  namespace SETTING
  {
    bool FST_TPC = false;
    bool SUPPORTCYL = false;
    
    //supporting structure parameters
    const int n_service_layers = 1;               //Number of service cable service_layers to generate 
    double service_layer_start_radius[] = {55.1}; //Inner radius of where the cables begin [cm] 
    int n_staves_service_layer[] = {1};//{48};    //Number of staves associated to each service layer 
    double service_barrel_radius = 55.6;          // [cm] From final design review 
    double service_barrel_start = -125.0*2.0;           //[cm] Approx.
    double service_barrel_length = 125.0*2.0;           // [cm] length of service barrel ~(to patch panel)

    double single_stave_service_copper_area = 0.0677;   //Cross-sectional area of copper for 1 stave [cm^2]
    double single_stave_service_water_area = 0.0098;    //Cross-sectional area of water for 1 stave [cm^2]
    double single_stave_service_plastic_area = 0.4303;  //Cross-sectional area of plastic for 1 stave [cm^2]
  }  // namespace SETTING
}  // namespace G4FST

//-----------------------------------------------------------------------------------//
void FST_Init()
{
  if (!Enable::FST_OVERLAPCHECK) 
  {
    Enable::FST_OVERLAPCHECK=Enable::OVERLAPCHECK;
  }

  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 48.0);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 125.0);
  if (G4FST::SETTING::SUPPORTCYL)
  {
    BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, 127.0);
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
    Enable::FST_OVERLAPCHECK = Enable::OVERLAPCHECK;
  }
  //bool OverlapCheck = Enable::OVERLAPCHECK || Enable::FST_OVERLAPCHECK;

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
    double tSupport = 0.5;                                                        //cm
    //make_supportCyl("FSTSupportCyl", g4Reco, 50.1 + gap, tSupport, 125.0 * 2.0);  //cm
    barrelService(g4Reco, 125);
  }


  
}
//-----------------------------------------------------------------------------------//
int make_LANL_FST_station(string name, PHG4Reco *g4Reco,
                          double zpos, double Rmin, double Rmax, double tSilicon, double pitch)  //silicon thickness
{

  double min_polar_angle =atan2(Rmin, zpos) ;
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
  fst->OverlapCheck(Enable::FST_OVERLAPCHECK);  //true);//overlapcheck);

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
  PHG4CylinderSubsystem *cyl = new PHG4CylinderSubsystem(name, 5);
  cyl->set_double_param("radius", r);
  cyl->set_double_param("length", length);
  cyl->set_string_param("material", "G4_GRAPHITE");
  cyl->set_double_param("thickness", t);
  cyl->set_double_param("place_x", 0.);
  cyl->set_double_param("place_y", 0.);
  cyl->set_double_param("place_z", 0);
  cyl->SetActive(0);
  //cyl->SuperDetector("");
  cyl->OverlapCheck(Enable::FST_OVERLAPCHECK);  //OverlapCheck);

  g4Reco->registerSubsystem(cyl);
  return 0;
}
//-----------------------------------------------------------------------------------/
//modified from MVTX macro (April 5th 2021 version)
double barrelService(PHG4Reco* g4Reco, double radius)
{
  // Note, cables are all south
  // Setup service_layers
  bool AbsorberActive = Enable::ABSORBER || Enable::FST_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::FST_OVERLAPCHECK;
  int verbosity = std::max(Enable::VERBOSITY, Enable::FST_VERBOSITY);

  double copper_OR[G4FST::SETTING::n_service_layers]; //Objects for material outer radii
  double water_OR[G4FST::SETTING::n_service_layers];
  double plastic_OR[G4FST::SETTING::n_service_layers];  //Objects for material outer radii
  int subsystem_service_layer = 0;
  std::string copper_name, water_name, plastic_name;
  PHG4CylinderSubsystem* cyl;

  for (int i = 0; i < G4FST::SETTING::n_service_layers; ++i)  //Build a service_layer of copper, then water, then plastic
    {
      calculateMaterialBoundaries(i, copper_OR[i], water_OR[i], plastic_OR[i]);

      copper_name = "barrel_Service_copper_service_layer_" + std::to_string(i);
      water_name = "barrel_Service_water_service_layer_" + std::to_string(i);
      plastic_name = "barrel_Service_plastic_service_layer_" + std::to_string(i);

      cyl = new PHG4CylinderSubsystem(copper_name, subsystem_service_layer);
      cyl->set_double_param("place_z", -1 * (G4FST::SETTING::service_barrel_length + G4FST::SETTING::service_barrel_start) - no_overlapp);
      cyl->set_double_param("radius", G4FST::SETTING::service_layer_start_radius[i]);
      cyl->set_int_param("lengthviarapidity", 0);
      cyl->set_double_param("length", G4FST::SETTING::service_barrel_length);
      cyl->set_string_param("material", "G4_Cu");
      cyl->set_double_param("thickness", copper_OR[i] - G4FST::SETTING::service_layer_start_radius[i]);
      cyl->SuperDetector("barrelSERVICE");
      if (AbsorberActive) cyl->SetActive();
      cyl->OverlapCheck(Enable::FST_OVERLAPCHECK);
      g4Reco->registerSubsystem(cyl);
      subsystem_service_layer += 1;

      cyl = new PHG4CylinderSubsystem(water_name, subsystem_service_layer);
      cyl->set_double_param("place_z", -1 * (G4FST::SETTING::service_barrel_length + G4FST::SETTING::service_barrel_start) - no_overlapp);
      cyl->set_double_param("radius", copper_OR[i]);
      cyl->set_int_param("lengthviarapidity", 0);
      cyl->set_double_param("length", G4FST::SETTING::service_barrel_length);
      cyl->set_string_param("material", "G4_WATER");
      cyl->set_double_param("thickness", water_OR[i] - copper_OR[i]);
      cyl->SuperDetector("barrelSERVICE");
      if (AbsorberActive) cyl->SetActive();
      cyl->OverlapCheck(Enable::FST_OVERLAPCHECK);
      g4Reco->registerSubsystem(cyl);
      subsystem_service_layer += 1;

      cyl = new PHG4CylinderSubsystem(plastic_name, subsystem_service_layer);
      cyl->set_double_param("place_z", -1 * (G4FST::SETTING::service_barrel_length + G4FST::SETTING::service_barrel_start) - no_overlapp);
      cyl->set_double_param("radius", water_OR[i]);
      cyl->set_int_param("lengthviarapidity", 0);
      cyl->set_double_param("length", G4FST::SETTING::service_barrel_length);
      cyl->set_string_param("material", "G4_POLYETHYLENE");
      cyl->set_double_param("thickness", plastic_OR[i] - water_OR[i]);
      cyl->SuperDetector("barrelSERVICE");
      if (AbsorberActive) cyl->SetActive();
      cyl->OverlapCheck(Enable::FST_OVERLAPCHECK);
      g4Reco->registerSubsystem(cyl);
      subsystem_service_layer += 1;
    }

  cyl = new PHG4CylinderSubsystem("barrel_Service_shell_service_layer", subsystem_service_layer);
  cyl->set_double_param("place_z", -1 * (G4FST::SETTING::service_barrel_length + G4FST::SETTING::service_barrel_start) - no_overlapp);
  cyl->set_double_param("radius", G4FST::SETTING::service_barrel_radius);
  cyl->set_int_param("lengthviarapidity", 0);
  cyl->set_double_param("length", G4FST::SETTING::service_barrel_length);
  cyl->set_string_param("material", "PEEK");  //Service barrel is carbon fibre (peek?)
  cyl->set_double_param("thickness", 0.1);    //Service barrel is 1mm thick
  cyl->SuperDetector("barrelSERVICE");
  if (AbsorberActive) cyl->SetActive();
  cyl->OverlapCheck(Enable::FST_OVERLAPCHECK);
  g4Reco->registerSubsystem(cyl);

  radius = G4FST::SETTING::service_barrel_radius;

  if (verbosity > 0)
    {
      cout << "=========================== barrel Service Barrel =============================" << endl;
      cout << " barrel Service Material Description:" << endl;

      //cout << "  Single stave copper area  = " << G4FST::SETTING::single_stave_service_copper_area << " cm^2" << endl;
      //cout << "  Single stave water area   = " << G4FST::SETTING::single_stave_service_water_area << " cm^2" << endl;
      //cout << "  Single stave plastic area = " << G4FST::SETTING::single_stave_service_plastic_area << " cm^2" << endl;

      for (int j = 0; j < G4FST::SETTING::n_service_layers; ++j)
	{
	  cout << "  service_layer " << j << " starts at " << G4FST::SETTING::service_layer_start_radius[j] << " cm" << endl;
	  cout << "  service_layer " << j << " services " << G4FST::SETTING::n_staves_service_layer[j] << " staves" << endl;
	}

      cout << "  Service barrel radius = " << G4FST::SETTING::service_barrel_radius << " cm" << endl;
      cout << "  Service barrel start = " << G4FST::SETTING::service_barrel_start << " cm" << endl;
      cout << "  Service barrel length = " << G4FST::SETTING::service_barrel_length << " cm" << endl;
      cout << "===============================================================================" << endl;
    }

  radius += no_overlapp;

  return radius;
}
//-----------------------------------------------------------------------------------//
void calculateMaterialBoundaries(int& service_layer_ID, double& outer_copper_radius, 
				 double& outer_water_radius, double& outer_plastic_radius)  
//Calculate where the transition between each material occurs
{
  outer_copper_radius = calculateOR(G4FST::SETTING::service_layer_start_radius[service_layer_ID], 
				    G4FST::SETTING::n_staves_service_layer[service_layer_ID] 
				    * G4FST::SETTING::single_stave_service_copper_area);
  outer_water_radius = calculateOR(outer_copper_radius, 
				   G4FST::SETTING::n_staves_service_layer[service_layer_ID] 
				   * G4FST::SETTING::single_stave_service_water_area);
  outer_plastic_radius = calculateOR(outer_water_radius, 
				     G4FST::SETTING::n_staves_service_layer[service_layer_ID] 
				     * G4FST::SETTING::single_stave_service_plastic_area);
}
//-----------------------------------------------------------------------------------//
double calculateOR(double inner_radius, double area)  
//Calculate the outer radius of a disk, knowing the inner radius and the area
{
  return std::sqrt(area / M_PI + std::pow(inner_radius, 2));
}
#endif
