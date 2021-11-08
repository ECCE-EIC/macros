#ifndef MACRO_G4TrackingService_C
#define MACRO_G4TrackingService_C

#include <GlobalVariables.C>
#include <QA.C>

#include <g4detectors/PHG4ConeSubsystem.h>
#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4main/PHG4Reco.h>

#include <qa_modules/QAG4SimulationMvtx.h>

#include <fun4all/Fun4AllServer.h>

#include <cmath>
#include <vector>

//ECCE Tracking Services
//Should be 5 barrels and 4 cones

using namespace std;

class ServiceProperties
{
 public:
  ServiceProperties();

  explicit ServiceProperties(const string &name,
                             const double &rad_len_copper,
                             const double &rad_len_aluminum,
                             const double &rad_len_water,
                             const double &rad_len_plastic,
                             const double &rad_len_carbon,
                             const double &rad_len_iron,
                             const double &z_south,
                             const double &z_north,
                             const double &r_south,
                             const double &r_north);

  virtual ~ServiceProperties(){};

  const string get_name();
  const double get_rad_len_copper();
  const double get_rad_len_aluminum();
  const double get_rad_len_water();
  const double get_rad_len_plastic();
  const double get_rad_len_carbon();
  const double get_rad_len_iron();
  const double get_z_south();
  const double get_z_north();
  const double get_r_south();
  const double get_r_north();

 private:
  const string m_name = "service";
  const double m_rad_len_copper = 0.0;
  const double m_rad_len_aluminum = 0.0;
  const double m_rad_len_water = 0.0;
  const double m_rad_len_plastic = 0.0;
  const double m_rad_len_carbon = 0.0;
  const double m_rad_len_iron = 0.0;
  const double m_z_south = 0.0;
  const double m_z_north = 0.0;
  const double m_r_south = 0.0;
  const double m_r_north = 0.0;
};
//8, 0, 0.42, 0.27, shellX0, 0
ServiceProperties::ServiceProperties(const string &name,
                                     const double &rad_len_copper,
                                     const double &rad_len_aluminum,
                                     const double &rad_len_water,
                                     const double &rad_len_plastic,
                                     const double &rad_len_carbon,
                                     const double &rad_len_iron,
                                     const double &z_south,
                                     const double &z_north,
                                     const double &r_south,
                                     const double &r_north)
  : m_name(name)
  , m_rad_len_copper(rad_len_copper)
  , m_rad_len_aluminum(rad_len_aluminum)
  , m_rad_len_water(rad_len_water)
  , m_rad_len_plastic(rad_len_plastic)
  , m_rad_len_carbon(rad_len_carbon)
  , m_rad_len_iron(rad_len_iron)
  , m_z_south(z_south)
  , m_z_north(z_north)
  , m_r_south(r_south)
  , m_r_north(r_north)
{
}

const string ServiceProperties::get_name() { return m_name; }
const double ServiceProperties::get_rad_len_copper() { return m_rad_len_copper; }
const double ServiceProperties::get_rad_len_aluminum() { return m_rad_len_aluminum; }
const double ServiceProperties::get_rad_len_water() { return m_rad_len_water; }
const double ServiceProperties::get_rad_len_plastic() { return m_rad_len_plastic; }
const double ServiceProperties::get_rad_len_carbon() { return m_rad_len_carbon; }
const double ServiceProperties::get_rad_len_iron() { return m_rad_len_iron; }
const double ServiceProperties::get_z_south() { return m_z_south; }
const double ServiceProperties::get_z_north() { return m_z_north; }
const double ServiceProperties::get_r_south() { return m_r_south; }
const double ServiceProperties::get_r_north() { return m_r_north; }

namespace Enable
{
  bool TrackingService = false;
  bool TrackingService_ABSORBER = false;
  bool TrackingService_OVERLAPCHECK = false;
  int TrackingService_VERBOSITY = 0;

}  // namespace Enable

namespace G4TrackingService
{  //List materials and radiation length in cm
  const int nMaterials = 6;
  pair<string, double> materials[nMaterials] = {make_pair("G4_Cu", 1.436)
          				       ,make_pair("G4_Al", 8.897)
          				       ,make_pair("G4_WATER", 36.08)
          				       ,make_pair("G4_POLYETHYLENE", 50.31)
          				       ,make_pair("PEEK", 30.00)
          				       ,make_pair("G4_Fe", 1.757)};

  double GlobalOffset = 0.0;
  double ShellThickness = 0.3;  //Thickness in cm
  int subsysID = 0;
}  // namespace G4TrackingService

vector<double> get_thickness(ServiceProperties *object)
{
  vector<double> thickness = {(object->get_rad_len_copper() / 100) * G4TrackingService::materials[0].second
                             ,(object->get_rad_len_aluminum() / 100) * G4TrackingService::materials[1].second
                             ,(object->get_rad_len_water() / 100) * G4TrackingService::materials[2].second
                             ,(object->get_rad_len_plastic() / 100) * G4TrackingService::materials[3].second
                             ,(object->get_rad_len_carbon() / 100) * G4TrackingService::materials[4].second
                             ,(object->get_rad_len_iron() / 100) * G4TrackingService::materials[5].second};
  return thickness;
}

void TrackingServiceInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 280.);
  // extends only to -z
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -450.);
}

double TrackingServiceCone(ServiceProperties *object, PHG4Reco *g4Reco, double radius)
{
  bool AbsorberActive = Enable::ABSORBER || Enable::TrackingService_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::TrackingService_OVERLAPCHECK;
  int verbosity = max(Enable::VERBOSITY, Enable::TrackingService_VERBOSITY);

  PHG4ConeSubsystem *cone;

  double innerRadiusSouth = object->get_r_south();
  double innerRadiusNorth = object->get_r_north();
  double length = abs(object->get_z_north() - object->get_z_south());
  vector<double> thickness = get_thickness(object);

  for (int i = 0; i < G4TrackingService::nMaterials; ++i)
  {
    if (thickness[i] == 0) continue;
    cone = new PHG4ConeSubsystem(object->get_name(), G4TrackingService::subsysID);
    cone->Verbosity(verbosity);
    cone->SetR1(innerRadiusSouth, innerRadiusSouth + thickness[i]);
    cone->SetR2(innerRadiusNorth, innerRadiusNorth + thickness[i]);
    cone->SetPlaceZ(object->get_z_south() + length / 2 + G4TrackingService::GlobalOffset);
    cone->SetZlength(length / 2);
    cone->SetMaterial(G4TrackingService::materials[i].first);
    cone->SuperDetector("TrackingService");
    if (AbsorberActive) cone->SetActive();
    cone->OverlapCheck(OverlapCheck);
    g4Reco->registerSubsystem(cone);
    ++G4TrackingService::subsysID;
    innerRadiusSouth += thickness[i];
    innerRadiusNorth += thickness[i];
  }
  radius = max(innerRadiusSouth, innerRadiusNorth);

  return radius;
}

double TrackingServiceCylinder(ServiceProperties *object, PHG4Reco *g4Reco, double radius)
{
  bool AbsorberActive = Enable::ABSORBER || Enable::TrackingService_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::TrackingService_OVERLAPCHECK;
  int verbosity = max(Enable::VERBOSITY, Enable::TrackingService_VERBOSITY);

  PHG4CylinderSubsystem *cyl;

  double innerRadius = object->get_r_south();
  double length = abs(object->get_z_north() - object->get_z_south());
  vector<double> thickness = get_thickness(object);

  for (int i = 0; i < G4TrackingService::nMaterials; ++i)
  {
    if (thickness[i] == 0) continue;
    cyl = new PHG4CylinderSubsystem(object->get_name(), G4TrackingService::subsysID);
    cyl->Verbosity(verbosity);
    cyl->set_double_param("place_z", object->get_z_south() + length / 2 + G4TrackingService::GlobalOffset);
    cyl->set_double_param("radius", innerRadius);
    cyl->set_double_param("length", length);
    cyl->set_string_param("material", G4TrackingService::materials[i].first);
    cyl->set_double_param("thickness", thickness[i]);
    cyl->SuperDetector("TrackingService");
    if (AbsorberActive) cyl->SetActive();
    cyl->OverlapCheck(OverlapCheck);
    g4Reco->registerSubsystem(cyl);
    ++G4TrackingService::subsysID;
    innerRadius += thickness[i];
  }
  radius = innerRadius;

  return radius;
}

double TrackingService(PHG4Reco *g4Reco, double radius)
{
  vector<ServiceProperties *> cylinders, cones;

  double shellX0 = 100 * G4TrackingService::ShellThickness / G4TrackingService::materials[4].second;
  double disk_cone_radii = 50.;  // max radius for the inner tracker Thickness incorporated in here
  double avg_thickness_inner = 1.0;  // for 13, 0, 0.56, 0.48, shellX0
  double vtx_support_radius = 9.0;
  double vtx_e_length = 19.44;
  double vtx_h_length = 19.44;
  double sagitta_support_radius = 17.0 + 1.5;  // 17.0 cms is the radiius of the second sagitta layer. 1.5cms is the gap given for support structure
  double sagitta_support_e_length = 37.7861;  // +1.5 is the gap given for the support structure
  double sagitta_support_h_length = 37.7861;
  double plateau_length = 5.1;

  /* Inclination 2 for e-going and h-going directions*/

  double inner_uRwell_radius = 33.12;  // Radius of the inner uRwell
  double gap = 1. * avg_thickness_inner;
  inner_uRwell_radius -= gap;  // subtract the thickness
  double inner_uRwell_e_length = 67.54;
  double inner_uRwell_h_length = 67.54;
  double inner_uRwell_length = inner_uRwell_e_length + inner_uRwell_h_length;

  double e_cone_ends = 110.;
  double h_cone_ends = 110.;

  /* Below are the support structures beyond EGEM*/
  cylinders.push_back(new ServiceProperties("ETrackingCylinderService_1", 9, 0, 0.42, 0.32, shellX0, 0, -400, -310, 270, 0));
  cones.push_back(new ServiceProperties("ETrackingConeService_1", 9, 0, 0.42, 0.32, shellX0, 0, -310, -300, 270, 68));
  cylinders.push_back(new ServiceProperties("ETrackingCylinderService_2", 17, 0, 0.56, 0.64, shellX0, 0, -300, -200, 68, 0));
  cylinders.push_back(new ServiceProperties("ETrackingCylinderService_3", 15, 0, 0.56, 0.56, shellX0, 0, -200, -152, 68, 0));

  cones.push_back(new ServiceProperties("ETrackingConeService_2", 13, 0, 0.56, 0.48, shellX0, 0, -152, -132.1, 68, 63.2));
  cones.push_back(new ServiceProperties("ETrackingConeService_3", 13, 0, 0.56, 0.48, shellX0, 0, -132.1, -120., 63.2, disk_cone_radii));

  /* Support structure thickness*/

  double CuThickness = 13.;          // 0.18668 cms
  double AlThickness = 0.;           // 0. cms
  double WaterThickness = 0.70;      // 0.25256 cms
  double PlasticThickness = 0.48;    // 0.241488 cms
  double CarbonThickness = shellX0;  // 0.3 cms

  // Cylinder from end Disk to EGEM
  cylinders.push_back(new ServiceProperties("ETrackingCyl_EGEMToDisk4", CuThickness, AlThickness, WaterThickness, PlasticThickness, CarbonThickness, 0, -120.0, -1 * e_cone_ends, disk_cone_radii, 0));

  // Cone service from the end Disk to the uRwell1 radius
  cones.push_back(new ServiceProperties("ETrackingCone_Disk4TouRwell", CuThickness, AlThickness, WaterThickness, PlasticThickness, CarbonThickness, 0, -1 * e_cone_ends, -1 * inner_uRwell_e_length - plateau_length, disk_cone_radii, inner_uRwell_radius));

  // The cylindrical plateau structure for Electron uRwell side.
  cylinders.push_back(new ServiceProperties("ETrackingCyl_uRWellPlateau", CuThickness, AlThickness, WaterThickness, PlasticThickness, CarbonThickness, 0, -1 * inner_uRwell_e_length - plateau_length, -1 * inner_uRwell_e_length, inner_uRwell_radius, 0));

  // Cone service from electron side uRwell1 to vertex support.
  cones.push_back(new ServiceProperties("ETrackingCone_uRwellToVertex", CuThickness, AlThickness, WaterThickness, PlasticThickness, CarbonThickness, 0, -1 * inner_uRwell_e_length, -1 * vtx_e_length, inner_uRwell_radius, vtx_support_radius));

  // Sagitta Cylindrical Support Structure
  cylinders.push_back(new ServiceProperties("BTrackingCyl_Sagitta", 0, 0, 0, 0, 0.1, 0, -1 * sagitta_support_e_length + gap, sagitta_support_h_length - gap, sagitta_support_radius + gap / 4, 0));

  // Vertex Cylindrical Support Structure
  cylinders.push_back(new ServiceProperties("BTrackingCyl_Vertex", 0, 0, 0, 0, shellX0, 0, -1 * vtx_e_length, vtx_h_length, vtx_support_radius, 0));

  // Cone service in H-region from vertex to inner uRwell
  cones.push_back(new ServiceProperties("HTrackingCone_VertexTouRwell", CuThickness, AlThickness, WaterThickness, PlasticThickness, CarbonThickness, 0, vtx_h_length, inner_uRwell_h_length, vtx_support_radius, inner_uRwell_radius));

  // Cylinder service to rest the uRwell in H region Plateau
  cylinders.push_back(new ServiceProperties("HTrackingRWellPlateau", CuThickness, AlThickness, WaterThickness, PlasticThickness, CarbonThickness, 0, inner_uRwell_h_length, inner_uRwell_h_length + plateau_length, inner_uRwell_radius, 0));

  // Cone service from uRwell to Disk 5 in h-region
  cones.push_back(new ServiceProperties("HTrackingCone_uRwellToDisk5", CuThickness, AlThickness, WaterThickness, PlasticThickness, CarbonThickness, 0, inner_uRwell_h_length + plateau_length, h_cone_ends, inner_uRwell_radius, disk_cone_radii));

  // Cylinder service from Disk 5 to 137 in h direction. This 137 cms is from 2nd campaign and keeping the outer tracker intact.
  cylinders.push_back(new ServiceProperties("HSidemRWellSupportCyl", CuThickness, AlThickness, WaterThickness, PlasticThickness, CarbonThickness, 0, h_cone_ends, 124., disk_cone_radii, 0));

  // Cone service from 137. to FGEM
  //cones.push_back(new ServiceProperties("HTrackingConeForFGEM", CuThickness, AlThickness, WaterThickness, PlasticThickness, CarbonThickness, 0, 137., 157.9, disk_cone_radii, 69.2));
  cones.push_back(new ServiceProperties("HTrackingConeForFGEM", CuThickness, AlThickness, WaterThickness, PlasticThickness, CarbonThickness, 0, 124., 173, disk_cone_radii, 69.2));

  // Supports beyond FGEM
  //cylinders.push_back(new ServiceProperties("HTrackingCylinderService_1", 15, 0, 0.84, 0.56, shellX0, 0, 157.9, 173, 69.2, 0));
  cones.push_back(new ServiceProperties("HTrackingConeService_7", 13, 0, 0.70, 0.48, shellX0, 0, 173, 180, 69.2, 85));
  cones.push_back(new ServiceProperties("HTrackingConeService_8", 13, 0, 0.70, 0.48, shellX0, 0, 180, 195, 85, 100));

  cylinders.push_back(new ServiceProperties("EEMCalSupport", 0, 0, 0, 0, 0, 171, -200, -197, 62, 0));

  for (ServiceProperties *cylinder : cylinders) radius += TrackingServiceCylinder(cylinder, g4Reco, radius);
  for (ServiceProperties *cone : cones) radius += TrackingServiceCone(cone, g4Reco, radius);

  return radius;
}

#endif
