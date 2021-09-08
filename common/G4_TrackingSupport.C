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
    const double m_z_south = 0.0;
    const double m_z_north = 0.0;
    const double m_r_south = 0.0;
    const double m_r_north = 0.0;
};

ServiceProperties::ServiceProperties(const string &name,
                                     const double &rad_len_copper,
                                     const double &rad_len_aluminum,
                                     const double &rad_len_water,
                                     const double &rad_len_plastic,
                                     const double &rad_len_carbon,
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
  , m_z_south(z_south)
  , m_z_north(z_north)
  , m_r_south(r_south)
  , m_r_north(r_north)
{}

const string ServiceProperties::get_name() { return m_name; }
const double ServiceProperties::get_rad_len_copper() { return m_rad_len_copper; }
const double ServiceProperties::get_rad_len_aluminum() { return m_rad_len_aluminum; }
const double ServiceProperties::get_rad_len_water() { return m_rad_len_water; }
const double ServiceProperties::get_rad_len_plastic() { return m_rad_len_plastic; }
const double ServiceProperties::get_rad_len_carbon() { return m_rad_len_carbon; }
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
{ //List materials and radiation length in cm
  const int nMaterials = 5;
  pair<string, double> materials[nMaterials] = { make_pair("G4_Cu", 1.436)
                                               , make_pair("G4_Al",  8.897)
                                               , make_pair("G4_WATER",  36.08)
                                               , make_pair("G4_POLYETHYLENE", 50.31)
                                               , make_pair("PEEK", 30.00) };
  
  double GlobalOffset = 0.0;
  int subsysID = 0;
}  // namespace G4TrackingService

vector<double> get_thickness(ServiceProperties *object)
{
  vector<double> thickness = {(object->get_rad_len_copper()/100)*G4TrackingService::materials[0].second
                             ,(object->get_rad_len_aluminum()/100)*G4TrackingService::materials[1].second
                             ,(object->get_rad_len_water()/100)*G4TrackingService::materials[2].second
                             ,(object->get_rad_len_plastic()/100)*G4TrackingService::materials[3].second
                             ,(object->get_rad_len_carbon()/100)*G4TrackingService::materials[4].second};
  return thickness;
}

void TrackingServiceInit()
{
}

double TrackingServiceCone(ServiceProperties *object, PHG4Reco* g4Reco, double radius)
{
  bool AbsorberActive = Enable::ABSORBER || Enable::TrackingService_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::TrackingService_OVERLAPCHECK;
  int verbosity = max(Enable::VERBOSITY, Enable::TrackingService_VERBOSITY);

  PHG4ConeSubsystem* cone;

  double innerRadiusSouth = object->get_r_south();
  double innerRadiusNorth = object->get_r_north();
  double length = abs(object->get_z_north() - object->get_z_south());
  vector<double> thickness = get_thickness(object);

  for (int i = 0; i < G4TrackingService::nMaterials; ++i)
  {
    if (thickness[i] == 0) continue;
    cone = new PHG4ConeSubsystem(object->get_name(), G4TrackingService::subsysID);
    cone->SetR1(innerRadiusSouth, innerRadiusSouth + thickness[i]);
    cone->SetR2(innerRadiusNorth, innerRadiusNorth + thickness[i]);
    cone->SetPlaceZ(object->get_z_south() + length/2 + G4TrackingService::GlobalOffset);
    cone->SetZlength(length);
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

double TrackingServiceCylinder(ServiceProperties *object, PHG4Reco* g4Reco, double radius)
{
  bool AbsorberActive = Enable::ABSORBER || Enable::TrackingService_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::TrackingService_OVERLAPCHECK;
  int verbosity = max(Enable::VERBOSITY, Enable::TrackingService_VERBOSITY);

  PHG4CylinderSubsystem* cyl;

  double innerRadius = object->get_r_south();
  double length = abs(object->get_z_north() - object->get_z_south());
  vector<double> thickness = get_thickness(object);

  for (int i = 0; i < G4TrackingService::nMaterials; ++i)
  {
    if (thickness[i] == 0) continue;
    cyl = new PHG4CylinderSubsystem(object->get_name(), G4TrackingService::subsysID);
    cyl->set_double_param("place_z", object->get_z_south() + length/2 + G4TrackingService::GlobalOffset);
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

double TrackingService(PHG4Reco* g4Reco, double radius)
{
  vector<ServiceProperties*> cylinders, cones;

  cylinders.push_back(new ServiceProperties("ETrackingCylinderService", 1, 0, 0.01, 0.1, 1, -80, -50, 10.75, 10.75));
  cylinders.push_back(new ServiceProperties("BTrackingCylinderService_1", 0.1, 0, 0.001, 0.01, 0.1, -50, 50, 10.75, 10.75));
  cylinders.push_back(new ServiceProperties("BTrackingCylinderService_2", 0.1, 0, 0.001, 0.01, 0.1, -40, 40, 7, 7));
  cylinders.push_back(new ServiceProperties("BTrackingCylinderService_3", 0.1, 0, 0.001, 0.01, 0.1, -30, 30, 5, 5));
  cylinders.push_back(new ServiceProperties("HTrackingCylinderService", 1, 0, 0.01, 0.1, 1, 50, 80, 10.75, 10.75));
 
  cones.push_back(new ServiceProperties("ETrackingConeService_1", 1, 0, 0.01, 0.1, 1, -50, -30, 10.75, 7));
  cones.push_back(new ServiceProperties("ETrackingConeService_2", 1, 0, 0.01, 0.1, 1, -30, -20, 7, 5));
  cones.push_back(new ServiceProperties("HTrackingConeService_1", 1, 0, 0.01, 0.1, 1, 30, 50, 7, 10.75));
  cones.push_back(new ServiceProperties("HTrackingConeService_2", 1, 0, 0.01, 0.1, 1, 20, 30, 5, 7));

  for (ServiceProperties *cylinder : cylinders) radius += TrackingServiceCylinder(cylinder, g4Reco, radius); 
  for (ServiceProperties *cone : cones) radius += TrackingServiceCone(cone, g4Reco, radius);

  return radius;
}

#endif
