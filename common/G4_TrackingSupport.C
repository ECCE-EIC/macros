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
                               const int &n_staves,
                               const double &z_south,
                               const double &z_north,
                               const double &r_south,
                               const double &r_north);

    virtual ~ServiceProperties(){};

    const string get_name();
    const int get_n_staves();
    const double get_z_south();
    const double get_z_north();
    const double get_r_south();
    const double get_r_north();
  
  private:
    const string m_name = "service";
    const int m_n_staves = 48;
    const double m_z_south = 0.0;
    const double m_z_north = 0.0;
    const double m_r_south = 0.0;
    const double m_r_north = 0.0;
};

ServiceProperties::ServiceProperties(const string &name,
                                     const int &n_staves,
                                     const double &z_south,
                                     const double &z_north,
                                     const double &r_south,
                                     const double &r_north)
  : m_name(name)
  , m_n_staves(n_staves)
  , m_z_south(z_south)
  , m_z_north(z_north)
  , m_r_south(r_south)
  , m_r_north(r_north)
{}

const string ServiceProperties::get_name() { return m_name; }
const int ServiceProperties::get_n_staves() { return m_n_staves; }
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
{
  std::string materials[] = {"G4_Cu", "G4_WATER", "G4_POLYETHYLENE", "PEEK"};

  double single_stave_service_copper_area = 0.0677;   //Cross-sectional area of copper for 1 stave [cm^2]
  double single_stave_service_water_area = 0.0098;    //Cross-sectional area of water for 1 stave [cm^2]
  double single_stave_service_plastic_area = 0.4303;  //Cross-sectional area of plastic for 1 stave [cm^2]
}  // namespace G4TrackingService

void TrackingServiceInit()
{
}

double calculateOR(double inner_radius, double area)  //Calculate the outer radius of a disk, knowing the inner radius and the area
{
  return std::sqrt(area / M_PI + std::pow(inner_radius, 2));
}

void calculateMaterialBoundaries(ServiceProperties *properties, double& outer_copper_radius, double& outer_water_radius, double& outer_plastic_radius, bool isSouth)  //Calculate where the transition between each material occurs
{
  double start_radius = isSouth ? properties->get_r_south() : properties->get_r_north();

  outer_copper_radius = calculateOR(start_radius, properties->get_n_staves() * G4TrackingService::single_stave_service_copper_area);
  outer_water_radius = calculateOR(outer_copper_radius, properties->get_n_staves() * G4TrackingService::single_stave_service_water_area);
  outer_plastic_radius = calculateOR(outer_water_radius, properties->get_n_staves() * G4TrackingService::single_stave_service_plastic_area);
}

double TrackingServiceCone(ServiceProperties *object, PHG4Reco* g4Reco, double radius)
{
  bool AbsorberActive = Enable::ABSORBER || Enable::TrackingService_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::TrackingService_OVERLAPCHECK;
  int verbosity = std::max(Enable::VERBOSITY, Enable::TrackingService_VERBOSITY);

  PHG4ConeSubsystem* cone;

  double materialBoundariesSouth[5] = {object->get_r_south(), 0, 0, 0, 0.1};
  calculateMaterialBoundaries(object, materialBoundariesSouth[1], materialBoundariesSouth[2], materialBoundariesSouth[3], true);
  materialBoundariesSouth[4] += materialBoundariesSouth[3];

  double materialBoundariesNorth[5] = {object->get_r_north(), 0, 0, 0, 0.1};
  calculateMaterialBoundaries(object, materialBoundariesNorth[1], materialBoundariesNorth[2], materialBoundariesNorth[3], false);
  materialBoundariesNorth[4] += materialBoundariesNorth[3];

  for (int i = 0; i < 4; ++i)
  {
    cone = new PHG4ConeSubsystem(object->get_name(), i);
    cone->SetR1(materialBoundariesSouth[i], materialBoundariesSouth[i+1]);
    cone->SetR2(materialBoundariesNorth[i], materialBoundariesNorth[i+1]);
    cone->SetPlaceZ(object->get_z_south());
    cone->SetZlength(abs(object->get_z_north() - object->get_z_south()));
    cone->SetMaterial(G4TrackingService::materials[i]);
    cone->SuperDetector("TrackingService");
    if (AbsorberActive) cone->SetActive();
    cone->OverlapCheck(OverlapCheck);
    g4Reco->registerSubsystem(cone);
  }
  radius = max(materialBoundariesSouth[4], materialBoundariesNorth[4]);

  return radius;
}

double TrackingServiceCylinder(ServiceProperties *object, PHG4Reco* g4Reco, double radius)
{
  bool AbsorberActive = Enable::ABSORBER || Enable::TrackingService_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::TrackingService_OVERLAPCHECK;
  int verbosity = std::max(Enable::VERBOSITY, Enable::TrackingService_VERBOSITY);

  PHG4CylinderSubsystem* cyl;

  double materialBoundaries[5] = {object->get_r_south(), 0, 0, 0, 0.1};
  calculateMaterialBoundaries(object, materialBoundaries[1], materialBoundaries[2], materialBoundaries[3], true);
  materialBoundaries[4] += materialBoundaries[3];

  for (int i = 0; i < 4; ++i)
  {
     cyl = new PHG4CylinderSubsystem(object->get_name(), i);
     cyl->set_double_param("place_z", object->get_z_south());
     cyl->set_double_param("radius", materialBoundaries[i]);
     cyl->set_double_param("length", abs(object->get_z_north() - object->get_z_south()));
     cyl->set_string_param("material", G4TrackingService::materials[i]);
     cyl->set_double_param("thickness", materialBoundaries[i + 1] - materialBoundaries[i]);
     cyl->SuperDetector("TrackingService");
     if (AbsorberActive) cyl->SetActive();
     cyl->OverlapCheck(OverlapCheck);
     g4Reco->registerSubsystem(cyl);
  }

  radius = materialBoundaries[4];

  return radius;
}

double TrackingService(PHG4Reco* g4Reco, double radius)
{
  std::vector<ServiceProperties*> cylinders, cones;

  cylinders.push_back(new ServiceProperties("ETrackingCylinderService", 48, -200, -50, 10.75, 10.75));

  cones.push_back(new ServiceProperties("ETrackingConeService", 48, -50, -20, 10.75, 5));

  for (ServiceProperties *cylinder : cylinders)
  {
    radius += TrackingServiceCylinder(cylinder, g4Reco, radius);
  }
  
  for (ServiceProperties *cone : cones)
  {
    radius += TrackingServiceCone(cone, g4Reco, radius);
  }

  return radius;
}

#endif
