#ifndef MACRO_G4mRWELL_C
#define MACRO_G4mRWELL_C

#include <GlobalVariables.C>

#include <fun4all/Fun4AllServer.h>
#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4detectors/PHG4SectorSubsystem.h>
#include <g4main/PHG4Reco.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4testbench.so)
R__LOAD_LIBRARY(libg4detectors.so)

namespace Enable
{
  bool RWELL = false;
  bool RWELL_OVERLAPCHECK = false;
}  // namespace Enable

namespace RWELL
{
  //All units specified in cm unless stated otherwise
  //  const int n_layer = 2;  //tracker layers
  //  //const double nom_radius[RWELL::n_layer] = {79.5,90.0}; //77 to not hit DIRC
  //  //  const double nom_radius[RWELL::n_layer] = {78.67, 90.0};  //77 to not hit DIRC
  //  const double nom_radius[RWELL::n_layer] = {69 - 1.6 - 2.6, 78.67};  //77 to not hit DIRC
  //  const double nom_driftgap[RWELL::n_layer] = {0.4, 0.4};
  //  const double nom_length[RWELL::n_layer] = {300.0, 300.0};

  const int n_layer = 1;  //tracker layers
  const double nom_radius[RWELL::n_layer] = {44, 47, 69 - 1.6 - 2.6};
  const double nom_driftgap[RWELL::n_layer] = {0.4, 0.4, 0.4};
  const double nom_length[RWELL::n_layer] = {150, 150, 300.0};
}  //namespace RWELL

void RWellInit(int verbosity = 0)
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, RWELL::nom_radius[RWELL::n_layer - 1] / 10. + 0.7);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, RWELL::nom_length[RWELL::n_layer - 1] / 2.0);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -RWELL::nom_length[RWELL::n_layer - 1] / 2.0);
}

double Build_G4_RWell_Bare(PHG4Reco* g4Reco,
                           double rwellrad = 80.0,
                           double driftgap = 1.5,
                           double length = 200.0,
                           int index = 0)
{
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::RWELL_OVERLAPCHECK;

  gSystem->Load("libfun4all");
  gSystem->Load("libg4detectors.so");
  gSystem->Load("libg4testbench.so");
  gSystem->Load("libg4trackfastsim.so");

  //double driftgap = 3.0;
  //double rwellrad = radius;
  //double length = length;
  double rsum = 0.0;

  //MPGD parameters (units in cm)
  double kapton_thickness = 0.0175;  //cm
  double cu_thickness = 0.002;       //cm
  double pcb_thickness = 0.010;      //cm
  double prepreg_thickness = 0.005;  //cm

  //Suppoort parameters (units in cm)
  //string supMat = "PEEK"; //support material
  //string supMat = "G4_Galactic"; //support material
  // making carbon fiber epoxy
  // G4Material *cfrp_intt = new G4Material("CFRP_INTT", density = 1.69 * g / cm3, ncomponents = 3);
  string supMat = "CFRP_INTT";
  //inner tube
  double support_01_thickness = 0.50;
  double support_01_length = 7.2;
  //inner ring
  double support_02_thickness = 1.6;
  double support_02_length = 1.2;
  //outer ring
  double support_03_thickness = 0.50;
  double support_03_length = 1.2;

  PHG4CylinderSubsystem* rwell_cyl(nullptr);

  // here is our uRwell:
  //Gass layer
  rwell_cyl = new PHG4CylinderSubsystem(Form("RWELL_%d", index), 0);
  rwell_cyl->set_double_param("radius", rwellrad);
  rwell_cyl->set_string_param("material", "G4_Ar");
  rwell_cyl->set_double_param("thickness", driftgap);
  rwell_cyl->SetActive(1);
  rwell_cyl->set_int_param("lengthviarapidity", 0);
  rwell_cyl->set_double_param("length", length);
  rwell_cyl->SuperDetector(Form("RWELL_%d", index));
  rwell_cyl->SetActive(1);
  g4Reco->registerSubsystem(rwell_cyl);
  rwell_cyl->OverlapCheck(OverlapCheck);

  //Kapton
  rwell_cyl = new PHG4CylinderSubsystem(Form("RWELL_Kapton_%d", index), 0);
  rwell_cyl->set_double_param("radius", rwellrad - kapton_thickness);
  rwell_cyl->set_string_param("material", "G4_KAPTON");
  rwell_cyl->set_double_param("thickness", kapton_thickness);
  rwell_cyl->set_int_param("lengthviarapidity", 0);
  rwell_cyl->set_double_param("length", length);
  rwell_cyl->SuperDetector(Form("RWELL_Kapton_%d", index));
  rwell_cyl->SetActive(0);
  g4Reco->registerSubsystem(rwell_cyl);
  rwell_cyl->OverlapCheck(OverlapCheck);
  //Cu
  rsum = rwellrad + driftgap;
  rwell_cyl = new PHG4CylinderSubsystem(Form("RWELL_Cu_%d", index), 0);
  rwell_cyl->set_double_param("radius", rsum);
  rwell_cyl->set_string_param("material", "G4_Cu");
  rwell_cyl->set_double_param("thickness", cu_thickness);
  rwell_cyl->set_int_param("lengthviarapidity", 0);
  rwell_cyl->set_double_param("length", length);
  rwell_cyl->SuperDetector(Form("RWELL_Cu_%d", index));
  rwell_cyl->SetActive(0);
  g4Reco->registerSubsystem(rwell_cyl);
  rwell_cyl->OverlapCheck(OverlapCheck);

  //Prepreg
  rsum += cu_thickness;
  rwell_cyl = new PHG4CylinderSubsystem(Form("RWELL_PrePreg_%d", index), 0);
  rwell_cyl->set_double_param("radius", rsum);
  rwell_cyl->set_string_param("material", "NOMEX");
  rwell_cyl->set_double_param("thickness", prepreg_thickness);
  rwell_cyl->set_int_param("lengthviarapidity", 0);
  rwell_cyl->set_double_param("length", length);
  rwell_cyl->SuperDetector(Form("RWELL_PrePreg_%d", index));
  rwell_cyl->SetActive(0);
  g4Reco->registerSubsystem(rwell_cyl);
  rwell_cyl->OverlapCheck(OverlapCheck);

  //PCB
  rsum += prepreg_thickness;
  rwell_cyl = new PHG4CylinderSubsystem(Form("RWELL_PCB_%d", index), 0);
  rwell_cyl->set_double_param("radius", rsum);
  rwell_cyl->set_string_param("material", "FR4");
  rwell_cyl->set_double_param("thickness", pcb_thickness);
  rwell_cyl->set_int_param("lengthviarapidity", 0);
  rwell_cyl->set_double_param("length", length);
  rwell_cyl->SuperDetector(Form("RWELL_PCB_%d", index));
  rwell_cyl->SetActive(0);
  g4Reco->registerSubsystem(rwell_cyl);
  rwell_cyl->OverlapCheck(OverlapCheck);

  return rwellrad;
}

double Build_G4_RWell_Sup01(PHG4Reco* g4Reco,
                            double rwellrad = 80.0,
                            double driftgap = 1.5,
                            double length = 200.0,
                            int index = 0)
{
  gSystem->Load("libfun4all");
  gSystem->Load("libg4detectors.so");
  gSystem->Load("libg4testbench.so");

  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::RWELL_OVERLAPCHECK;

  //double driftgap = gap;
  //double rwellrad = radius;
  //double length = length;
  double rsum = 0.0;

  //MPGD parameters (units in cm)
  double kapton_thickness = 0.0175;  //cm
  double cu_thickness = 0.002;       //cm
  double pcb_thickness = 0.010;      //cm
  double prepreg_thickness = 0.005;  //cm

  //Suppoort parameters (units in cm)
  //string supMat = "PEEK"; //support material
  //string supMat = "G4_Galactic"; //support material
  // making carbon fiber epoxy
  // G4Material *cfrp_intt = new G4Material("CFRP_INTT", density = 1.69 * g / cm3, ncomponents = 3);
  string supMat = "CFRP_INTT";
  //inner tube
  double support_01_thickness = 0.50;
  double support_01_length = 7.2;
  //inner ring
  double support_02_thickness = 1.6;
  double support_02_length = 1.2;
  //outer ring
  double support_03_thickness = 0.50;
  double support_03_length = 1.2;

  PHG4CylinderSubsystem* rwell_cyl(nullptr);

  // here is our uRwell:
  //Gass layer
  rwell_cyl = new PHG4CylinderSubsystem(Form("RWELL_%d", index), 0);
  rwell_cyl->set_double_param("radius", rwellrad);
  rwell_cyl->set_string_param("material", "G4_METHANE");
  rwell_cyl->set_double_param("thickness", driftgap);
  rwell_cyl->SetActive(1);
  rwell_cyl->set_int_param("lengthviarapidity", 0);
  rwell_cyl->set_double_param("length", length);
  rwell_cyl->SuperDetector(Form("RWELL_%d", index));
  rwell_cyl->SetActive(1);
  g4Reco->registerSubsystem(rwell_cyl);
  rwell_cyl->OverlapCheck(OverlapCheck);

  //Kapton
  rwell_cyl = new PHG4CylinderSubsystem(Form("RWELL_Kapton_%d", index), 0);
  rwell_cyl->set_double_param("radius", rwellrad - kapton_thickness);
  rwell_cyl->set_string_param("material", "G4_KAPTON");
  rwell_cyl->set_double_param("thickness", kapton_thickness);
  rwell_cyl->set_int_param("lengthviarapidity", 0);
  rwell_cyl->set_double_param("length", length);
  rwell_cyl->SuperDetector(Form("RWELL_Kapton_%d", index));
  rwell_cyl->SetActive(0);
  g4Reco->registerSubsystem(rwell_cyl);
  rwell_cyl->OverlapCheck(OverlapCheck);
  //Cu
  rsum = rwellrad + driftgap;
  rwell_cyl = new PHG4CylinderSubsystem(Form("RWELL_Cu_%d", index), 0);
  rwell_cyl->set_double_param("radius", rsum);
  rwell_cyl->set_string_param("material", "G4_Cu");
  rwell_cyl->set_double_param("thickness", cu_thickness);
  rwell_cyl->set_int_param("lengthviarapidity", 0);
  rwell_cyl->set_double_param("length", length);
  rwell_cyl->SuperDetector(Form("RWELL_Cu_%d", index));
  rwell_cyl->SetActive(0);
  g4Reco->registerSubsystem(rwell_cyl);
  rwell_cyl->OverlapCheck(OverlapCheck);

  //Prepreg
  rsum += cu_thickness;
  rwell_cyl = new PHG4CylinderSubsystem(Form("RWELL_PrePreg_%d", index), 0);
  rwell_cyl->set_double_param("radius", rsum);
  rwell_cyl->set_string_param("material", "NOMEX");
  rwell_cyl->set_double_param("thickness", prepreg_thickness);
  rwell_cyl->set_int_param("lengthviarapidity", 0);
  rwell_cyl->set_double_param("length", length);
  rwell_cyl->SuperDetector(Form("RWELL_PrePreg_%d", index));
  rwell_cyl->SetActive(0);
  g4Reco->registerSubsystem(rwell_cyl);
  rwell_cyl->OverlapCheck(OverlapCheck);

  //PCB
  rsum += prepreg_thickness;
  rwell_cyl = new PHG4CylinderSubsystem(Form("RWELL_PCB_%d", index), 0);
  rwell_cyl->set_double_param("radius", rsum);
  rwell_cyl->set_string_param("material", "FR4");
  rwell_cyl->set_double_param("thickness", pcb_thickness);
  rwell_cyl->set_int_param("lengthviarapidity", 0);
  rwell_cyl->set_double_param("length", length);
  rwell_cyl->SuperDetector(Form("RWELL_PCB_%d", index));
  rwell_cyl->SetActive(0);
  g4Reco->registerSubsystem(rwell_cyl);
  rwell_cyl->OverlapCheck(OverlapCheck);

  //---Support structure--
  //tube
  rsum += pcb_thickness;
  rwell_cyl = new PHG4CylinderSubsystem(Form("RWELL_Support01_0_%d", index), 0);  //RWELL_<support type>_<location>_<index>
                                                                                  //Support 01 = tube
                                                                                  //Support 02 = inner ring
                                                                                  //Support 03 = outer ring
                                                                                  //location 0 = z < 0
                                                                                  //location 1 = z > 0
  rwell_cyl->set_double_param("radius", rsum);
  rwell_cyl->set_string_param("material", supMat);
  rwell_cyl->set_double_param("thickness", support_01_thickness);
  rwell_cyl->set_int_param("lengthviarapidity", 0);
  rwell_cyl->set_double_param("place_z", -length / 2);
  rwell_cyl->set_double_param("length", support_01_length);
  rwell_cyl->SuperDetector(Form("RWELL_Support01_0_%d", index));
  rwell_cyl->SetActive(0);
  g4Reco->registerSubsystem(rwell_cyl);
  rwell_cyl->OverlapCheck(OverlapCheck);
  //tube 2
  rwell_cyl = new PHG4CylinderSubsystem(Form("RWELL_Support01_1_%d", index), 0);
  rwell_cyl->set_double_param("radius", rsum);
  rwell_cyl->set_string_param("material", supMat);
  rwell_cyl->set_double_param("thickness", support_01_thickness);
  rwell_cyl->set_int_param("lengthviarapidity", 0);
  rwell_cyl->set_double_param("place_z", length / 2);
  rwell_cyl->set_double_param("length", support_01_length);
  rwell_cyl->SuperDetector(Form("RWELL_Support01_1_%d", index));
  rwell_cyl->SetActive(0);
  g4Reco->registerSubsystem(rwell_cyl);
  rwell_cyl->OverlapCheck(OverlapCheck);
  //inner ring
  rsum += support_01_thickness;
  rwell_cyl = new PHG4CylinderSubsystem(Form("RWELL_Support02_0_%d", index), 0);
  rwell_cyl->set_double_param("radius", rsum);
  rwell_cyl->set_string_param("material", supMat);
  rwell_cyl->set_double_param("thickness", support_02_thickness);
  rwell_cyl->set_int_param("lengthviarapidity", 0);
  rwell_cyl->set_double_param("place_z", -length / 2 + support_01_length / 2);
  rwell_cyl->set_double_param("length", support_02_length);
  rwell_cyl->SuperDetector(Form("RWELL_Support02_0_%d", index));
  rwell_cyl->SetActive(0);
  g4Reco->registerSubsystem(rwell_cyl);
  rwell_cyl->OverlapCheck(OverlapCheck);
  //inner ring 2
  rwell_cyl = new PHG4CylinderSubsystem(Form("RWELL_Support02_1_%d", index), 0);
  rwell_cyl->set_double_param("radius", rsum);
  rwell_cyl->set_string_param("material", supMat);
  rwell_cyl->set_double_param("thickness", support_02_thickness);
  rwell_cyl->set_int_param("lengthviarapidity", 0);
  rwell_cyl->set_double_param("place_z", length / 2 - support_01_length / 2);
  rwell_cyl->set_double_param("length", support_02_length);
  rwell_cyl->SuperDetector(Form("RWELL_Support02_1_%d", index));
  rwell_cyl->SetActive(0);
  g4Reco->registerSubsystem(rwell_cyl);
  rwell_cyl->OverlapCheck(OverlapCheck);
  //outer ring
  rsum += support_02_thickness;
  rwell_cyl = new PHG4CylinderSubsystem(Form("RWELL_Support03_0_%d", index), 0);
  rwell_cyl->set_double_param("radius", rsum);
  rwell_cyl->set_string_param("material", supMat);
  rwell_cyl->set_double_param("thickness", support_03_thickness);
  rwell_cyl->set_int_param("lengthviarapidity", 0);
  rwell_cyl->set_double_param("length", support_03_length);
  rwell_cyl->SuperDetector(Form("RWELL_Support03_0_%d", index));
  rwell_cyl->set_double_param("place_z", -length / 2 + support_01_length / 2);
  rwell_cyl->SetActive(0);
  g4Reco->registerSubsystem(rwell_cyl);
  rwell_cyl->OverlapCheck(OverlapCheck);
  //outer ring 2
  rwell_cyl = new PHG4CylinderSubsystem(Form("RWELL_Support03_1_%d", index), 0);
  rwell_cyl->set_double_param("radius", rsum);
  rwell_cyl->set_string_param("material", supMat);
  rwell_cyl->set_double_param("thickness", support_03_thickness);
  rwell_cyl->set_int_param("lengthviarapidity", 0);
  rwell_cyl->set_double_param("length", support_03_length);
  rwell_cyl->SuperDetector(Form("RWELL_Support03_1_%d", index));
  rwell_cyl->set_double_param("place_z", length / 2 - support_01_length / 2);
  rwell_cyl->SetActive(0);
  g4Reco->registerSubsystem(rwell_cyl);
  rwell_cyl->OverlapCheck(OverlapCheck);

  return rwellrad;
}

//! type selects the RWell material
//   0: bare RWell (no support structure)
//   1: Implimentation of FIT support rings
double RWellSetup(PHG4Reco* g4Reco,
                  int type = 1)
{
  double radius = 0;

  for (int ilyr = 0; ilyr < RWELL::n_layer; ilyr++)  //RWELL trackers are registered in Build_RWELL macro
  {
    if (type == 0)
    {
      radius = Build_G4_RWell_Bare(g4Reco,                     //returns RWELL radiaus
                                   RWELL::nom_radius[ilyr],    //radius
                                   RWELL::nom_driftgap[ilyr],  //driftgap,
                                   RWELL::nom_length[ilyr],    //length
                                   ilyr);                      //index
    }
    if (type == 1)
    {
      radius = Build_G4_RWell_Sup01(g4Reco,                     //returns RWELL radiaus
                                    RWELL::nom_radius[ilyr],    //radius
                                    RWELL::nom_driftgap[ilyr],  //driftgap,
                                    RWELL::nom_length[ilyr],    //length
                                    ilyr);                      //index
    }

    // sourav: For spatial resolution the mRwell I will use about 55 microns
    // (usually it is in between 40-60 microns depending on the angle of incidence of
    // primary tracks when mRwell are used in microTPC mode i.e drift gap of 3-4 mm) .
    if (TRACKING::FastKalmanFilter)
      TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + string(Form("RWELL_%d", ilyr)),  //      const std::string& phg4hitsNames,
                                               PHG4TrackFastSim::Cylinder,                         //      const DETECTOR_TYPE phg4dettype,
                                               1. / sqrt(12.),                                     //      const float radres,
                                               55e-4,                                              //      const float phires,
                                               55e-4,                                              //      const float lonres,
                                               1,                                                  //      const float eff,
                                               0);                                                 //      const float noise
  }
  return radius;  //cm
}

// Central detector cell reco is disabled as EIC setup use the fast tracking sim for now
void RWell_Cells(int verbosity = 0)
{
  // runs the cellularization of the energy deposits (g4hits)
  // into detector hits (g4cells)

  //---------------
  // Load libraries
  //---------------

  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4detectors.so");

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer* se = Fun4AllServer::instance();

  //-----------
  // SVTX cells
  //-----------

  return;
}
// Central detector  reco is disabled as EIC setup use the fast tracking sim for now
void RWell_Reco(int verbosity = 0)
{
  //---------------
  // Load libraries
  //---------------
  gSystem->Load("libfun4all.so");

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer* se = Fun4AllServer::instance();

  return;
}

#endif
