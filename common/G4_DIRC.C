/*!
 * \file G4_DIRC.C
 * \brief Setup the gas DIRC detector full sim
 * \author Cameron Dean <cdean@bnl.gov>
 * \version $Revision: 0.1 $
 * \date $Date: 2021/08/17 $
 */
#ifndef MACRO_G4DIRC_C
#define MACRO_G4DIRC_C

#include <GlobalVariables.C>

#include <g4eicdirc/G4EicDircSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>

#include <g4main/PHG4Reco.h>

R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4eicdirc.so)

namespace Enable
{
  bool DIRC = false;
  bool DIRC_OVERLAPCHECK = false;
  int DIRC_VERBOSITY = 0;
  double DIRC_SCALE = 10; //DIRC class is in mm, ECCE is in cm
}  // namespace Enable

void DIRCInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 210.);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 280.);
}

void DIRCSetup(PHG4Reco* g4Reco)
{
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::DIRC_OVERLAPCHECK;
  int verbosity = std::max(Enable::VERBOSITY, Enable::DIRC_VERBOSITY);

  G4EicDircSubsystem *dircSubsys = new G4EicDircSubsystem("hpDIRC");
  dircSubsys->SuperDetector("hpDIRC");
  dircSubsys->set_double_param("rMin", 74.1 * Enable::DIRC_SCALE);
  dircSubsys->set_double_param("rMin_inner", 60.0 * Enable::DIRC_SCALE);
  dircSubsys->set_double_param("length", (287 + 168) * Enable::DIRC_SCALE);
  dircSubsys->set_double_param("NBars", 11);
  dircSubsys->set_double_param("Radius", 75.0 * Enable::DIRC_SCALE);
  dircSubsys->set_double_param("Prizm_width", 38.65 * Enable::DIRC_SCALE);
  dircSubsys->set_double_param("Prizm_length", 30.0 * Enable::DIRC_SCALE);
  dircSubsys->set_double_param("Prizm_height_at_lens", 5 * Enable::DIRC_SCALE); // override internal default
  dircSubsys->set_double_param("Bar_thickness", 1.7 * Enable::DIRC_SCALE);
  dircSubsys->set_double_param("Bar_width", 3.5 * Enable::DIRC_SCALE);
  dircSubsys->set_double_param("BarL_length", 122.5 * Enable::DIRC_SCALE);
  dircSubsys->set_double_param("BarS_length", 56.0 * Enable::DIRC_SCALE);
  dircSubsys->set_double_param("Mirror_height", 2.0 * Enable::DIRC_SCALE);
  dircSubsys->set_double_param("z_shift", -45 * Enable::DIRC_SCALE);
  dircSubsys->set_int_param("Geom_type", 0); // 0-whole DIRC, 1-one bar box
  dircSubsys->set_int_param("Lens_id", 3); // 3- 3-layer spherical lens
  dircSubsys->set_int_param("MCP_rows", 6);
  dircSubsys->set_int_param("MCP_columns", 4);
  dircSubsys->set_int_param("NBoxes",12);
  dircSubsys->set_int_param("Bar_pieces", 4);

  dircSubsys->OverlapCheck(OverlapCheck);
  dircSubsys->Verbosity(verbosity);
  dircSubsys->SetActive();

  g4Reco->registerSubsystem(dircSubsys);

  if (TRACKING::FastKalmanFilter)
  {
    // project to an reference plane at z=170 cm
    TRACKING::FastKalmanFilter-> add_zplane_state("DIRC", 185);
    TRACKING::ProjectionNames.insert("DIRC");
  }

}
#endif
