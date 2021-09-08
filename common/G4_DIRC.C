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
  //dircSubsys->set_double_param("place_z", z + dz * 0.5);// relative position to mother vol.
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
