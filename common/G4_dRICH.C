/*!
 * \file G4_RICH.C
 * \brief Setup the gas RICH detector as designed in ePHENIX LOI
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.2 $
 * \date $Date: 2013/10/09 01:08:17 $
 */
#ifndef MACRO_G4dRICH_C
#define MACRO_G4dRICH_C

#include <GlobalVariables.C>

#include <g4drich/EICG4dRICHSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>

#include <g4main/PHG4Reco.h>

R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libEICG4dRICH.so)

namespace Enable
{
  bool RICH = false;
  bool RICH_OVERLAPCHECK = false;
  int RICH_VERBOSITY = 0;
}  // namespace Enable

void RICHInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 210.);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 280.);
}

//! Fast geometry model for dRICH
//Grzegorz Kalicy <gkalicy@jlab.org>
//- it starts 180 cm from IP
//- radius of aerogel part starts at 110 cm at rises up to 120cm over 20 cm of length.
//- at 175 cm from IP it rapidly grows radially to radius 210cm at stays constant like a cylinder until end at 280cm from the IP
void RICHSetup(PHG4Reco* g4Reco)
{
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::RICH_OVERLAPCHECK;
  int verbosity = std::max(Enable::VERBOSITY, Enable::RICH_VERBOSITY);

  double z = 185;
  double dz = 20;

  EICG4dRICHSubsystem *drichSubsys = new EICG4dRICHSubsystem("dRICh");
  //drichSubsys->SetGeometryFile(string(getenv("CALIBRATIONROOT")) + "/dRICH/mapping/drich-g4model.txt");
  drichSubsys->SetGeometryFile(string(getenv("ECCE")) + "/fun4all_eiccalibrations/dRICH/mapping/drich-g4model.txt");
  drichSubsys->set_double_param("place_z", z + dz * 0.5);// relative position to mother vol.
  drichSubsys->OverlapCheck(OverlapCheck);
  drichSubsys->Verbosity(verbosity);
  drichSubsys->SetActive();

  g4Reco->registerSubsystem(drichSubsys);

  if (TRACKING::FastKalmanFilter)
  {
    // project to an reference plane at z=170 cm
    TRACKING::FastKalmanFilter-> add_zplane_state("RICH", 185);
    TRACKING::ProjectionNames.insert("RICH");
  }

}
#endif
