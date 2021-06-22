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

#include <g4detectors/PHG4ConeSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>

#include <g4main/PHG4Reco.h>

R__LOAD_LIBRARY(libg4detectors.so)

namespace Enable
{
  bool RICH = false;
  bool RICH_OVERLAPCHECK = false;
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

  double z = 185;
  double dz = 0;

  dz = 20;
  PHG4ConeSubsystem* coneAeroGel = new PHG4ConeSubsystem("dRICHconeAeroGel", 0);
  coneAeroGel->SetR1(11, 110);
  coneAeroGel->SetR2(11, 120);
  coneAeroGel->SetZlength(dz * 0.5);
  coneAeroGel->SetPlaceZ(z + dz * 0.5);
  coneAeroGel->set_string_param("material", "ePHENIX_AeroGel");
  coneAeroGel->set_color(0.5, 0.5, 0.8);
  coneAeroGel->SetActive();
  coneAeroGel->SuperDetector("RICH");
  coneAeroGel->OverlapCheck(OverlapCheck);
  g4Reco->registerSubsystem(coneAeroGel);

  z += dz;
  dz = 5;
  PHG4ConeSubsystem* coneGas = new PHG4ConeSubsystem("dRICHconeGas", 1);
  coneGas->SetR1(11, 120);
  coneGas->SetR2(11, 170); // <- reduce from 210cm to avoid overlap with the HCal
  coneGas->SetZlength(dz * 0.5);
  coneGas->SetPlaceZ(z + dz * 0.5);
  coneGas->set_string_param("material", "C4F10");
  coneGas->set_color(0.5, 0.5, 0.5);
  coneGas->SetActive();
  coneGas->SuperDetector("RICH");
  coneGas->OverlapCheck(OverlapCheck);
  g4Reco->registerSubsystem(coneGas);

  z += dz;
  dz = 75;
  coneGas = new PHG4ConeSubsystem("dRICHconeGas", 2);
  coneGas->SetR1(11, 170);
  coneGas->SetR2(15, 170);
  coneGas->SetZlength(dz * 0.5);
  coneGas->SetPlaceZ(z + dz * 0.5);
  coneGas->set_string_param("material", "C4F10");
  coneGas->set_color(0.5, 0.5, 0.5);
  coneGas->SetActive();
  coneGas->SuperDetector("RICH");
  coneGas->OverlapCheck(OverlapCheck);
  g4Reco->registerSubsystem(coneGas);

  z += dz;


  if (TRACKING::FastKalmanFilter)
  {
    // project to an reference plane at z=170 cm
    TRACKING::FastKalmanFilter-> add_zplane_state("RICH", 185);
    TRACKING::ProjectionNames.insert("RICH");
  }

}
#endif
