#pragma once

#include "GlobalVariables.C"

#include <g4main/PHG4Reco.h>

#include <g4mrich/PHG4mRICHSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <eccefastpidreco/ECCEFastPIDReco.h>
#include <eccefastpidreco/ECCEmRICHFastPIDMap.h>
/*!
 * \file G4_mRICH.C
 * \brief Aerogel mRICH for EIC detector
 * \author Murad Sarsour
 * \date $Date: 2020/7/2  $
 */

R__LOAD_LIBRARY(libECCEFastPIDReco.so)

namespace Enable
{
  bool mRICH = false;
  bool mRICH_RECO = false;
  bool mRICH_OVERLAPCHECK = false;
  int mRICH_VERBOSITY = 0;
}  // namespace Enable

void mRICHInit()
{
}

//-1: single module
// 0: build h-side sectors and e-side wall
// 1: build h-side sectors
// 2: build e-side wall
// 3: build h-side wall
// 4: build h-side wall and e-side wall
// 5: build b-side sectors
// 6: build projected e-side wall

void mRICHSetup(PHG4Reco *g4Reco, const int detectorSetup = 1,  //1: full setup; 0:skeleton
                const int mRICHsystemSetup = 6                  //5//2//-1//
)
{
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::mRICH_OVERLAPCHECK;

  PHG4mRICHSubsystem *mRICH = new PHG4mRICHSubsystem("mRICH", 1);
  mRICH->set_int_param("detectorSetup", detectorSetup);
  mRICH->set_int_param("subsystemSetup", mRICHsystemSetup);
  mRICH->UseCalibFiles(PHG4DetectorSubsystem::xml);
  mRICH->SetCalibrationFileDir(string(getenv("CALIBRATIONROOT")) + string("/mRICH/GeometryV2/"));
  mRICH->OverlapCheck(OverlapCheck);
  //  mRICH->Verbosity(5);

  g4Reco->registerSubsystem(mRICH);

  if (TRACKING::FastKalmanFilter)
  {
    // project to an reference plane at z=170 cm
    TRACKING::FastKalmanFilter-> add_zplane_state("mRICH", -125);
    TRACKING::ProjectionNames.insert("mRICH");
  }
}


void mRICHReco(int verbosity = 0)
{
  verbosity = std::max(Enable::VERBOSITY, verbosity);
  verbosity = std::max(Enable::mRICH_VERBOSITY, verbosity);

  Fun4AllServer *se = Fun4AllServer::instance();

  ECCEmRICHFastPIDMap * pidmap = new ECCEmRICHFastPIDMap();
  pidmap->Verbosity(verbosity);

  ECCEFastPIDReco * reco = new ECCEFastPIDReco(pidmap, EICPIDDefs::mRICH, "ECCEFastPIDReco-mRICH");
  reco->Verbosity(verbosity);

  se->registerSubsystem(reco);
}

void mRICH_Eval(std::string outputfile, int verbosity = 0)
{
  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4eval.so");
  Fun4AllServer *se = Fun4AllServer::instance();

  return;
}
