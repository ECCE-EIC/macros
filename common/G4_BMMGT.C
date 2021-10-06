#ifndef MACRO_G4BMMGT_C
#define MACRO_G4BMMGT_C

#include <GlobalVariables.C>

#include <fun4all/Fun4AllServer.h>
#include <g4barrelmmg/PHG4CylinderStripSubsystem.h>
#include <g4barrelmmg/CreateCZHitContainer.h>
#include <g4main/PHG4Reco.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>

R__LOAD_LIBRARY(libg4barrelmmg.so)
/*
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)
*/
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4testbench.so)
R__LOAD_LIBRARY(libg4detectors.so)

namespace Enable
{
  bool BMMG = false;
  bool use_2Dreadout = true;
  bool MMG_OVERLAPCHECK = true;
  int BMMG_VERBOSITY = 1;
  
}

namespace BMMG
{
  
  const int n_layer  = 3; 
  
  const double rad[BMMG::n_layer] = {60.2, 65.4, 70.4}; // approximate radial location
  const double len[BMMG::n_layer] = {260, 270, 280.0}; 
  
  
  
}

void BMMGInit(int verbosity = 1)
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, BMMG::rad[BMMG::n_layer - 1] / 10. + 0.7);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, BMMG::len[BMMG::n_layer - 1] / 2.0+1.0);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -BMMG::len[BMMG::n_layer - 1] / 2.0-10.);
}

void BMMGSetup(PHG4Reco *g4Reco)
{
  
  Fun4AllServer* se = Fun4AllServer::instance();
  se->Verbosity(INT_MAX-10);
  
  bool OverlapCheck = Enable::OVERLAPCHECK||Enable::MMG_OVERLAPCHECK;
  
  gSystem->Load("libfun4all");
  gSystem->Load("libg4detectors.so");
  gSystem->Load("libg4testbench.so");
  gSystem->Load("libg4trackfastsim.so");

  PHG4CylinderStripSubsystem *barrel_mmg;
  double gap_betweenCZ = 1.5;
  double Gap_betweenlayer = 1.5;
  double thickness = 0.36499;
  int nCZlayer = 2;
  if (Enable::use_2Dreadout) {
    gap_betweenCZ = 0;
    nCZlayer = 1;
  }
  
  
  const double prapidity = 1;
  
  for (int ilayer = 0; ilayer< BMMG::n_layer; ilayer++){
    barrel_mmg = new PHG4CylinderStripSubsystem(Form("BMT_%d", ilayer),ilayer);
    barrel_mmg->set_double_param("radius", BMMG::rad[ilayer]);
    barrel_mmg->set_string_param("gas", "myMMGas");
    //barrel_mmg->set_double_param("steplimits", 300e-4);
    barrel_mmg->set_double_param("phi0", 15*ilayer);
    barrel_mmg->set_double_param("gap", gap_betweenCZ);
    barrel_mmg->SetActive();
    barrel_mmg->SuperDetector("BMT");
    barrel_mmg->set_int_param("lengthviarapidity",0);
    barrel_mmg->set_double_param("length", BMMG::len[ilayer]);
    barrel_mmg->set_double_param("deadzone", 0.2);
    barrel_mmg->set_int_param("nhit", 1);
    barrel_mmg->OverlapCheck(true);
    barrel_mmg->set_int_param("use_2Dreadout",Enable::use_2Dreadout);
    g4Reco->registerSubsystem(barrel_mmg);
    //barrel_mmg->Print();
    
    //
    /*
    if(TRACKING::FastKalmanFilter)
      {
	if(Enable::use_2Dreadout){
	  TRACKING::FastKalmanFilter->add_phg4hits(
						   "G4HIT_BMT",                 //      const std::string& phg4hitsNames,
						   PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
						   2.5/2/sqrt(12),//1./sqrt(12),                      //       radial-resolution , only used for Vertical Plane Detector Type
						   100e-4,//150e-4,                       //        azimuthal-resolution [cm]
						   100e-4,//150e-4,                           //      z-resolution [cm]
						   1,                           //      efficiency,
						   0);                          //      noise hits
					       
      
	}
	else {
      
	  TRACKING::FastKalmanFilter->add_phg4hits(
						   "G4HIT_CZBMT",                 //      const std::string& phg4hitsNames,
						   PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
						   2.5/2/sqrt(12),//1/sqrt(12),                      //       radial-resolution [cm], only used for Vertical Plane Detector Type
						   100e-4,//150e-4,                       //        azimuthal-resolution [cm]
						   100e-4,//150e-4,                           //      z-resolution [cm]
						   1,                           //      efficiency,
						   0);                           //      noise hits
					  
      
      
	}
    
      }

    if (TRACKING::FastKalmanFilterInnerTrack)
      {
	if(Enable::use_2Dreadout){
	  TRACKING::FastKalmanFilterInnerTrack->add_phg4hits(
							     "G4HIT_BMT",  //      const std::string& phg4hitsNames,
							     PHG4TrackFastSim::Cylinder,                         //      const DETECTOR_TYPE phg4dettype,
							     2.5/2/sqrt(12),                                     //      const float radres,
							     100e-4,                                              //      const float phires,
							     100e-4,                                              //      const float lonres,
							     1,                                                  //      const float eff,
							     0);							                        //      const float noise
	}
	else {
	  TRACKING::FastKalmanFilterInnerTrack->add_phg4hits(
							     "G4HIT_CZBMT",                //      const std::string& phg4hitsNames,
							     PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
							     2.5/2/sqrt(12),//1/sqrt(12),                      //       radial-resolution [cm], only used for Vertical Plane Detector Type
							     100e-4,//150e-4,                       //        azimuthal-resolution [cm]
							     100e-4,//150e-4,                           //      z-resolution [cm]
							     1,                           //      efficiency,
							     0);                           //      noise hits
							    
	}
       
      }

  
    if (TRACKING::FastKalmanFilterSiliconTrack and BMMG::rad[ilayer] < 50)
      {
	if(Enable::use_2Dreadout){
	  TRACKING::FastKalmanFilterSiliconTrack->add_phg4hits(
							       "G4HIT_BMT",  //      const std::string& phg4hitsNames,
							       PHG4TrackFastSim::Cylinder,                         //      const DETECTOR_TYPE phg4dettype,
							       2.5/2/sqrt(12) ,                                     //      const float radres,
							       100e-4,                                              //      const float phires,
							       100e-4,                                              //      const float lonres,
							       1,                                                  //      const float eff,
							       0);            
	}
	else {
	  TRACKING::FastKalmanFilterSiliconTrack->add_phg4hits(
							       "G4HIT_CZBMT",  //      const std::string& phg4hitsNames,
							       PHG4TrackFastSim::Cylinder,                         //      const DETECTOR_TYPE phg4dettype,
							       2.5/2/sqrt(12) ,                                     //      const float radres,
							       100e-4,                                              //      const float phires,
							       100e-4,                                              //      const float lonres,
							       1,                                                  //      const float eff,
							       0);            


	}
 
      }
*/

  }// ilayer loop
  

   if(TRACKING::FastKalmanFilter)
      {
	if(Enable::use_2Dreadout){
	  TRACKING::FastKalmanFilter->add_phg4hits(
						   "G4HIT_BMT",                 //      const std::string& phg4hitsNames,
						   PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
						   2.5/2/sqrt(12),//1./sqrt(12),                      //       radial-resolution , only used for Vertical Plane Detector Type
						   75e-4,//150e-4,                       //        azimuthal-resolution [cm]
						   75e-4,//150e-4,                           //      z-resolution [cm]
						   1,                           //      efficiency,
						   0);                          //      noise hits
					       
      
	}
	else {
      
	  TRACKING::FastKalmanFilter->add_phg4hits(
						   "G4HIT_CZBMT",                 //      const std::string& phg4hitsNames,
						   PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
						   2.5/2/sqrt(12),//1/sqrt(12),                      //       radial-resolution [cm], only used for Vertical Plane Detector Type
						   75e-4,//150e-4,                       //        azimuthal-resolution [cm]
						   75e-4,//150e-4,                           //      z-resolution [cm]
						   1,                           //      efficiency,
						   0);                           //      noise hits
					  
      
      
	}
    
      }

   

  return;
 
}

void BMMGT_Reco()
{
  
  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4detectors.so");
  
  
  int verbosity = std::max(Enable::VERBOSITY, Enable::BMMG_VERBOSITY);
  Fun4AllServer* se = Fun4AllServer::instance();
  //se->Verbosity(INT_MAX-10);
  
  return;
  
}
#endif
