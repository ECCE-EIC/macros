#ifndef MACRO_G4EToF_C
#define MACRO_G4EToF_C

#include <GlobalVariables.C>

#include <fun4all/Fun4AllServer.h>
#include <g4etof/PHG4ECAPToFSubsystem.h>

#include <g4main/PHG4Reco.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>

R__LOAD_LIBRARY(libg4etof.so)
//R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)

namespace Enable
{
  bool ETOF = false;
  bool ETOF_GAS = false;
  bool ETOF_OVERLAPCHECK = false;
  int ETOF_VERBOSITY = 0;
}  // namespace Enable

namespace ETOF
{
  int f_gas_lyr = 6.;   // total number of layers
  int f_mrpc_lyr = 7.;  //total number of layers
  int b_gas_lyr = 6.;   // total number of layers
  int b_mrpc_lyr = 7.;  //total number of layers
  int pcb_lyr = 3.;
  int mylar_lyr = 4.;
  int cu_lyr = 4.;
  int carbon_lyr = 4.;
  int honeycomb_lyr = 2.;
  double z_start = -165.;         // cm , from left side of -ve Z axis moving towards origin
  double R_in = 6.5;              // cm
  double R_out = 68.;             //cm
  double gas_gap = 0.0220;        // 220 microns
  double mrpc_thick = 0.04;       // 400 microns
  double pcb_thick = 0.06;        // 600 microns
  double cu_thick = 0.003;        // 30 microns, layer over pcb, (1 each on outer pcb and 2 on central pcb = 4)
  double carbon_thick = 0.01;     // 100 microns , 2 layers
  double mylar_thick = 0.04;      // 400 microns, 4 layers
  double honeycomb_thick = 0.75;  // 7.5 mm, 2 honeycomb

  double tof_width = (f_gas_lyr + b_gas_lyr) * gas_gap + (f_mrpc_lyr + b_mrpc_lyr) * mrpc_thick + carbon_lyr * carbon_thick + pcb_lyr * pcb_thick + cu_lyr * cu_thick + mylar_lyr * mylar_thick + honeycomb_lyr * honeycomb_thick;

  double z_end = (z_start + tof_width);

}  // namespace ETOF

void ETOFInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, ETOF::R_out);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, ETOF::z_end - 5);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, ETOF::z_start + 5);
}

void ETOFSetup(PHG4Reco* g4Reco)
{
  bool GasActive = Enable::ABSORBER || Enable::ETOF_GAS;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::ETOF_OVERLAPCHECK;
  int verbosity = std::max(Enable::VERBOSITY, Enable::ETOF_VERBOSITY);

  Fun4AllServer* se = Fun4AllServer::instance();
  se->Verbosity(verbosity);

  PHG4ECAPToFSubsystem* eTOF = new PHG4ECAPToFSubsystem("eTOF", 1);
  eTOF->Verbosity(verbosity);

  eTOF->set_int_param("n_fgas_layer", ETOF::f_gas_lyr);
  eTOF->set_int_param("n_bgas_layer", ETOF::b_gas_lyr);
  eTOF->set_double_param("gas_gap", ETOF::gas_gap);
  eTOF->set_double_param("glass_thick", ETOF::mrpc_thick);
  eTOF->set_double_param("Carbon_thick", ETOF::carbon_thick);
  eTOF->set_double_param("pcb_thick", ETOF::pcb_thick);
  eTOF->set_double_param("cu_thick", ETOF::cu_thick);
  eTOF->set_double_param("honeycomb_thick", ETOF::honeycomb_thick);
  eTOF->set_double_param("mylar_thick", ETOF::mylar_thick);
  eTOF->set_double_param("Rin", ETOF::R_in);
  eTOF->set_double_param("Rout", ETOF::R_out);
  eTOF->set_double_param("z_begin", ETOF::z_start);
  eTOF->set_int_param("use_g4steps", 1);
  eTOF->SetActive(1);
  eTOF->SuperDetector("ETOF");
  if (GasActive)
  {
    eTOF->SetAbsorberActive(1);
  }
  eTOF->OverlapCheck(OverlapCheck);

  g4Reco->registerSubsystem(eTOF);

  //trd_hcap->OverlapCheck(1);

  if (verbosity > 1) cout << " ETOF gas layer  :" << ETOF::f_gas_lyr << endl;

  if (TRACKING::FastKalmanFilter)
  {
    /*
       TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + string(Form("ACTIVEGAS_ETOF")),  //      const std::string& phg4hitsNames,
                                               PHG4TrackFastSim::Vertical_Plane,                         //      const DETECTOR_TYPE phg4dettype,
                                               1, //1. / sqrt(12.),                                     //      const float radres,
                                               5.0e-1,//55e-4,                                              //      const float phires,
                                               5.0e-1,//55e-4,                                              //      const float lonres,
                                               1,                                                  //      const float eff,
                                               0);                                                 //      const float noise

       */

    //Reference plane projection at initial R of ToF
    TRACKING::FastKalmanFilter->add_zplane_state(string("G4HIT_") + string(Form("ACTIVEGAS_ETOF")), ETOF::z_start);
    TRACKING::ProjectionNames.insert(string("G4HIT_") + string(Form("ACTIVEGAS_ETOF")));
  }
  return;
}

void ETOF_Reco()
{
  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4detectors.so");

  int verbosity = std::max(Enable::VERBOSITY, Enable::ETOF_VERBOSITY);
  Fun4AllServer* se = Fun4AllServer::instance();
  //se->Verbosity(INT_MAX-10);

  return;
}

#endif
