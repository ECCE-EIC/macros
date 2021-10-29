#ifndef MACRO_G4HToF_C
#define MACRO_G4HToF_C

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
  bool HTOF = false;
  bool HTOF_GAS = false;
  bool HTOF_OVERLAPCHECK = true;
  int HTOF_VERBOSITY = 0;
}  // namespace Enable

namespace HTOF
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
  double z_start = 287.;          //cm, starting point from left on +ve Z axis moving away from origin
  double R_in = 15.0;             // cm
  double R_out = 170.;            //cm
  double gas_gap = 0.0220;        // 220 microns
  double mrpc_thick = 0.04;       // 400 microns
  double pcb_thick = 0.06;        // 600 microns
  double cu_thick = 0.003;        // 30 microns, layer over pcb, (1 each on outer pcb and 2 on central pcb = 4)
  double carbon_thick = 0.01;     // 100 microns , 2 layers
  double mylar_thick = 0.04;      // 400 microns, 4 layers
  double honeycomb_thick = 0.75;  // 7.5 mm, 2 honeycomb

  double tof_width = (f_gas_lyr + b_gas_lyr) * gas_gap + (f_mrpc_lyr + b_mrpc_lyr) * mrpc_thick + carbon_lyr * carbon_thick + pcb_lyr * pcb_thick + cu_lyr * cu_thick + mylar_lyr * mylar_thick + honeycomb_lyr * honeycomb_thick;

  double z_end = (z_start + tof_width);

}  // namespace HTOF

void HTOFInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, HTOF::R_out);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, HTOF::z_end + 5);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, HTOF::z_start - 5);
}

void HTOFSetup(PHG4Reco* g4Reco)
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::HTOF_VERBOSITY);
  bool GasActive = Enable::ABSORBER || Enable::HTOF_GAS;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::HTOF_OVERLAPCHECK;

  Fun4AllServer* se = Fun4AllServer::instance();
  se->Verbosity(verbosity);

  PHG4ECAPToFSubsystem* hTOF = new PHG4ECAPToFSubsystem("hTOF", 1);
  hTOF->Verbosity(verbosity);
  hTOF->set_int_param("n_fgas_layer", HTOF::f_gas_lyr);
  hTOF->set_int_param("n_bgas_layer", HTOF::b_gas_lyr);
  hTOF->set_double_param("gas_gap", HTOF::gas_gap);
  hTOF->set_double_param("glass_thick", HTOF::mrpc_thick);
  hTOF->set_double_param("Carbon_thick", HTOF::carbon_thick);
  hTOF->set_double_param("pcb_thick", HTOF::pcb_thick);
  hTOF->set_double_param("cu_thick", HTOF::cu_thick);
  hTOF->set_double_param("honeycomb_thick", HTOF::honeycomb_thick);
  hTOF->set_double_param("mylar_thick", HTOF::mylar_thick);
  hTOF->set_double_param("Rin", HTOF::R_in);
  hTOF->set_double_param("Rout", HTOF::R_out);
  hTOF->set_double_param("z_begin", HTOF::z_start);
  hTOF->set_int_param("use_g4steps", 1);
  hTOF->SetActive(1);
  hTOF->SuperDetector("HTOF");
  if (GasActive)
  {
    hTOF->SetAbsorberActive(1);
  }
  hTOF->OverlapCheck(OverlapCheck);

  g4Reco->registerSubsystem(hTOF);

  //trd_hcap->OverlapCheck(1);

  if (verbosity > 1) cout << " HTOF gas layer  :" << HTOF::f_gas_lyr << endl;

  if (TRACKING::FastKalmanFilter)
  {
    /*
       TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + string(Form("ACTIVEGAS_HTOF")),  //      const std::string& phg4hitsNames,
                                               PHG4TrackFastSim::Vertical_Plane,                         //      const DETECTOR_TYPE phg4dettype,
                                               1, //1. / sqrt(12.),                                     //      const float radres,
                                               5.0e-1,//55e-4,                                              //      const float phires,
                                               5.0e-1,//55e-4,                                              //      const float lonres,
                                               1,                                                  //      const float eff,
                                               0);                                                 //      const float noise

       */
    //Reference plane projection at initial R of ToF
    TRACKING::FastKalmanFilter->add_zplane_state(string("G4HIT_") + string(Form("ACTIVEGAS_HTOF")), HTOF::z_end);
    TRACKING::ProjectionNames.insert(string("G4HIT_") + string(Form("ACTIVEGAS_HTOF")));
  }
  return;
}

void HTOF_Reco()
{
  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4detectors.so");

  int verbosity = std::max(Enable::VERBOSITY, Enable::HTOF_VERBOSITY);
  Fun4AllServer* se = Fun4AllServer::instance();
  //se->Verbosity(INT_MAX-10);

  return;
}

#endif
