#ifndef MACRO_G4BToF_C
#define MACRO_G4BToF_C

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
  bool BTOF = false;
  bool BTOF_OVERLAPCHECK = false;
  int BTOF_VERBOSITY = 0;
}  // namespace Enable

namespace BTOF
{
  const int gas_lyr = 6;       // 1/2 of the total number of layers
  const int mrpc_inn_lyr = 5;  // 1/2 of the total number of layers
  const double rad = 76.5;     // cm
  const double zpos = -60.0;   // cm
  const double length = 400.;  //cm
  int subsysID = 0;
}  // namespace BTOF

void BToFInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 85.);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, BTOF::length / 2.0 + 10.0);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -BTOF::length / 2.0 - 10.0);
}

double Build_G4_BTof(PHG4Reco* g4Reco,
                     double tof_rad = 82.0,
                     double zpos = -60.0,
                     double tof_length = 400.0)
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::BTOF_VERBOSITY);
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::BTOF_OVERLAPCHECK;

  gSystem->Load("libfun4all");
  gSystem->Load("libg4detectors.so");
  gSystem->Load("libg4trackfastsim.so");

  double rsum = 0.0;

  //mRPC TOF parameters
  double gas_gap = 0.0220;            // 220 microns , 12 gas gaps
  double mrpc_in_thick = 0.04;        // 400 microns, 10 think glass
  double mrpc_out_thick = 0.07;       // 700 microns, 4 thick glass
  double pcb_thickness = 0.06;        // 600 microns, 3 PCBs
  double cu_thickness = 0.003;        // 30 microns, layer over pcb, (1 each on outer pcb and 2 on central pcb = 4)
  double carbon_thickness = 0.01;     // 100 microns , 2 layers
  double mylar_thickness = 0.04;      // 400 microns, 4 layers
  double honeycomb_thickness = 0.75;  // 7.5 mm, 2 honeycomb

  //-------------------------------
  // mRPC material needs to change to plate glass (rad length 10.69 cm).
  //Currently implemented G4_Si has rad length of 9.37 cm.
  //Changing material to plate glass will increase radiation length by 10% w.r.t. to 1.X0(~10%)
  // and the new radiation length will be 1.1X0 (~ 11%).Assuming G4_Si for mRPC is reasonable for the time being
  //-----------------------------------

  PHG4CylinderSubsystem* tof_cyl;

  //Honeycomb
  tof_cyl = new PHG4CylinderSubsystem("ToF_honeycomb_bottom", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", tof_rad);
  tof_cyl->set_string_param("material", "NOMEX");
  tof_cyl->set_double_param("thickness", honeycomb_thickness);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  g4Reco->registerSubsystem(tof_cyl);
  ++BTOF::subsysID;
  if (verbosity > 1) cout << " bottom HC :" << tof_rad << endl;

  //PCB
  rsum = tof_rad + honeycomb_thickness;
  tof_cyl = new PHG4CylinderSubsystem("ToF_pcb_bottom", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "FR4");
  tof_cyl->set_double_param("thickness", pcb_thickness);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  g4Reco->registerSubsystem(tof_cyl);
  ++BTOF::subsysID;
  if (verbosity > 1) cout << " bototm PCB " << rsum << endl;

  //PCB Cu
  rsum += pcb_thickness;
  tof_cyl = new PHG4CylinderSubsystem("ToF_bottompcb_cu", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "G4_Cu");
  tof_cyl->set_double_param("thickness", cu_thickness);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);  // Length restricted to active area
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  g4Reco->registerSubsystem(tof_cyl);
  ++BTOF::subsysID;
  if (verbosity > 1) cout << " PCB Cu :" << rsum << endl;

  //Mylar
  rsum += cu_thickness;
  tof_cyl = new PHG4CylinderSubsystem("ToF_mylar_bottom", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "G4_MYLAR");
  tof_cyl->set_double_param("thickness", mylar_thickness);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  g4Reco->registerSubsystem(tof_cyl);
  ++BTOF::subsysID;
  if (verbosity > 1) cout << " Mylar : " << rsum << endl;

  //Carbon layer
  rsum += mylar_thickness;
  tof_cyl = new PHG4CylinderSubsystem("ToF_Carbon_bottom", BTOF::subsysID);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "G4_C");
  tof_cyl->set_double_param("thickness", carbon_thickness);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  g4Reco->registerSubsystem(tof_cyl);
  ++BTOF::subsysID;
  if (verbosity > 1) cout << " Carbon :" << rsum << endl;

  //Outside Glass layer
  rsum += carbon_thickness;
  tof_cyl = new PHG4CylinderSubsystem("ToF_glass_bottom", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "G4_Si");
  tof_cyl->set_double_param("thickness", mrpc_out_thick);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  ++BTOF::subsysID;
  g4Reco->registerSubsystem(tof_cyl);

  if (verbosity > 1) cout << " Glass outside :" << rsum << endl;

  //Active  gas gap abd mRPC inner glass layers
  double rsum_gasin_begin = rsum + mrpc_out_thick;
  double rsum_innglass_begin = rsum_gasin_begin + gas_gap;
  double rsum_gasin;
  double rsum_innglass;
  for (int layer = 0; layer < BTOF::gas_lyr; layer++)
  {
    rsum_gasin = rsum_gasin_begin + layer * gas_gap + layer * mrpc_in_thick;
    tof_cyl = new PHG4CylinderSubsystem(Form("ToF_gas_%d", layer), BTOF::subsysID);
    tof_cyl->Verbosity(verbosity);
    tof_cyl->set_double_param("radius", rsum_gasin);
    tof_cyl->set_string_param("material", "G4_Ar");
    tof_cyl->set_double_param("thickness", gas_gap);
    tof_cyl->set_int_param("lengthviarapidity", 0);
    tof_cyl->set_double_param("place_z", zpos);
    tof_cyl->set_double_param("length", tof_length);
    tof_cyl->SuperDetector("bTOF");
    tof_cyl->SetActive(1);
    tof_cyl->OverlapCheck(OverlapCheck);
    ++BTOF::subsysID;
    g4Reco->registerSubsystem(tof_cyl);
    if (verbosity > 1) cout << " gas inner" << layer << " :  " << rsum_gasin << endl;

    if (layer < BTOF::mrpc_inn_lyr)
    {
      rsum_innglass = rsum_gasin_begin + layer * mrpc_in_thick + (layer + 1) * gas_gap;
      tof_cyl = new PHG4CylinderSubsystem(Form("ToF_inner_glass_%d", layer), BTOF::subsysID);
      tof_cyl->Verbosity(verbosity);
      tof_cyl->set_double_param("radius", rsum_innglass);
      tof_cyl->set_string_param("material", "G4_Si");
      tof_cyl->set_double_param("thickness", mrpc_in_thick);
      tof_cyl->set_int_param("lengthviarapidity", 0);
      tof_cyl->set_double_param("place_z", zpos);
      tof_cyl->set_double_param("length", tof_length);
      tof_cyl->SuperDetector("bTOF");
      tof_cyl->SetActive(0);
      tof_cyl->OverlapCheck(OverlapCheck);
      ++BTOF::subsysID;
      g4Reco->registerSubsystem(tof_cyl);
      if (verbosity > 1) cout << " inner mrpc" << layer << " :  " << rsum_innglass << endl;
    }
  }
  rsum = rsum_gasin;

  //Middle thick  Glass layer
  rsum += gas_gap;
  tof_cyl = new PHG4CylinderSubsystem("ToF_glass_bot_mid", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "G4_Si");
  tof_cyl->set_double_param("thickness", mrpc_out_thick);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  ++BTOF::subsysID;
  g4Reco->registerSubsystem(tof_cyl);
  if (verbosity > 1) cout << " mid bottom glass  :" << rsum << endl;

  //Middle Carbon layer
  rsum += mrpc_out_thick;
  tof_cyl = new PHG4CylinderSubsystem("ToF_Carbon_bot_mid", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "G4_C");
  tof_cyl->set_double_param("thickness", carbon_thickness);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  ++BTOF::subsysID;
  g4Reco->registerSubsystem(tof_cyl);
  if (verbosity > 1) cout << " mid bottom Carbon :" << rsum << endl;

  //Middle mylar layer

  rsum += carbon_thickness;
  tof_cyl = new PHG4CylinderSubsystem("ToF_mylar_bot_mid", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "G4_MYLAR");
  tof_cyl->set_double_param("thickness", mylar_thickness);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  ++BTOF::subsysID;
  g4Reco->registerSubsystem(tof_cyl);
  if (verbosity > 1) cout << " Mid bottom Mylar : " << rsum << endl;

  //Middle PCB bottom Cu
  rsum += mylar_thickness /*cu_thickness*/;
  tof_cyl = new PHG4CylinderSubsystem("ToF_midpcb_cub", BTOF::subsysID);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "G4_Cu");
  tof_cyl->set_double_param("thickness", cu_thickness);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);  // Length restricted to active area
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  g4Reco->registerSubsystem(tof_cyl);
  if (verbosity > 1) cout << " mid PCB bottom cu :" << rsum << endl;

  //Mid PCB
  rsum += cu_thickness /*pcb_thickness*/;
  tof_cyl = new PHG4CylinderSubsystem("ToF_pcb_mid", BTOF::subsysID);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "FR4");
  tof_cyl->set_double_param("thickness", pcb_thickness);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  ++BTOF::subsysID;
  g4Reco->registerSubsystem(tof_cyl);
  if (verbosity > 1) cout << " mid PCB " << rsum << endl;

  //Middle PCB top Cu
  rsum += pcb_thickness;
  tof_cyl = new PHG4CylinderSubsystem("ToF_midpcb_tcu", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "G4_Cu");
  tof_cyl->set_double_param("thickness", cu_thickness);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);  // Length restricted to active area
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  ++BTOF::subsysID;
  g4Reco->registerSubsystem(tof_cyl);
  if (verbosity > 1) cout << " mid PCB top cu :" << rsum << endl;

  //Middle top mylar layer

  rsum += cu_thickness;
  tof_cyl = new PHG4CylinderSubsystem("ToF_mylar_topmid", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "G4_MYLAR");
  tof_cyl->set_double_param("thickness", mylar_thickness);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  ++BTOF::subsysID;
  g4Reco->registerSubsystem(tof_cyl);
  if (verbosity > 1) cout << " Mid top  Mylar : " << rsum << endl;

  //Middle top carbon layer
  rsum += mylar_thickness;
  tof_cyl = new PHG4CylinderSubsystem("ToF_Carbon_top_mid", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "G4_C");
  tof_cyl->set_double_param("thickness", carbon_thickness);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  ++BTOF::subsysID;
  g4Reco->registerSubsystem(tof_cyl);
  if (verbosity > 1) cout << " mid top Carbon :" << rsum << endl;

  //Middle top thick glass
  rsum += carbon_thickness;
  tof_cyl = new PHG4CylinderSubsystem("ToF_glass_top_mid", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "G4_Si");
  tof_cyl->set_double_param("thickness", mrpc_out_thick);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  ++BTOF::subsysID;
  g4Reco->registerSubsystem(tof_cyl);
  if (verbosity > 1) cout << " mid top glass  :" << rsum << endl;

  //Upper half gas gaps and mRPCs
  //Active  gas gap abd mRPC inner glass layers
  double rsum_gasin_begin_uphalf = (rsum + mrpc_out_thick);
  double rsum_innglass_begin_uphalf = (rsum_gasin_begin_uphalf + gas_gap);
  double rsum_gasin_uphalf;
  double rsum_innglass_uphalf;
  for (int layer = 0; layer < BTOF::gas_lyr; layer++)
  {
    rsum_gasin_uphalf = rsum_gasin_begin_uphalf + layer * gas_gap + layer * mrpc_in_thick;
    tof_cyl = new PHG4CylinderSubsystem(Form("ToF_gas_%d", layer + 6), BTOF::subsysID);  //stupid way to implement indexing of active gas layer for the upper half of ToF
    tof_cyl->Verbosity(verbosity);
    tof_cyl->set_double_param("radius", rsum_gasin_uphalf);
    tof_cyl->set_string_param("material", "G4_Ar");
    tof_cyl->set_double_param("thickness", gas_gap);
    tof_cyl->set_int_param("lengthviarapidity", 0);
    tof_cyl->set_double_param("place_z", zpos);
    tof_cyl->set_double_param("length", tof_length);
    tof_cyl->SuperDetector("bTOF");
    tof_cyl->SetActive(1);
    tof_cyl->OverlapCheck(OverlapCheck);
    ++BTOF::subsysID;
    g4Reco->registerSubsystem(tof_cyl);
    if (verbosity > 1) cout << " gas inner" << layer << " :  " << rsum_gasin_uphalf << endl;

    if (layer < BTOF::mrpc_inn_lyr)
    {
      rsum_innglass_uphalf = rsum_gasin_begin_uphalf + layer * mrpc_in_thick + (layer + 1) * gas_gap;
      tof_cyl = new PHG4CylinderSubsystem(Form("ToF_inner_up_half_glass_%d", layer), BTOF::subsysID);
      tof_cyl->Verbosity(verbosity);
      tof_cyl->set_double_param("radius", rsum_innglass_uphalf);
      tof_cyl->set_string_param("material", "G4_Si");
      tof_cyl->set_double_param("thickness", mrpc_in_thick);
      tof_cyl->set_int_param("lengthviarapidity", 0);
      tof_cyl->set_double_param("place_z", zpos);
      tof_cyl->set_double_param("length", tof_length);
      tof_cyl->SuperDetector("bTOF");
      tof_cyl->SetActive(0);
      tof_cyl->OverlapCheck(OverlapCheck);
      ++BTOF::subsysID;
      g4Reco->registerSubsystem(tof_cyl);
      if (verbosity > 1) cout << " inner mrpc" << layer << " :  " << rsum_innglass_uphalf << endl;
    }
  }
  rsum = rsum_gasin_uphalf;

  //top outer glass
  rsum += gas_gap;
  tof_cyl = new PHG4CylinderSubsystem("ToF_glass_top", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "G4_Si");
  tof_cyl->set_double_param("thickness", mrpc_out_thick);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  ++BTOF::subsysID;
  g4Reco->registerSubsystem(tof_cyl);
  if (verbosity > 1) cout << " top glass  :" << rsum << endl;

  // Carbon layer
  rsum += mrpc_out_thick;
  tof_cyl = new PHG4CylinderSubsystem("ToF_Carbon_top", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "G4_C");
  tof_cyl->set_double_param("thickness", carbon_thickness);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  ++BTOF::subsysID;
  g4Reco->registerSubsystem(tof_cyl);
  if (verbosity > 1) cout << " Carbon :" << rsum << endl;

  //Mylar
  rsum += carbon_thickness;
  tof_cyl = new PHG4CylinderSubsystem("ToF_mylar_top", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "G4_MYLAR");
  tof_cyl->set_double_param("thickness", mylar_thickness);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  ++BTOF::subsysID;
  g4Reco->registerSubsystem(tof_cyl);
  if (verbosity > 1) cout << "Mylar :" << rsum << endl;

  //PCB Cu
  rsum += mylar_thickness;
  tof_cyl = new PHG4CylinderSubsystem("ToF_toppcb_cu", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "G4_Cu");
  tof_cyl->set_double_param("thickness", cu_thickness);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);  // Length restricted to active area
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  ++BTOF::subsysID;
  g4Reco->registerSubsystem(tof_cyl);
  if (verbosity > 1) cout << " PCB cu :" << rsum << endl;

  //Top PCB
  rsum += cu_thickness;
  tof_cyl = new PHG4CylinderSubsystem("ToF_pcb_top", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", rsum);
  tof_cyl->set_string_param("material", "FR4");
  tof_cyl->set_double_param("thickness", pcb_thickness);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  ++BTOF::subsysID;
  g4Reco->registerSubsystem(tof_cyl);
  if (verbosity > 1) cout << " top PCB " << rsum << endl;

  //Honeycomb
  rsum += pcb_thickness;
  tof_cyl = new PHG4CylinderSubsystem("ToF_honeycomb_top", BTOF::subsysID);
  tof_cyl->Verbosity(verbosity);
  tof_cyl->set_double_param("radius", tof_rad);
  tof_cyl->set_string_param("material", "NOMEX");
  tof_cyl->set_double_param("thickness", honeycomb_thickness);
  tof_cyl->SetActive(0);
  tof_cyl->set_int_param("lengthviarapidity", 0);
  tof_cyl->set_double_param("place_z", zpos);
  tof_cyl->set_double_param("length", tof_length);
  tof_cyl->SuperDetector("bTOF");
  tof_cyl->SetActive(0);
  tof_cyl->OverlapCheck(OverlapCheck);
  ++BTOF::subsysID;
  g4Reco->registerSubsystem(tof_cyl);
  if (verbosity > 1) cout << " Top honeycomb :" << rsum << endl;

  if (verbosity > 1) cout << " tof thickness :" << tof_rad << endl;
  return tof_rad;
}

double BToFSetup(PHG4Reco* g4Reco)
{
  double radius = 0;

  radius = Build_G4_BTof(g4Reco,
                         BTOF::rad,
                         BTOF::zpos,
                         BTOF::length);

  for (int igap = 0; igap < BTOF::gas_lyr + 6; igap++)
  {
    //Uncomment below if one wants tracking evaluation
    /*
      if (TRACKING::FastKalmanFilter)
      TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + string(Form("ToF_gas_%d", igap)),  //      const std::string& phg4hitsNames,
      PHG4TrackFastSim::Cylinder,                         //      const DETECTOR_TYPE phg4dettype,
      1, //1. / sqrt(12.),                                     //      const float radres,
      5.0e-1,//55e-4,                                              //      const float phires,
      5.0e-1,//55e-4,                                              //      const float lonres,
      1,                                                  //      const float eff,
      0);                                                 //      const float noise

    */

    //Reference plane projection at initial R of ToF
    TRACKING::FastKalmanFilter->add_cylinder_state(string("G4HIT_") + string(Form("ACTIVEGAS_BTOF_%d", igap)), 82.0);
    TRACKING::ProjectionNames.insert(string("G4HIT_") + string(Form("ACTIVEGAS_BTOF_%d", igap)));
  }
  return radius;  // cm
}

void BToF_Reco(int verbosity = 0)
{
  //---------------
  // Load libraries
  //---------------
  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4detectors.so");

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer* se = Fun4AllServer::instance();

  return;
}

#endif
