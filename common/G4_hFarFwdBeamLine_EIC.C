#ifndef MACRO_G4HFARFWDBEAMLINE_EIC_C
#define MACRO_G4HFARFWDBEAMLINE_EIC_C

#include <GlobalVariables.C>

#include <g4detectors/BeamLineMagnetSubsystem.h>
#include <g4detectors/PHG4BlockSubsystem.h>
#include <g4detectors/PHG4ConeSubsystem.h>
#include <g4detectors/PHG4CylinderSubsystem.h>

#include <eicg4zdc/EICG4ZDCHitTree.h>
#include <eicg4zdc/EICG4ZDCNtuple.h>
#include <eicg4zdc/EICG4ZDCSubsystem.h>

#include <eiceval/FarForwardEvaluator.h>

#include <g4main/PHG4Reco.h>

#include <TSystem.h>

#include <fun4all/Fun4AllServer.h>

R__LOAD_LIBRARY(libg4detectors.so)

float PosFlip(float pos);
float AngleFlip(float angle);
float MagFieldFlip(float Bfield);

// This creates the Enable Flag to be used in the main steering macro
namespace Enable
{
  bool HFARFWD_MAGNETS = false;
  bool HFARFWD_VIRTUAL_DETECTORS = false;

  bool HFARFWD_PIPE = false;
  bool HFARFWD_OVERLAPCHECK = false;
  int HFARFWD_VERBOSITY = 0;
  bool ZDC_DISABLE_BLACKHOLE = false;

  //enabled automatically in hFarFwdBeamLineInit(), unless overridden by user
  bool HFARFWD_MAGNETS_IP6 = false;
  bool HFARFWD_MAGNETS_IP8 = false;

  //enabled automatically in hFarFwdBeamLineInit(), unless overridden by user
  bool HFARFWD_VIRTUAL_DETECTORS_IP6 = false;
  bool HFARFWD_VIRTUAL_DETECTORS_IP8 = false;

  float HFARFWD_ION_ENERGY = 0;

  bool FFR_EVAL = false;

}  // namespace Enable

namespace hFarFwdBeamLine
{
  double starting_z = 450;  //cm as center-forward interface
  double enclosure_z_max = NAN;
  double enclosure_r_max = NAN;
  double enclosure_center = NAN;

  PHG4CylinderSubsystem *hFarFwdBeamLineEnclosure(nullptr);

  BeamLineMagnetSubsystem *B0Magnet = (nullptr);
}  // namespace hFarFwdBeamLine

void hFarFwdBeamLineInit()
{
  Enable::HFARFWD_MAGNETS_IP6 |= Enable::HFARFWD_MAGNETS and Enable::IP6;
  Enable::HFARFWD_MAGNETS_IP8 |= Enable::HFARFWD_MAGNETS and Enable::IP8;

  Enable::HFARFWD_VIRTUAL_DETECTORS_IP6 |= Enable::HFARFWD_VIRTUAL_DETECTORS and Enable::IP6;
  Enable::HFARFWD_VIRTUAL_DETECTORS_IP8 |= Enable::HFARFWD_VIRTUAL_DETECTORS and Enable::IP8;

  if (Enable::HFARFWD_MAGNETS_IP6 && Enable::HFARFWD_MAGNETS_IP8)
  {
    cout << "You cannot have magnets for both IP6 and IP8 ON at the same time" << endl;
    gSystem->Exit(1);
  }

  if (Enable::HFARFWD_MAGNETS_IP6)
  {
    hFarFwdBeamLine::enclosure_z_max = 4500.;
    BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, hFarFwdBeamLine::starting_z);
    hFarFwdBeamLine::enclosure_r_max = 200.;
  }

  if (Enable::HFARFWD_MAGNETS_IP8)
  {
    hFarFwdBeamLine::enclosure_z_max = 4500.;
    BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, hFarFwdBeamLine::starting_z);
    hFarFwdBeamLine::enclosure_r_max = 200.;
  }

  hFarFwdBeamLine::enclosure_center = 0.5 * (hFarFwdBeamLine::starting_z + hFarFwdBeamLine::enclosure_z_max);

  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, hFarFwdBeamLine::enclosure_z_max);
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, hFarFwdBeamLine::enclosure_r_max);
}

void hFarFwdDefineMagnets(PHG4Reco *g4Reco)
{
  bool overlapCheck = Enable::OVERLAPCHECK || Enable::HFARFWD_OVERLAPCHECK;
  int verbosity = std::max(Enable::VERBOSITY, Enable::HFARFWD_VERBOSITY);

  hFarFwdBeamLine::hFarFwdBeamLineEnclosure = new PHG4CylinderSubsystem("hFarFwdBeamLineEnclosure");
  hFarFwdBeamLine::hFarFwdBeamLineEnclosure->set_double_param("place_z", hFarFwdBeamLine::enclosure_center);
  hFarFwdBeamLine::hFarFwdBeamLineEnclosure->set_double_param("radius", 0);
  hFarFwdBeamLine::hFarFwdBeamLineEnclosure->set_double_param("thickness", hFarFwdBeamLine::enclosure_r_max);  // This is intentionally made large 25cm radius
  hFarFwdBeamLine::hFarFwdBeamLineEnclosure->set_double_param("length", hFarFwdBeamLine::enclosure_z_max - hFarFwdBeamLine::starting_z);
  hFarFwdBeamLine::hFarFwdBeamLineEnclosure->set_string_param("material", "G4_Galactic");
  hFarFwdBeamLine::hFarFwdBeamLineEnclosure->set_color(.5, .5, .5, 0.2);
  hFarFwdBeamLine::hFarFwdBeamLineEnclosure->OverlapCheck(overlapCheck);
  if (verbosity) hFarFwdBeamLine::hFarFwdBeamLineEnclosure->Verbosity(verbosity);
  g4Reco->registerSubsystem(hFarFwdBeamLine::hFarFwdBeamLineEnclosure);

  if (verbosity > 0)
  {
    std::cout << "hFarFwdBeamLine::hFarFwdBeamLineEnclosure CanBeMotherSubsystem = " << hFarFwdBeamLine::hFarFwdBeamLineEnclosure->CanBeMotherSubsystem() << std::endl;
  }

  string magFile;
  if (Enable::HFARFWD_MAGNETS_IP6)
    magFile = string(getenv("CALIBRATIONROOT")) + "/Beam/ip6_h_farFwdBeamLineMagnets_v2.0.dat";
  else if (Enable::HFARFWD_MAGNETS_IP8)
    magFile = string(getenv("CALIBRATIONROOT")) + "/Beam/ip8_35mrad_h_farFwdBeamLineMagnets.dat";
  else
  {
    cout << " You have to enable either the IP6 or IP8 Magnet configuration to define magnets! " << endl;
    gSystem->Exit(1);
  }

  // make magnet active volume if you want to study the hits
  bool magnet_active = false;
  int absorberactive = 0;

  // if you insert numbers it only displays those magnets, do not comment out the set declaration
  set<int> magnetlist;
  //magnetlist.insert(7);

  BeamLineMagnetSubsystem *bl = nullptr;
  std::ifstream infile(magFile);
  if (infile.is_open())
  {
    double biggest_z = 0.;
    int imagnet = 0;
    std::string line;
    while (std::getline(infile, line))
    {
      if (!line.compare(0, 1, "B") ||
          !line.compare(0, 1, "Q") ||
          !line.compare(0, 1, "S"))
      {
        std::istringstream iss(line);
        string magname;
        double x;
        double y;
        double z;
        double inner_radius_zin;
        double inner_radius_zout;
        double outer_magnet_diameter;
        double length;
        double angle;
        double dipole_field_x;
        double fieldgradient;
        if (!(iss >> magname >> x >> y >> z >> inner_radius_zin >> inner_radius_zout >> outer_magnet_diameter >> length >> angle >> dipole_field_x >> fieldgradient))
        {
          cout << "could not decode " << line << endl;
          gSystem->Exit(1);
        }
        else
        {
	  //------------------------
	  //Select only the magnet component in the far forward region
	  if (z < 0.0)
		continue;

          string magtype;
          if (inner_radius_zin != inner_radius_zout)
          {
            cout << "inner radius at front of magnet " << inner_radius_zin
                 << " not equal radius at back of magnet " << inner_radius_zout
                 << " needs change in code (replace tube by cone for beamline)" << endl;
            gSystem->Exit(1);
          }
          if (verbosity > 0)
          {
            cout << endl
                 << endl
                 << "\tID number " << imagnet << endl;
            cout << "magname: " << magname << endl;
            cout << "x: " << x << endl;
            cout << "y: " << y << endl;
            cout << "z: " << z << endl;
            cout << "inner_radius_zin: " << inner_radius_zin << endl;
            cout << "inner_radius_zout: " << inner_radius_zout << endl;
            cout << "outer_magnet_diameter: " << outer_magnet_diameter << endl;
            cout << "length: " << length << endl;
            cout << "angle: " << angle << endl;
            cout << "dipole_field_x: " << dipole_field_x << endl;
            cout << "fieldgradient: " << fieldgradient << endl;
          }
          if (!magname.compare(0, 1, "B"))
          {
            magtype = "DIPOLE";
          }
          else if (!magname.compare(0, 1, "Q"))
          {
            magtype = "QUADRUPOLE";
          }
          else if (!magname.compare(0, 1, "S"))
          {
            magtype = "SEXTUPOLE";
          }
          else
          {
            cout << "cannot decode magnet name " << magname << endl;
            gSystem->Exit(1);
          }
          // convert to our units (cm, deg)
          x *= 100.;
          y *= 100.;
          z *= 100.;
          length *= 100.;
          inner_radius_zin *= 100.;
          outer_magnet_diameter *= 100.;
          angle = (angle / TMath::Pi() * 180.) / 1000.;  // given in mrad

	  //------------------------
	  // Linearly scaling down the magnetic field for lower energy proton
	  if( Enable::HFARFWD_ION_ENERGY != 275 ) {
             float scaleFactor = Enable::HFARFWD_ION_ENERGY / 275. ;
	     dipole_field_x = dipole_field_x*scaleFactor;
   	  }

          if (magnetlist.empty() || magnetlist.find(imagnet) != magnetlist.end())
          {
            bl = new BeamLineMagnetSubsystem("BEAMLINEMAGNET", imagnet);
            bl->set_double_param("field_y", MagFieldFlip(dipole_field_x));
            bl->set_double_param("fieldgradient", MagFieldFlip(fieldgradient));
            bl->set_string_param("magtype", magtype);
            bl->set_double_param("length", length);
            bl->set_double_param("place_x", PosFlip(x));// relative position to mother vol.
            bl->set_double_param("place_y", y);// relative position to mother vol.
            bl->set_double_param("place_z", z - hFarFwdBeamLine::enclosure_center);// relative position to mother vol.
            bl->set_double_param("field_global_position_x", PosFlip(x));// abs. position to world for field manager
            bl->set_double_param("field_global_position_y", y);// abs. position to world for field manager
            bl->set_double_param("field_global_position_z", z);// abs. position to world for field manager
            bl->set_double_param("rot_y", AngleFlip(angle));
            bl->set_double_param("field_global_rot_y", AngleFlip(angle));// abs. rotation to world for field manager
            bl->set_double_param("inner_radius", inner_radius_zin);
            bl->set_double_param("outer_radius", outer_magnet_diameter / 2.);
            bl->SetActive(magnet_active);
            bl->BlackHole();
            bl->SetMotherSubsystem(hFarFwdBeamLine::hFarFwdBeamLineEnclosure);
            if (absorberactive)
            {
              bl->SetAbsorberActive();
            }
            bl->OverlapCheck(overlapCheck);
            bl->SuperDetector("BEAMLINEMAGNET");
            if (verbosity) bl->Verbosity(verbosity);
            g4Reco->registerSubsystem(bl);

            // rag the B0 magnet
            if (imagnet == 0)
              hFarFwdBeamLine::B0Magnet = bl;
          }
          imagnet++;
          if (fabs(z) + length > biggest_z)
          {
            biggest_z = fabs(z) + length;
          }
        }
      }
    }
    infile.close();
  }
}

void hFarFwdDefineDetectorsIP6(PHG4Reco *g4Reco)
{
  bool overlapCheck = Enable::OVERLAPCHECK || Enable::HFARFWD_OVERLAPCHECK;
  if (Enable::HFARFWD_VIRTUAL_DETECTORS_IP6 && Enable::HFARFWD_VIRTUAL_DETECTORS_IP8)
  {
    cout << "You cannot have detectors enabled for both IP6 and IP8 ON at the same time" << endl;
    gSystem->Exit(1);
  }

  int verbosity = std::max(Enable::VERBOSITY, Enable::HFARFWD_VERBOSITY);

  auto *detZDCsurrogate = new PHG4BlockSubsystem("zdcTruth");
  const double detZDCsurrogate_size_z = 0.1;
  detZDCsurrogate->SuperDetector("ZDCsurrogate");
  detZDCsurrogate->set_double_param("place_x", PosFlip(-96.24));
  detZDCsurrogate->set_double_param("place_y", 0);
  detZDCsurrogate->set_double_param("place_z", 3700 - hFarFwdBeamLine::enclosure_center);
  detZDCsurrogate->set_double_param("rot_y", AngleFlip(0.025 * TMath::RadToDeg()));
  detZDCsurrogate->set_double_param("size_x", 60);
  detZDCsurrogate->set_double_param("size_y", 60);
  detZDCsurrogate->set_double_param("size_z", detZDCsurrogate_size_z);
  detZDCsurrogate->set_string_param("material", "G4_Si");
  detZDCsurrogate->SetActive();
  detZDCsurrogate->set_color(1, 0, 0, 0.5);
  detZDCsurrogate->OverlapCheck(overlapCheck);
  if (!Enable::ZDC_DISABLE_BLACKHOLE) detZDCsurrogate->BlackHole();
  if (verbosity) detZDCsurrogate->Verbosity(verbosity);
  detZDCsurrogate->SetMotherSubsystem(hFarFwdBeamLine::hFarFwdBeamLineEnclosure);
  g4Reco->registerSubsystem(detZDCsurrogate);

  if (Enable::ZDC_DISABLE_BLACKHOLE)
  {
    EICG4ZDCSubsystem *detZDC = new EICG4ZDCSubsystem("EICG4ZDC");
    detZDC->SetActive();
    detZDC->set_double_param("place_z", 3700. + detZDCsurrogate_size_z - hFarFwdBeamLine::enclosure_center);
    detZDC->set_double_param("place_x", PosFlip(-96.24));
    detZDC->set_double_param("rot_y", AngleFlip(0.025));
    detZDC->OverlapCheck(overlapCheck);
    detZDC->SetMotherSubsystem(hFarFwdBeamLine::hFarFwdBeamLineEnclosure);
    g4Reco->registerSubsystem(detZDC);
  }

  const int offMomDetNr = 2;
  const double om_zCent[offMomDetNr] = {3450, 3650};
  const double om_xCent[offMomDetNr] = {-162, -171};
  for (int i = 0; i < offMomDetNr; i++)
  {
    auto *detOM = new PHG4BlockSubsystem(Form("offMomTruth_%d", i), i);
    detOM->SuperDetector("offMomTruth");
    detOM->set_double_param("place_x", PosFlip(om_xCent[i]));
    detOM->set_double_param("place_y", 0);
    detOM->set_double_param("place_z", om_zCent[i] - hFarFwdBeamLine::enclosure_center);
    detOM->set_double_param("rot_y", AngleFlip(0.045 * TMath::RadToDeg()));
    detOM->set_double_param("size_x", 50);
    detOM->set_double_param("size_y", 35);
    detOM->set_double_param("size_z", 0.03);
    detOM->set_string_param("material", "G4_Si");
    detOM->SetActive();
    if (verbosity) detOM->Verbosity(verbosity);
    detOM->OverlapCheck(overlapCheck);
    detOM->SetMotherSubsystem(hFarFwdBeamLine::hFarFwdBeamLineEnclosure);
    g4Reco->registerSubsystem(detOM);
  }

  const int rpDetNr = 2;
  const double rp_zCent[rpDetNr] = {2600, 2800};
  const double rp_xCent[rpDetNr] = {-83.22, -92.20};
  for (int i = 0; i < rpDetNr; i++)
  {
    ////*********************
    //// Square design
    //// 25 cm in x
    //
    //    auto *detRP = new PHG4BlockSubsystem(Form("rpTruth_%d",i));
    ////    detRP->SuperDetector("RomanPots");
    //    detRP->SuperDetector(Form("RomanPots_%d",i));
    //    detRP->set_double_param("place_x",rp_xCent[i]);
    //    detRP->set_double_param("place_y",0);
    //    detRP->set_double_param("place_z",rp_zCent[i]);
    //    detRP->set_double_param("rot_y",-0.025*TMath::RadToDeg());
    //    detRP->set_double_param("size_x",25);
    //    detRP->set_double_param("size_y",10);
    //    detRP->set_double_param("size_z",0.03);
    //    detRP->set_string_param("material","G4_Si");

    ////*********************
    //// Disk design
    //// 50 cm in x

    auto *detRP = new PHG4CylinderSubsystem(Form("rpTruth_%d", i), i);
    detRP->SuperDetector("rpTruth");
    detRP->set_double_param("place_x", PosFlip(rp_xCent[i]));
    detRP->set_double_param("place_y", 0);
    detRP->set_double_param("place_z", rp_zCent[i] - hFarFwdBeamLine::enclosure_center);
    detRP->set_double_param("rot_y", AngleFlip(0.047 * TMath::RadToDeg()));
    detRP->set_double_param("radius", 0);
    detRP->set_double_param("thickness", 25);  // This is intentionally made large 25cm radius
    detRP->set_double_param("length", 0.03);
    detRP->set_string_param("material", "G4_Si");
    detRP->OverlapCheck(overlapCheck);
    detRP->SetMotherSubsystem(hFarFwdBeamLine::hFarFwdBeamLineEnclosure);

    detRP->SetActive();
    if (verbosity) detRP->Verbosity(verbosity);
    g4Reco->registerSubsystem(detRP);
  }

  const int b0DetNr = 4;

  // Sep 09 2021 by Bill: 
  // B0 magnet center location in z: 640
  // B0 place location in z after 50cm shift: 592, 616, 640, 664
  // B0 layers has the same x coordinate: -14.57
 
  const double b0Mag_zCent = 640;
  const double b0Mag_zLen = 120;

  for (int i = 0; i < b0DetNr; i++)
  {
    auto *detB0 = new PHG4CylinderSubsystem(Form("b0Truth_%d", i), i);
    detB0->SuperDetector("b0Truth");
    detB0->set_double_param("radius", 0);
    detB0->set_double_param("thickness", 20);
    detB0->set_double_param("length", 0.1);
    detB0->set_string_param("material", "G4_Si");
    detB0->set_double_param("place_z", b0Mag_zLen / (b0DetNr + 1) * (i - b0DetNr / 2));  // relative to B0 magnet
    detB0->SetActive(true);
    if (verbosity) detB0->Verbosity(verbosity);
    detB0->OverlapCheck(overlapCheck);

    detB0->SetMotherSubsystem(hFarFwdBeamLine::B0Magnet);

    g4Reco->registerSubsystem(detB0);
  }
}

void hFarFwdDefineDetectorsIP8(PHG4Reco *g4Reco)
{

//--------------------------------------------------------
// The IP8 detector position is implemented by Wenliang Li (billlee@jlab.org)
// on July 07, 2021
// Reference of this implementation: https://indico.bnl.gov/event/10974/contributions/51160/

  bool overlapCheck = Enable::OVERLAPCHECK || Enable::HFARFWD_OVERLAPCHECK;
  if (Enable::HFARFWD_VIRTUAL_DETECTORS_IP6 && Enable::HFARFWD_VIRTUAL_DETECTORS_IP8)
  {
    cout << "You cannot have detectors enabled for both IP6 and IP8 ON at the same time" << endl;
    gSystem->Exit(1);
  }

  int verbosity = std::max(Enable::VERBOSITY, Enable::HFARFWD_VERBOSITY);

  const int offMomDetNr = 2;
  const double om_xCent[offMomDetNr] = {46, 49};
  const double om_zCent[offMomDetNr] = {3250, 3450};

  for (int i = 0; i < offMomDetNr; i++)
  {
    auto *detOM = new PHG4BlockSubsystem(Form("offMomTruth_%d", i), i);
    detOM->SuperDetector("offMomTruth");
    detOM->set_double_param("place_x", PosFlip(om_xCent[i]));
    detOM->set_double_param("place_y", 0);
    detOM->set_double_param("place_z", PosFlip(om_zCent[i] - hFarFwdBeamLine::enclosure_center));
    detOM->set_double_param("rot_y", AngleFlip(-0.045 * TMath::RadToDeg()));
    detOM->set_double_param("size_x", 40);  // Original design specification
    detOM->set_double_param("size_y", 35);  // Original design specification
    detOM->set_double_param("size_z", 0.03);
    detOM->set_string_param("material", "G4_Si");
    detOM->OverlapCheck(overlapCheck);
    detOM->SetMotherSubsystem(hFarFwdBeamLine::hFarFwdBeamLineEnclosure);
    detOM->SetActive();
    detOM->set_color(0, 0, 1, 0.5);
    if (verbosity) detOM->Verbosity(verbosity);
    g4Reco->registerSubsystem(detOM);
  }

  auto *detZDCsurrogate = new PHG4BlockSubsystem("zdcTruth");
  const double detZDCsurrogate_size_z = 0.1;
  detZDCsurrogate->SuperDetector("ZDCsurrogate");
  detZDCsurrogate->set_double_param("place_x", PosFlip(120));
  detZDCsurrogate->set_double_param("place_y", 0);
  detZDCsurrogate->set_double_param("place_z", 3350 - hFarFwdBeamLine::enclosure_center);
  detZDCsurrogate->set_double_param("rot_y", AngleFlip(-0.035 * TMath::RadToDeg()));
  detZDCsurrogate->set_double_param("size_x", 60);
  detZDCsurrogate->set_double_param("size_y", 60);
  detZDCsurrogate->set_double_param("size_z", detZDCsurrogate_size_z);
  detZDCsurrogate->set_string_param("material", "G4_Si");
  detZDCsurrogate->SetActive();
  detZDCsurrogate->OverlapCheck(overlapCheck);
  detZDCsurrogate->set_color(1, 0, 0, 0.5);
  if (!Enable::ZDC_DISABLE_BLACKHOLE) detZDCsurrogate->BlackHole();
  if (verbosity)
    detZDCsurrogate->Verbosity(verbosity);
  detZDCsurrogate->SetMotherSubsystem(hFarFwdBeamLine::hFarFwdBeamLineEnclosure);
  g4Reco->registerSubsystem(detZDCsurrogate);

  if (Enable::ZDC_DISABLE_BLACKHOLE)
  {

    EICG4ZDCSubsystem *detZDC = new EICG4ZDCSubsystem("EICG4ZDC");
    detZDC->SetActive();
    detZDC->set_double_param("place_z", 3350. + detZDCsurrogate_size_z - hFarFwdBeamLine::enclosure_center);
    detZDC->set_double_param("place_x", PosFlip(120.0));
    detZDC->set_double_param("rot_y", AngleFlip(-0.035));
    detZDC->SetMotherSubsystem(hFarFwdBeamLine::hFarFwdBeamLineEnclosure);
    detZDC->OverlapCheck(overlapCheck);
    g4Reco->registerSubsystem(detZDC);

  }

  //------------------
  // Roman pot set #1
  const int rpDetNr = 2;

  const double rp_xCent[rpDetNr] = {75.6, 78.15};
  const double rp_zCent[rpDetNr] = {2600, 2800};

  for (int i = 0; i < rpDetNr; i++)
  {
    auto *detRP = new PHG4BlockSubsystem(Form("rpTruth_%d", i), i);
    detRP->SuperDetector("rpTruth");
    detRP->set_double_param("place_x", PosFlip(rp_xCent[i]));
    detRP->set_double_param("place_y", 0);
    detRP->set_double_param("place_z", rp_zCent[i] - hFarFwdBeamLine::enclosure_center);
    detRP->set_double_param("rot_y", AngleFlip(-0.035 * TMath::RadToDeg()));
    detRP->set_double_param("size_x", 25);  // Original design specification
    detRP->set_double_param("size_y", 20);  // Original design specification
    detRP->set_double_param("size_z", 0.03);
    detRP->set_string_param("material", "G4_Si");
    detRP->OverlapCheck(overlapCheck);
    detRP->SetMotherSubsystem(hFarFwdBeamLine::hFarFwdBeamLineEnclosure);
    detRP->SetActive();
    if (verbosity)
      detRP->Verbosity(verbosity);
    g4Reco->registerSubsystem(detRP);
  }

  //------------------
  // Roman pot set #2 before and after the secondary focus

  const int rp2ndDetNr = 2;
  const double rp_2nd_xCent[rp2ndDetNr] = {101.94, 106.94};
  const double rp_2nd_zCent[rp2ndDetNr] = {4300, 4450};

  for (int i = 0; i < rp2ndDetNr; i++)
  {
    auto *detRP_2nd = new PHG4BlockSubsystem(Form("rpTruth2_%d", i), i);
    detRP_2nd->SuperDetector("rpTruth2");
    detRP_2nd->set_double_param("place_x", PosFlip(rp_2nd_xCent[i]));
    detRP_2nd->set_double_param("place_y", 0);
    detRP_2nd->set_double_param("place_z", rp_2nd_zCent[i] - hFarFwdBeamLine::enclosure_center);
    detRP_2nd->set_double_param("rot_y", AngleFlip(-0.029 * TMath::RadToDeg()));
    detRP_2nd->set_double_param("size_x", 25);
    detRP_2nd->set_double_param("size_y", 20);
    detRP_2nd->set_double_param("size_z", 0.03);
    detRP_2nd->set_string_param("material", "G4_Si");
    detRP_2nd->OverlapCheck(overlapCheck);
    detRP_2nd->SetMotherSubsystem(hFarFwdBeamLine::hFarFwdBeamLineEnclosure);
    detRP_2nd->SetActive();
    if (verbosity)
      detRP_2nd->Verbosity(verbosity);
    g4Reco->registerSubsystem(detRP_2nd);
  }

  if (verbosity > 0)
  {
    std::cout << "B0Magnet can be mother = " << hFarFwdBeamLine::B0Magnet->CanBeMotherSubsystem() << std::endl;
  }

  const int b0DetNr = 4;
  const double b0Mag_zCent = 610;
  const double b0Mag_zLen = 120;
  for (int i = 0; i < b0DetNr; i++)
  {
    auto *detB0 = new PHG4CylinderSubsystem(Form("b0Truth_%d", i), i);
    detB0->SuperDetector("b0Truth");
    detB0->set_double_param("radius", 0);
    detB0->set_double_param("thickness", 20);
    detB0->set_double_param("length", 0.1);
    detB0->set_string_param("material", "G4_Si");
    detB0->set_double_param("place_y", 0);
    detB0->set_double_param("place_z", b0Mag_zLen / (b0DetNr + 1) * (i - b0DetNr / 2));  
    detB0->OverlapCheck(overlapCheck);
    detB0->SetMotherSubsystem(hFarFwdBeamLine::B0Magnet);
    detB0->SetActive(true);
    if (verbosity)
      detB0->Verbosity(verbosity);
    g4Reco->registerSubsystem(detB0);
  }
}

void hFarFwdDefineBeamPipe(PHG4Reco *g4Reco)
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::HFARFWD_VERBOSITY);
  //exit window
  PHG4CylinderSubsystem *exitWin = new PHG4CylinderSubsystem("exitWin", 0);
  exitWin->set_double_param("radius", 3.2);
  exitWin->set_double_param("thickness", 11.8);
  exitWin->set_double_param("length", 0.15);
  exitWin->set_double_param("rot_y", AngleFlip(-0.025 * TMath::RadToDeg()));
  exitWin->set_string_param("material", "G4_STAINLESS-STEEL");
  exitWin->set_double_param("place_x", PosFlip(12.5));
  exitWin->set_double_param("place_y", 0);
  exitWin->set_double_param("place_z", 500);
  exitWin->SetActive(false);
  g4Reco->registerSubsystem(exitWin);

  //B0 magnet pipe
  PHG4CylinderSubsystem *pipeB0 = new PHG4CylinderSubsystem("beamPipeB0", 0);
  pipeB0->set_double_param("radius", 2.8);
  pipeB0->set_double_param("thickness", 0.25);
  pipeB0->set_double_param("length", 195);
  pipeB0->set_double_param("rot_y", AngleFlip(-0.025 * TMath::RadToDeg()));
  pipeB0->set_string_param("material", "G4_Al");
  pipeB0->set_double_param("place_x", PosFlip(14.748));
  pipeB0->set_double_param("place_y", 0);
  pipeB0->set_double_param("place_z", 590);
  pipeB0->SetActive(false);
  g4Reco->registerSubsystem(pipeB0);

  //Quad pipes
  const int nSecQ = 5;  //B0apf, Q1apf, Q1bpf, Q2pf, B1pf
  const string nm[nSecQ] = {"B0apf", "Q1apf", "Q1bpf", "Q2pf", "B1pf"};
  const double qlen[nSecQ] = {160, 150, 220, 440, 330};
  const double qir[nSecQ] = {4, 5.1, 7, 12, 12.2};
  const double qor[nSecQ] = {4.1, 5.2, 7.2, 12.2, 12.4};
  const double qrot[nSecQ] = {25, 19.5, 15, 15, 34};  //mrad
  const double qxC[nSecQ] = {19.8, 24.47, 30.05, 39.5, 48};
  const double qyC[nSecQ] = {0, 0, 0, 0, 0};
  const double qzC[nSecQ] = {770, 922.8, 1106.3, 1416.7, 1806.7};
  for (int i = 0; i < nSecQ; i++)
  {
    PHG4CylinderSubsystem *pipe = new PHG4CylinderSubsystem(Form("beamPipe%s", nm[i].c_str()), 0);
    pipe->set_double_param("radius", qir[i]);
    pipe->set_double_param("thickness", qor[i] - qir[i]);
    pipe->set_double_param("length", qlen[i]);
    pipe->set_double_param("rot_y", AngleFlip(-qrot[i] / 1000 * TMath::RadToDeg()));
    pipe->set_string_param("material", "G4_Al");
    pipe->set_double_param("place_x", PosFlip(qxC[i]));
    pipe->set_double_param("place_y", qyC[i]);
    pipe->set_double_param("place_z", qzC[i]);
    pipe->SetActive(false);
    g4Reco->registerSubsystem(pipe);
  }

  //Electron pipe
  PHG4CylinderSubsystem *pipeElectron = new PHG4CylinderSubsystem("beamPipeElectron", 0);
  pipeElectron->set_double_param("radius", 1);
  pipeElectron->set_double_param("thickness", 1);
  pipeElectron->set_double_param("length", 3000);
  pipeElectron->set_double_param("rot_y", AngleFlip(-0.025 * TMath::RadToDeg()));
  pipeElectron->set_string_param("material", "G4_Al");
  pipeElectron->set_double_param("place_x", PosFlip(0));
  pipeElectron->set_double_param("place_y", 0);
  pipeElectron->set_double_param("place_z", 2000);
  pipeElectron->SetActive(false);
  //g4Reco->registerSubsystem(pipeElectron);

  //ZDC pipe
  PHG4CylinderSubsystem *pipeZDC = new PHG4CylinderSubsystem("beamPipeZDC", 0);
  pipeZDC->set_double_param("radius", 16.5);
  pipeZDC->set_double_param("thickness", 0.1);
  pipeZDC->set_double_param("length", 170);
  pipeZDC->set_double_param("rot_y", AngleFlip(-0.025 * TMath::RadToDeg()));
  pipeZDC->set_string_param("material", "G4_Al");
  pipeZDC->set_double_param("place_x", PosFlip(59));
  pipeZDC->set_double_param("place_y", 0);
  pipeZDC->set_double_param("place_z", 2041.59);
  pipeZDC->SetActive(false);
  g4Reco->registerSubsystem(pipeZDC);

  //Roman Pot pipe
  const int nSec = 2;
  const double len[nSec] = {850, 1150};
  const double ir1[nSec] = {17, 17};
  const double or1[nSec] = {17.1, 17.1};
  const double ir2[nSec] = {17, 7};
  const double or2[nSec] = {17.1, 7.1};
  const double xC[nSec] = {83, 130};
  const double yC[nSec] = {0, 0};
  const double zC[nSec] = {2550, 3550};
  for (int i = 0; i < nSec; i++)
  {
    PHG4ConeSubsystem *pipe = new PHG4ConeSubsystem(Form("beamPipeRP%d", i), 0);
    pipe->set_string_param("material", "G4_STAINLESS-STEEL");
    pipe->set_double_param("place_x", PosFlip(xC[i]));
    pipe->set_double_param("place_y", yC[i]);
    pipe->set_double_param("place_z", zC[i]);
    pipe->set_double_param("length", len[i] / 2);
    pipe->set_double_param("rmin1", ir1[i]);
    pipe->set_double_param("rmin2", ir2[i]);
    pipe->set_double_param("rmax1", or1[i]);
    pipe->set_double_param("rmax2", or2[i]);
    pipe->set_double_param("rot_y", AngleFlip(-0.047 * TMath::RadToDeg()));
    g4Reco->registerSubsystem(pipe);
  }
}

float PosFlip(float pos) {
  if(Enable::HFARFWD_MAGNETS_IP6) {
  	return pos;
  } else {
  	return pos;
  }
}

float AngleFlip(float angle){
  if(Enable::HFARFWD_MAGNETS_IP6) {
  	return angle;
  } else {
  	return angle;
  }
}

float MagFieldFlip(float Bfield){
  if(Enable::HFARFWD_MAGNETS_IP6) {
  	return Bfield;
  } else {
  	return Bfield;
  }
}


//------------------------------------------

void FFR_Eval(const std::string &outputfile)
{

  string ip_str;


  if(Enable::IP6) {
    ip_str = "IP6";
  } else {
    ip_str = "IP8";
  }

  int verbosity = std::max(Enable::VERBOSITY, Enable::EEMC_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();

  FarForwardEvaluator *eval = new FarForwardEvaluator("FARFORWARDEVALUATOR", "FFR", outputfile.c_str(), ip_str);

  eval->Verbosity(verbosity);
  se->registerSubsystem(eval);

  return;
}



#endif
