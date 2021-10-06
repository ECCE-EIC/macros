#ifndef MACRO_G4HFARBWDBEAMLINE_EIC_C
#define MACRO_G4HFARBWDBEAMLINE_EIC_C

#include <GlobalVariables.C>

#include <g4detectors/BeamLineMagnetSubsystem.h>
#include <g4detectors/PHG4BlockSubsystem.h>
#include <g4detectors/PHG4ConeSubsystem.h>
#include <g4detectors/PHG4CylinderSubsystem.h>

#include <eicg4zdc/EICG4ZDCHitTree.h>
#include <eicg4zdc/EICG4ZDCNtuple.h>
#include <eicg4zdc/EICG4ZDCSubsystem.h>

#include <g4main/PHG4Reco.h>

#include <TSystem.h>

R__LOAD_LIBRARY(libg4detectors.so)

// This creates the Enable Flag to be used in the main steering macro
namespace Enable
{
  bool HFARBWD_MAGNETS = false;
  bool HFARBWD_VIRTUAL_DETECTORS = false;

  bool HFARBWD_PIPE = false;
  bool HFARBWD_OVERLAPCHECK = false;
  int HFARBWD_VERBOSITY = 0;

  //enabled automatically in hFarBwdBeamLineInit(), unless overridden by user
  bool HFARBWD_MAGNETS_IP6 = false;
  bool HFARBWD_MAGNETS_IP8 = false;

  //enabled automatically in hFarBwdBeamLineInit(), unless overridden by user
  bool HFARBWD_VIRTUAL_DETECTORS_IP6 = false;
  bool HFARBWD_VIRTUAL_DETECTORS_IP8 = false;

}  // namespace Enable

namespace hFarBwdBeamLine
{
  double starting_z = -463;  //cm clear the backward beam chamber
  double enclosure_z_max = NAN;
  double enclosure_r_max = NAN;
  double enclosure_center = NAN;

  PHG4CylinderSubsystem *hFarBwdBeamLineEnclosure(nullptr);

  BeamLineMagnetSubsystem *B0Magnet = (nullptr);
}  // namespace hFarBwdBeamLine

void hFarBwdBeamLineInit()
{
  Enable::HFARBWD_MAGNETS_IP6 |= Enable::HFARBWD_MAGNETS and Enable::IP6;
  Enable::HFARBWD_MAGNETS_IP8 |= Enable::HFARBWD_MAGNETS and Enable::IP8;

  Enable::HFARBWD_VIRTUAL_DETECTORS_IP6 |= Enable::HFARBWD_VIRTUAL_DETECTORS and Enable::IP6;
  Enable::HFARBWD_VIRTUAL_DETECTORS_IP8 |= Enable::HFARBWD_VIRTUAL_DETECTORS and Enable::IP8;

  if (Enable::HFARBWD_MAGNETS_IP6 && Enable::HFARBWD_MAGNETS_IP8)
  {
    cout << "You cannot have magnets for both IP6 and IP8 ON at the same time" << endl;
    gSystem->Exit(1);
  }

  if (Enable::HFARBWD_MAGNETS_IP6)
  {
    hFarBwdBeamLine::enclosure_z_max = -4700.;
//    BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, hFarBwdBeamLine::starting_z);
    BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, hFarBwdBeamLine::enclosure_z_max);
    hFarBwdBeamLine::enclosure_r_max = 200.;
  }

  if (Enable::HFARBWD_MAGNETS_IP8)
  {
    hFarBwdBeamLine::enclosure_z_max = -4700.;
//    BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, hFarBwdBeamLine::starting_z);
    BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, hFarBwdBeamLine::enclosure_z_max);
    hFarBwdBeamLine::enclosure_r_max = 200.;
  }

  hFarBwdBeamLine::enclosure_center = 0.5 * (hFarBwdBeamLine::starting_z + hFarBwdBeamLine::enclosure_z_max);

//  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, hFarBwdBeamLine::enclosure_z_max);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, hFarBwdBeamLine::starting_z);

  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, hFarBwdBeamLine::enclosure_r_max);
}

void hFarBwdDefineMagnets(PHG4Reco *g4Reco)
{

//  cout << " HERE ???? " << endl;
//  gSystem->Exit(1);
//  exit(0);

  bool overlapCheck = Enable::OVERLAPCHECK || Enable::HFARBWD_OVERLAPCHECK;
  int verbosity = std::max(Enable::VERBOSITY, Enable::HFARBWD_VERBOSITY);

  hFarBwdBeamLine::hFarBwdBeamLineEnclosure = new PHG4CylinderSubsystem("hFarBwdBeamLineEnclosure");
//  hFarBwdBeamLine::hFarBwdBeamLineEnclosure->set_double_param("place_z", -2450);
  hFarBwdBeamLine::hFarBwdBeamLineEnclosure->set_double_param("place_z", hFarBwdBeamLine::enclosure_center);
  hFarBwdBeamLine::hFarBwdBeamLineEnclosure->set_double_param("radius", 0);
  hFarBwdBeamLine::hFarBwdBeamLineEnclosure->set_double_param("thickness", hFarBwdBeamLine::enclosure_r_max);  // This is intentionally made large 25cm radius
//  hFarBwdBeamLine::hFarBwdBeamLineEnclosure->set_double_param("length", 4000);
  hFarBwdBeamLine::hFarBwdBeamLineEnclosure->set_double_param("length", fabs(hFarBwdBeamLine::enclosure_z_max - hFarBwdBeamLine::starting_z));
  hFarBwdBeamLine::hFarBwdBeamLineEnclosure->set_string_param("material", "G4_Galactic");
  hFarBwdBeamLine::hFarBwdBeamLineEnclosure->set_color(.5, .5, .5, 0.2);
  hFarBwdBeamLine::hFarBwdBeamLineEnclosure->OverlapCheck(overlapCheck);
  if (verbosity)
    hFarBwdBeamLine::hFarBwdBeamLineEnclosure->Verbosity(verbosity);
  g4Reco->registerSubsystem(hFarBwdBeamLine::hFarBwdBeamLineEnclosure);

  string magFile;
  if (Enable::HFARFWD_MAGNETS_IP6)
    magFile = string(getenv("CALIBRATIONROOT")) + "/Beam/ip6_h_farFwdBeamLineMagnets.dat";
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
          cout << "coud not decode " << line << endl;
          gSystem->Exit(1);
        }
        else
        {
	  //------------------------
	  //Select only the magnet component in the far backward region
	  if (z > 0.0)
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
            bl = new BeamLineMagnetSubsystem("BWDBEAMLINEMAGNET", imagnet);
            bl->set_double_param("field_y", dipole_field_x);
            bl->set_double_param("fieldgradient", fieldgradient);
            bl->set_string_param("magtype", magtype);
            bl->set_double_param("length", length);
            bl->set_double_param("place_x", x);// relative position to mother vol.
            bl->set_double_param("place_y", y);// relative position to mother vol.
            bl->set_double_param("place_z", z - hFarBwdBeamLine::enclosure_center);// relative position to mother vol.
            bl->set_double_param("rot_y", angle);

            bl->set_double_param("field_global_position_x", x);// abs. position to world for field manager
            bl->set_double_param("field_global_position_y", y);// abs. position to world for field manager
            bl->set_double_param("field_global_position_z", z);// abs. position to world for field manager
            bl->set_double_param("field_global_rot_y", angle);// abs. rotation to world for field manager
            bl->set_double_param("inner_radius", inner_radius_zin);
            bl->set_double_param("outer_radius", outer_magnet_diameter / 2.);
            bl->SetActive(magnet_active);
            bl->BlackHole();
            bl->SetMotherSubsystem(hFarBwdBeamLine::hFarBwdBeamLineEnclosure);
            if (absorberactive)
            {
              bl->SetAbsorberActive();
            }
            bl->OverlapCheck(overlapCheck);
            bl->SuperDetector("BWDBEAMLINEMAGNET");
            if (verbosity)
              bl->Verbosity(verbosity);
            g4Reco->registerSubsystem(bl);

            // rag the B0 magnet
            if (imagnet == 0)
              hFarBwdBeamLine::B0Magnet = bl;
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

void hFarBwdDefineDetectorsIP6(PHG4Reco *g4Reco)
{
  bool overlapCheck = Enable::OVERLAPCHECK || Enable::HFARBWD_OVERLAPCHECK;
  if (Enable::HFARBWD_VIRTUAL_DETECTORS_IP6 && Enable::HFARBWD_VIRTUAL_DETECTORS_IP8)
  {
    cout << "You cannot have detectors enabled for both IP6 and IP8 ON at the same time" << endl;
    gSystem->Exit(1);
  }
 
   int verbosity = std::max(Enable::VERBOSITY, Enable::HFARBWD_VERBOSITY);

   auto *detBackward = new PHG4CylinderSubsystem("detBackward");
    detBackward->SuperDetector("backTruth");
    detBackward->set_double_param("place_x", 0);
    detBackward->set_double_param("place_y", 0);
//    detBackward->set_double_param("place_z", -500);
    detBackward->set_double_param("place_z", -500 - hFarBwdBeamLine::enclosure_center);
    detBackward->set_double_param("rot_y", 0);
    detBackward->set_double_param("radius", 0);
    detBackward->set_double_param("thickness", 30);  // This is intentionally made large 25cm radius
    detBackward->set_double_param("length", 0.03);
    detBackward->set_string_param("material", "G4_Si");

    detBackward->SetMotherSubsystem(hFarBwdBeamLine::hFarBwdBeamLineEnclosure);


    detBackward->SetActive();
    detBackward->OverlapCheck(overlapCheck);
    detBackward->set_color(1, 0, 0, 0.5);

    detBackward->BlackHole();
    if (verbosity) detBackward->Verbosity(verbosity);
    g4Reco->registerSubsystem(detBackward);

//  int verbosity = std::max(Enable::VERBOSITY, Enable::HFARFWD_VERBOSITY);

}

void hFarBwdDefineDetectorsIP8(PHG4Reco *g4Reco)
{
//  cout << __PRETTY_FUNCTION__ << " : IP8 setup is not yet validated!" << endl;
//  gSystem->Exit(1);

  bool overlapCheck = Enable::OVERLAPCHECK || Enable::HFARBWD_OVERLAPCHECK;
  if (Enable::HFARBWD_VIRTUAL_DETECTORS_IP6 && Enable::HFARBWD_VIRTUAL_DETECTORS_IP8)
  {
    cout << "You cannot have detectors enabled for both IP6 and IP8 ON at the same time" << endl;
    gSystem->Exit(1);
  }

   int verbosity = std::max(Enable::VERBOSITY, Enable::HFARBWD_VERBOSITY);

   auto *detBackward = new PHG4CylinderSubsystem("detBackward");
    detBackward->SuperDetector("backTruth");
    detBackward->set_double_param("place_x", 0);
    detBackward->set_double_param("place_y", 0);
//    detBackward->set_double_param("place_z", -500);
    detBackward->set_double_param("place_z", -500 - hFarBwdBeamLine::enclosure_center);
    detBackward->set_double_param("rot_y", 0);
    detBackward->set_double_param("radius", 0);
    detBackward->set_double_param("thickness", 30);  // This is intentionally made large 25cm radius
    detBackward->set_double_param("length", 0.03);
    detBackward->set_string_param("material", "G4_Si");

    detBackward->SetMotherSubsystem(hFarBwdBeamLine::hFarBwdBeamLineEnclosure);

    detBackward->SetActive();
    detBackward->OverlapCheck(overlapCheck);
    detBackward->set_color(1, 0, 0, 0.5);

    detBackward->BlackHole();
    if (verbosity) detBackward->Verbosity(verbosity);
    g4Reco->registerSubsystem(detBackward);

//  int verbosity = std::max(Enable::VERBOSITY, Enable::HFARFWD_VERBOSITY);

}

void hFarBwdDefineBeamPipe(PHG4Reco *g4Reco)
{
//  int verbosity = std::max(Enable::VERBOSITY, Enable::HFARFWD_VERBOSITY);
}





#endif
