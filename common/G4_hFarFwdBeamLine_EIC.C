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

#include <eicg4b0/EICG4B0Subsystem.h>
#include <eicg4b0ecal/EICG4B0ECALSubsystem.h>
#include <eicg4rp/EICG4RPSubsystem.h>

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
	
// Detector configuration options
  bool ZDC_DISABLE_BLACKHOLE = false;
  bool B0_DISABLE_HITPLANE = false;
  bool B0_FULLHITPLANE = false;
  bool B0_VAR_PIPE_HOLE = false;
  bool B0_CIRCLE_PIPE_HOLE = false;
  bool RP_DISABLE_HITPLANE = false;
  bool B0ECALTOWERS = true; //Set to 'false' for nice PackMan views. Set 'true' for physics studies.

  //enabled automatically in hFarFwdBeamLineInit(), unless overridden by user
  bool HFARFWD_MAGNETS_IP6 = false;
  bool HFARFWD_MAGNETS_IP8 = false;

  //enabled automatically in hFarFwdBeamLineInit(), unless overridden by user
  bool HFARFWD_VIRTUAL_DETECTORS_IP6 = false;
  bool HFARFWD_VIRTUAL_DETECTORS_IP8 = false;

  bool FFR_EVAL = false;

}  // namespace Enable

namespace hFarFwdBeamLine
{
  double starting_z = 500;  //cm as center-forward interface
  double enclosure_z_max = NAN;
  double enclosure_r_max = NAN;
  double enclosure_center = NAN;

  PHG4CylinderSubsystem *hFarFwdBeamLineEnclosure(nullptr);

  BeamLineMagnetSubsystem *B0Magnet = (nullptr);
  double B0Magnet_x = NAN;
  double B0Magnet_y = NAN;
  double B0Magnet_z = NAN;
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
  hFarFwdBeamLine::hFarFwdBeamLineEnclosure->SetActive();
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
  bool magnet_active = true;
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
             fieldgradient = fieldgradient * scaleFactor;
   	  }

          if (magnetlist.empty() || magnetlist.find(imagnet) != magnetlist.end())
          {
            bl = new BeamLineMagnetSubsystem("BEAMLINEMAGNET", imagnet);
            bl->set_double_param("field_y", MagFieldFlip(dipole_field_x));
            bl->set_double_param("fieldgradient", MagFieldFlip(fieldgradient));
            bl->set_string_param("magtype", magtype);
            bl->set_double_param("length", length);
            bl->set_double_param("place_x", PosFlip(x));				// relative position to mother vol.
            bl->set_double_param("place_y", y);						// relative position to mother vol.
            bl->set_double_param("place_z", z - hFarFwdBeamLine::enclosure_center);	// relative position to mother vol.
            bl->set_double_param("field_global_position_x", PosFlip(x));		// abs. position to world for field manager
            bl->set_double_param("field_global_position_y", y);				// abs. position to world for field manager
            bl->set_double_param("field_global_position_z", z);				// abs. position to world for field manager
            bl->set_double_param("rot_y", AngleFlip(angle));
            bl->set_double_param("field_global_rot_y", AngleFlip(angle));		// abs. rotation to world for field manager
            bl->set_double_param("inner_radius", inner_radius_zin);
            bl->set_double_param("outer_radius", outer_magnet_diameter / 2.);
            bl->SetActive(magnet_active);
            bl->SetAbsorberActive();
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
	    {	//To tell the B0 Calorimeter the global coordinates of the B0 Magnet
            	hFarFwdBeamLine::B0Magnet = bl;
		hFarFwdBeamLine::B0Magnet_x = PosFlip(x);
		hFarFwdBeamLine::B0Magnet_y = y;
		hFarFwdBeamLine::B0Magnet_z = z;
		}
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
  
  
  
  //----------------------
  // Roman Pots
  //----------------------
  
  if( ! Enable::RP_DISABLE_HITPLANE )
  {
	  string paramFile = string(getenv("CALIBRATIONROOT")) + "/RomanPots/RP_parameters_IP6.dat";
	  int Nlayers = GetParameterFromFile <int> (paramFile, "Number_layers");

	  for( int layer = 0; layer < Nlayers; layer++ ) {
		  auto *detRP = new EICG4RPSubsystem(Form("rpTruth_%d", layer), layer);
		  detRP->SuperDetector("rpTruth");
		  detRP->SetParameterFile( paramFile );
		  detRP->set_double_param("FFenclosure_center", hFarFwdBeamLine::enclosure_center );
		  detRP->set_int_param("layerNumber", layer + 1);    

		  detRP->OverlapCheck(overlapCheck);
		  detRP->SetMotherSubsystem(hFarFwdBeamLine::hFarFwdBeamLineEnclosure);
		  detRP->SetActive(true);
		  if (verbosity) detRP->Verbosity(verbosity);
		  g4Reco->registerSubsystem(detRP);
	  }
  }

 
   //---------------------------------
   // B0 implementation
   // Three choices: 1. Realistic detector; 2. Circulat plane; 3. hit plane with realistic detector goemetry
	
        double b0tr_z = 0; //Subsystem position relative to B0 magnet (for iterator)
        const int b0DetNr = 4;
        const double b0Mag_zCent = 640;
        const double b0Mag_zLen = 120;
	const double b0tr[4]={10,40,70,100};
//	const double b0tr[4]={10,45,80,115}; //Tracker layers when no ECAL
        const double b0Cu_zLen = .2; //B0 dead material length
        const double b0Si_zLen = .1; //B0 Si length
        const double b0Ecal_zLen = 10; //B0 Ecal length
        double pipe_hole_r = 3.5; //detector cut off for beam pipe
	double pipe_hole = 2.5;
	const double cable_hole = 2.0;
	const double cable_x = -17.0;
        double pipe_x = -1.; //pipe hole position
        const double d_radius = 7.0; //detector cut off Packman
        const double b0_radius = 19.0; //outer radius of B0-detector
        const double b0_magradius = 20.0; //inner radius of B0-magnet
        const double spanning_angle = 240; //spanning angle Packman
        const double b0Ecal_z = 48;//B0 ECal position (relative to the B0-magnet)
        double start_angle = 60; //start angle Packman
	const double cross_angle = 0.025;

    if (Enable::B0_DISABLE_HITPLANE) {

	// Choice 1 realistic detector
//	const double b0tr[4]={10,45,80,115};
	//const double b0tr[4]={0,30,60,90};
	//const double b0tr[5]={0,25,50,75,100};
	cout << "Realistic B0"<<endl;
        for (int i = 0; i < b0DetNr; i++)
        {
	if (Enable::B0_VAR_PIPE_HOLE){
		pipe_hole = b0tr[i]*cross_angle;
		pipe_x = - cross_angle*(b0Mag_zCent - b0Mag_zLen/2 + b0tr[i]/2) - hFarFwdBeamLine::B0Magnet_x;
	}
	else if (Enable::B0_CIRCLE_PIPE_HOLE){
		pipe_hole = 0.1;
		pipe_hole_r = pipe_hole_r + b0tr[b0DetNr-1]*cross_angle/2;
		pipe_x = - cross_angle*(b0Mag_zCent - b0Mag_zLen/2 + b0tr[b0DetNr-1]/2) - hFarFwdBeamLine::B0Magnet_x;
	}
	else {
		pipe_hole = b0tr[b0DetNr-1]*cross_angle;
		pipe_x = - cross_angle*(b0Mag_zCent - b0Mag_zLen/2 + b0tr[b0DetNr-1]/2) - hFarFwdBeamLine::B0Magnet_x;
	}
	cout <<"Starting B0 Tracker layer "<<i+1<<endl;
	cout <<"Pipe Hole: "<< pipe_hole<<"\t"<<pipe_x<<endl;
	  b0tr_z = b0tr[i] - b0Mag_zLen / 2;
          auto *detB0 = new EICG4B0Subsystem(Form("b0Truth_%d", i), i);
          detB0->SuperDetector(Form("b0Truth_%d", i));
          detB0->set_double_param("place_x", 0);
          detB0->set_double_param("place_y", 0);
      //  detB0->set_int_param("ispipe", 0); //for future pipe implementation
          detB0->set_double_param("pipe_hole", pipe_hole);
          detB0->set_double_param("cable_hole", cable_hole);
          detB0->set_double_param("outer_radius", b0_radius);
          detB0->set_double_param("d_radius", d_radius);
          detB0->set_double_param("length", b0Si_zLen);
          detB0->set_string_param("material", "G4_Si");
          detB0->set_double_param("startAngle",start_angle);
          detB0->set_double_param("spanningAngle",spanning_angle);
          detB0->set_double_param("detid",i);
          detB0->set_double_param("pipe_x", pipe_x);
          detB0->set_double_param("pipe_y", 0);
          detB0->set_double_param("pipe_z", 0);
          detB0->set_double_param("pipe_hole_r", pipe_hole_r);
          detB0->set_double_param("cable_x", cable_x);
          detB0->set_double_param("cable_y", 0);
          detB0->set_double_param("cable_z", 0);
          detB0->set_double_param("place_z", b0tr_z);  // relative to B0 magnet
          detB0->SetActive(true);
          if (verbosity)
          detB0->Verbosity(verbosity);
          detB0->OverlapCheck(overlapCheck);
          detB0->SetMotherSubsystem(hFarFwdBeamLine::B0Magnet);
          g4Reco->registerSubsystem(detB0);
// For B0 Tracking Implementation
          if (Enable::B0TRACKING){ 
		  if (B0TRACKING::FastKalmanFilter)
   		  {
   	 		B0TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + Form("b0Truth_%d", i) ,           //      const std::string& phg4hitsNames,
        	                             B0TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             G4B0TRACKING::PositionResolution,           //      const float radres,
                                             G4B0TRACKING::PositionResolution,           //      const float phires,
                                             0,              //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
	   		B0TRACKING::FastKalmanFilter->add_zplane_state(Form("b0Truth_%d", i), b0Mag_zCent+b0tr_z);
   		 	B0TRACKING::FastKalmanFilterB0Track->add_phg4hits(string("G4HIT_") + Form("b0Truth_%d", i) ,           //      const std::string& phg4hitsNames,
        	                             B0TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             G4B0TRACKING::PositionResolution,           //      const float radres,
                                             G4B0TRACKING::PositionResolution,           //      const float phires,
                                             0,              //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
   			 B0TRACKING::FastKalmanFilterB0Track->add_zplane_state(Form("b0Truth_%d", i), b0Mag_zCent+b0tr_z);
	   		 B0TRACKING::B0ProjectionNames.insert(Form("b0Truth_%d", i));
 		 }
	  }

          auto *detB0e = new EICG4B0Subsystem(Form("b0Dead_%d", i), i);
          detB0e->SuperDetector("b0Dead");
      //  detB0e->set_int_param("ispipe", 0); //for future pipe implementation
          detB0e->set_double_param("pipe_hole", pipe_hole);
          detB0e->set_double_param("place_x", 0);
          detB0e->set_double_param("place_y", 0);
          detB0e->set_double_param("d_radius", d_radius);
          detB0e->set_double_param("pipe_x", pipe_x);
          detB0e->set_double_param("pipe_y", 0);
          detB0e->set_double_param("pipe_z", 0);
          detB0e->set_double_param("pipe_hole_r", pipe_hole_r);
          detB0e->set_double_param("cable_x", cable_x);
          detB0e->set_double_param("cable_y", 0);
          detB0e->set_double_param("cable_z", 0);
          detB0e->set_double_param("outer_radius", b0_radius);
          detB0e->set_double_param("length", b0Cu_zLen);
          detB0e->set_string_param("material", "G4_Cu");
          detB0e->set_double_param("detid",i);
          detB0e->set_double_param("startAngle",start_angle);
          detB0e->set_double_param("spanningAngle",spanning_angle);
          detB0e->set_double_param("place_z", b0tr_z +(b0Cu_zLen+b0Si_zLen)/2) ;  // relative to B0 magnet
          detB0e->SetActive(false);
          if (verbosity)
            detB0e->Verbosity(verbosity);
          detB0e->OverlapCheck(overlapCheck);
          detB0e->SetMotherSubsystem(hFarFwdBeamLine::B0Magnet);
          g4Reco->registerSubsystem(detB0e);
        }
      
  if (Enable::B0ECAL) {
	pipe_hole = b0Mag_zLen*cross_angle;
	pipe_x = - cross_angle*b0Mag_zCent - hFarFwdBeamLine::B0Magnet_x;
	if (Enable::B0_CIRCLE_PIPE_HOLE){
		pipe_hole = 0.1;
		pipe_hole_r = pipe_hole_r + b0Mag_zLen*cross_angle/2;
	}
	cout <<"Starting B0 ECAL "<<endl;
	cout <<"Pipe Hole: "<< pipe_hole<<"\t"<<pipe_x<<endl;
        if (Enable::B0ECALTOWERS){				//Use this option to do physics studies
//	pipe_x=-1.25;
//	pipe_hole=3.0;
	cout << hFarFwdBeamLine::B0Magnet_x<<endl;
	  	ostringstream mapping_b0ecal;
		mapping_b0ecal << getenv("CALIBRATIONROOT") << "/B0Ecal/mapping/B0ECAL_mapping_v2.txt"; // Specify the mapping file for B0 ECal Towers here
	//	mapping_b0ecal << "B0ECAL_mapping_v2.txt"; // Specify the mapping file for B0 ECal Towers here
		//cout <<"Will use B0 mapping file "<< mapping_b0ecal.str()<<endl;  
	        auto *B0Ecal = new EICG4B0ECALSubsystem("B0ECAL");
		B0Ecal->SetTowerMappingFile(mapping_b0ecal.str());
	        B0Ecal->SuperDetector("B0ECAL");
	        B0Ecal->set_double_param("pipe_hole", pipe_hole);
       		B0Ecal->set_double_param("place_x", 0);
	        B0Ecal->set_double_param("place_y", 0);
        	B0Ecal->set_double_param("place_z", b0Ecal_z);
	        B0Ecal->set_double_param("pipe_x", pipe_x);
        	B0Ecal->set_double_param("pipe_y", 0);
	        B0Ecal->set_double_param("pipe_z", 0);
         	B0Ecal->set_double_param("pipe_hole_r", pipe_hole_r);
          	B0Ecal->set_double_param("cable_x", cable_x);
          	B0Ecal->set_double_param("cable_y", 0);
         	B0Ecal->set_double_param("cable_z", 0);
	        B0Ecal->set_double_param("length", b0Ecal_zLen);
	        B0Ecal->set_double_param("outer_radius", b0_radius);
	        B0Ecal->set_double_param("d_radius", d_radius);
	        B0Ecal->set_string_param("material", "G4_PbWO4");
	        B0Ecal->set_double_param("startAngle",start_angle);
	        B0Ecal->set_double_param("spanningAngle",spanning_angle);
	        B0Ecal->set_double_param("detid",0);
	        B0Ecal->set_double_param("global_x",hFarFwdBeamLine::B0Magnet_x);
	        B0Ecal->set_double_param("global_y",hFarFwdBeamLine::B0Magnet_y);
	        B0Ecal->set_double_param("global_z",hFarFwdBeamLine::B0Magnet_z);
		B0Ecal->set_int_param("lightyield",1); 		//Note additional parameter for storing Light Yield in B0 Ecal
		B0Ecal->SetActive(true);
	        if (verbosity)
        	  B0Ecal->Verbosity(verbosity);
	        B0Ecal->OverlapCheck(overlapCheck);
        	B0Ecal->SetMotherSubsystem(hFarFwdBeamLine::B0Magnet);
       		g4Reco->registerSubsystem(B0Ecal);
	}
	else {					//Use this option to have a circular packman-shape of the B0 ECal for plots.
	        auto *B0Ecal = new EICG4B0Subsystem(Form("b0Truth_%d", 2*b0DetNr), 2*b0DetNr);
        	B0Ecal->SuperDetector("b0Truth");
	        B0Ecal->set_double_param("pipe_hole", pipe_hole);
       		B0Ecal->set_double_param("place_x", 0);
	        B0Ecal->set_double_param("place_y", 0);
        	B0Ecal->set_double_param("place_z", b0Ecal_z);
	        B0Ecal->set_double_param("pipe_x", pipe_x);
        	B0Ecal->set_double_param("pipe_y", 0);
	        B0Ecal->set_double_param("pipe_z", 0);
         	B0Ecal->set_double_param("pipe_hole_r", pipe_hole_r);
          	B0Ecal->set_double_param("cable_x", cable_x);
          	B0Ecal->set_double_param("cable_y", 0);
         	B0Ecal->set_double_param("cable_z", 0);
	        B0Ecal->set_double_param("length", b0Ecal_zLen);
	        B0Ecal->set_double_param("outer_radius", b0_radius);
	        B0Ecal->set_double_param("d_radius", d_radius);
	        B0Ecal->set_string_param("material", "G4_PbWO4");
	        B0Ecal->set_double_param("startAngle",start_angle);
	        B0Ecal->set_double_param("spanningAngle",spanning_angle);
	        B0Ecal->set_double_param("detid",2*b0DetNr);
		B0Ecal->SetActive(true);
        	if (verbosity)
		  B0Ecal->Verbosity(verbosity);
	        B0Ecal->OverlapCheck(overlapCheck);
	        B0Ecal->SetMotherSubsystem(hFarFwdBeamLine::B0Magnet);
        	g4Reco->registerSubsystem(B0Ecal);
	}
      
        auto *B0Ecale = new EICG4B0Subsystem(Form("b0Dead_%d", b0DetNr), b0DetNr); //B0 ECal dead layer is the same subsystem as other four dead layers
        B0Ecale->SuperDetector("b0Dead");
      //  B0Ecale->set_int_param("ispipe", 0); //for future pipe implementation
        B0Ecale->set_double_param("pipe_hole", pipe_hole);
        B0Ecale->set_double_param("place_x", 0);
        B0Ecale->set_double_param("place_y", 0);
        B0Ecale->set_double_param("place_z", b0Ecal_z + (b0Ecal_zLen + b0Cu_zLen)/2);
        B0Ecale->set_double_param("pipe_x", pipe_x);
        B0Ecale->set_double_param("pipe_y", 0);
        B0Ecale->set_double_param("pipe_z", 0);
          B0Ecale->set_double_param("pipe_hole_r", pipe_hole_r);
          B0Ecale->set_double_param("cable_x", cable_x);
          B0Ecale->set_double_param("cable_y", 0);
          B0Ecale->set_double_param("cable_z", 0);
        B0Ecale->set_double_param("length", b0Cu_zLen);
        B0Ecale->set_double_param("d_radius", d_radius);
        B0Ecale->set_double_param("outer_radius", b0_radius);
        B0Ecale->set_string_param("material", "G4_Cu");
        B0Ecale->set_double_param("startAngle",start_angle);
        B0Ecale->set_double_param("spanningAngle",spanning_angle);
        B0Ecale->set_double_param("detid",b0DetNr+1);
        //B0Ecale->SetActive(true);
        B0Ecale->SetActive(false);
        if (verbosity)
          B0Ecale->Verbosity(verbosity);
        B0Ecale->OverlapCheck(overlapCheck);
        B0Ecale->SetMotherSubsystem(hFarFwdBeamLine::B0Magnet);
        g4Reco->registerSubsystem(B0Ecale);
	}
    } else {

       if (Enable::B0_FULLHITPLANE) {

	// Choice 2 circular hit planes
	cout << "Circular hit planes"<<endl;

       	    for (int i = 0; i < b0DetNr; i++)
       	    {
	      b0tr_z = b0tr[i] - b0Mag_zLen / 2;
       	      auto *detB0 = new PHG4CylinderSubsystem(Form("b0Truth_%d", i), i);
              detB0->SuperDetector("b0Truth");
              //detB0->SuperDetector(Form("b0Truth_%d", i));
       	      detB0->set_double_param("radius", 0);
       	      detB0->set_double_param("thickness", 20);
       	      detB0->set_double_param("length", 0.1);
       	      detB0->set_string_param("material", "G4_Si");
       	      detB0->set_double_param("place_z", b0tr_z);  // relative to B0 magnet
       	      detB0->SetActive(true);
       	      if (verbosity) detB0->Verbosity(verbosity);
       	      detB0->OverlapCheck(overlapCheck);
       	    
       	      detB0->SetMotherSubsystem(hFarFwdBeamLine::B0Magnet);
       	    
       	      g4Reco->registerSubsystem(detB0);
          if (Enable::B0TRACKING){ 
		  if (B0TRACKING::FastKalmanFilter)
   		  {
   	 		B0TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + Form("b0Truth_%d", i) ,           //      const std::string& phg4hitsNames,
        	                             B0TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             G4B0TRACKING::PositionResolution,           //      const float radres,
                                             G4B0TRACKING::PositionResolution,           //      const float phires,
                                             0,              //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
	   		B0TRACKING::FastKalmanFilter->add_zplane_state(Form("b0Truth_%d", i), b0Mag_zCent+b0tr_z);
   		 	B0TRACKING::FastKalmanFilterB0Track->add_phg4hits(string("G4HIT_") + Form("b0Truth_%d", i) ,           //      const std::string& phg4hitsNames,
        	                             B0TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             G4B0TRACKING::PositionResolution,           //      const float radres,
                                             G4B0TRACKING::PositionResolution,           //      const float phires,
                                             0,              //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
   			 B0TRACKING::FastKalmanFilterB0Track->add_zplane_state(Form("b0Truth_%d", i), b0Mag_zCent+b0tr_z);
	   		 B0TRACKING::B0ProjectionNames.insert(Form("b0Truth_%d", i));
 		 }
	  }
       	    
       	    }

	} else {

        /// Fun4All default B0 planes
	/// Choice 3 Hit planes with real detector geometry
	cout << "Realistic hit planes"<<endl;
	    
	    for (int i = 0; i < b0DetNr; i++) {
	if (Enable::B0_VAR_PIPE_HOLE){
		pipe_hole = b0tr[i]*cross_angle;
		pipe_x = - cross_angle*(b0Mag_zCent - b0Mag_zLen/2 + b0tr[i]/2) - hFarFwdBeamLine::B0Magnet_x;
	}
	else if (Enable::B0_CIRCLE_PIPE_HOLE){
		pipe_hole = 0.1;
		pipe_hole_r = pipe_hole_r + b0tr[b0DetNr-1]*cross_angle/2;
		pipe_x = - cross_angle*(b0Mag_zCent - b0Mag_zLen/2 + b0tr[b0DetNr-1]/2) - hFarFwdBeamLine::B0Magnet_x;
	}
	else {
		pipe_hole = b0tr[b0DetNr-1]*cross_angle;
		pipe_x = - cross_angle*(b0Mag_zCent - b0Mag_zLen/2 + b0tr[b0DetNr-1]/2) - hFarFwdBeamLine::B0Magnet_x;
	}
	cout <<"Starting B0 Tracker layer "<<i+1<<endl;
	cout <<"Pipe Hole: "<< pipe_hole<<"\t"<<pipe_x<<endl;
	  b0tr_z = b0tr[i] - b0Mag_zLen / 2;
          auto *detB0 = new EICG4B0Subsystem(Form("b0Truth_%d", i), i);
          detB0->SuperDetector(Form("b0Truth_%d", i));
          detB0->set_double_param("place_x", 0);
          detB0->set_double_param("place_y", 0);
      //  detB0->set_int_param("ispipe", 0); //for future pipe implementation
          detB0->set_double_param("pipe_hole", pipe_hole);
          detB0->set_double_param("cable_hole", cable_hole);
          detB0->set_double_param("outer_radius", b0_radius);
          detB0->set_double_param("d_radius", d_radius);
          detB0->set_double_param("length", b0Si_zLen);
          detB0->set_string_param("material", "G4_Si");
          detB0->set_double_param("startAngle",start_angle);
          detB0->set_double_param("spanningAngle",spanning_angle);
          detB0->set_double_param("detid",i);
          detB0->set_double_param("pipe_x", pipe_x);
          detB0->set_double_param("pipe_y", 0);
          detB0->set_double_param("pipe_z", 0);
          detB0->set_double_param("pipe_hole_r", pipe_hole_r);
          detB0->set_double_param("cable_x", cable_x);
          detB0->set_double_param("cable_y", 0);
          detB0->set_double_param("cable_z", 0);
          detB0->set_double_param("place_z", b0tr_z);  // relative to B0 magnet
	        detB0->SetActive(true);
	        if (verbosity)
	          detB0->Verbosity(verbosity);
	        detB0->OverlapCheck(overlapCheck);
	        detB0->SetMotherSubsystem(hFarFwdBeamLine::B0Magnet);
	        g4Reco->registerSubsystem(detB0);
          if (Enable::B0TRACKING){ 
		  if (B0TRACKING::FastKalmanFilter)
   		  {
   	 		B0TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + Form("b0Truth_%d", i) ,           //      const std::string& phg4hitsNames,
        	                             B0TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             G4B0TRACKING::PositionResolution,           //      const float radres,
                                             G4B0TRACKING::PositionResolution,           //      const float phires,
                                             0,              //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
	   		B0TRACKING::FastKalmanFilter->add_zplane_state(Form("b0Truth_%d", i), b0Mag_zCent+b0tr_z);
   		 	B0TRACKING::FastKalmanFilterB0Track->add_phg4hits(string("G4HIT_") + Form("b0Truth_%d", i) ,           //      const std::string& phg4hitsNames,
        	                             B0TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             G4B0TRACKING::PositionResolution,           //      const float radres,
                                             G4B0TRACKING::PositionResolution,           //      const float phires,
                                             0,              //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
   			 B0TRACKING::FastKalmanFilterB0Track->add_zplane_state(Form("b0Truth_%d", i), b0Mag_zCent+b0tr_z);
	   		 B0TRACKING::B0ProjectionNames.insert(Form("b0Truth_%d", i));
 		 }
	  }
	    }

	} 
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
  if (verbosity) detZDCsurrogate->Verbosity(verbosity);
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


  //----------------------
  // Roman Pots: Both sets before and near the secondary focus
  //----------------------

  if( ! Enable::RP_DISABLE_HITPLANE )
  {
	  string paramFile = string(getenv("CALIBRATIONROOT")) + "/RomanPots/RP_parameters_IP8.dat";
	  int Nlayers = GetParameterFromFile <int> (paramFile, "Number_layers");

	  for( int layer = 0; layer < Nlayers; layer++ ) {
		  auto *detRP = new EICG4RPSubsystem(Form("rpTruth_%d", layer), layer);
		  detRP->SuperDetector("rpTruth");
		  detRP->SetParameterFile( paramFile );
		  detRP->set_double_param("FFenclosure_center", hFarFwdBeamLine::enclosure_center );
		  detRP->set_int_param("layerNumber", layer + 1);    

		  detRP->OverlapCheck(overlapCheck);
		  detRP->SetMotherSubsystem(hFarFwdBeamLine::hFarFwdBeamLineEnclosure);
		  detRP->SetActive(true);
		  if (verbosity) detRP->Verbosity(verbosity);
		  g4Reco->registerSubsystem(detRP);
	  }
  }

  if (verbosity > 0)
  {
    std::cout << "B0Magnet can be mother = " << hFarFwdBeamLine::B0Magnet->CanBeMotherSubsystem() << std::endl;
  }

/*  const int b0DetNr = 4;
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
  }*/
   //---------------------------------
   // B0 implementation
   // Three choices: 1. Realistic detector; 2. Circulat plane; 3. hit plane with realistic detector goemetry
        double b0tr_z = 0; //Subsystem position relative to B0 magnet (for iterator)
        const int b0DetNr = 4;
        const double b0Mag_zCent = 610;
        const double b0Mag_zLen = 120;
	const double b0tr[4]={10,40,70,100};
        const double b0Cu_zLen = .2; //B0 dead material length
        const double b0Si_zLen = .1; //B0 Si length
        const double b0Ecal_zLen = 10; //B0 Ecal length
        double pipe_hole_r = 3.5; //detector cut off for beam pipe
	double pipe_hole = 2.5;
	const double cable_hole = 2.0;
	const double cable_x = 21.5;
        double pipe_x = -1.; //pipe hole position
        const double d_radius = 7.0; //detector cut off Packman
        const double b0_radius = 23.5; //outer radius of B0-detector
        const double b0_magradius = 24.5; //inner radius of B0-magnet
        const double spanning_angle = 240; //spanning angle Packman
        const double b0Ecal_z = 48;//B0 ECal position (relative to the B0-magnet)
        double start_angle = -120; //start angle Packman
	const double cross_angle = 0.035;

    if (Enable::B0_DISABLE_HITPLANE) {

	// Choice 1 realistic detector
//	const double b0tr[4]={10,45,80,115};
	//const double b0tr[4]={0,30,60,90};
	//const double b0tr[5]={0,25,50,75,100};
	cout << "Realistic B0"<<endl;
        for (int i = 0; i < b0DetNr; i++)
        {
	if (Enable::B0_VAR_PIPE_HOLE){
		pipe_hole = b0tr[i]*cross_angle;
		pipe_x =  cross_angle*(b0Mag_zCent - b0Mag_zLen/2 + b0tr[i]/2) - hFarFwdBeamLine::B0Magnet_x;
	}
	else if (Enable::B0_CIRCLE_PIPE_HOLE){
		pipe_hole = 0.1;
		pipe_hole_r = pipe_hole_r + b0tr[b0DetNr-1]*cross_angle/2;
		pipe_x =  cross_angle*(b0Mag_zCent - b0Mag_zLen/2 + b0tr[b0DetNr-1]/2) - hFarFwdBeamLine::B0Magnet_x;
	}
	else {
		pipe_hole = b0tr[b0DetNr-1]*cross_angle;
		pipe_x =  cross_angle*(b0Mag_zCent - b0Mag_zLen/2 + b0tr[b0DetNr-1]/2) - hFarFwdBeamLine::B0Magnet_x;
	}
	cout <<"Starting B0 Tracker layer "<<i+1<<endl;
	cout <<"Pipe Hole: "<< pipe_hole<<"\t"<<pipe_x<<endl;
	  b0tr_z = b0tr[i] - b0Mag_zLen / 2;
          auto *detB0 = new EICG4B0Subsystem(Form("b0Truth_%d", i), i);
          detB0->SuperDetector(Form("b0Truth_%d", i));
          detB0->set_double_param("place_x", 0);
          detB0->set_double_param("place_y", 0);
      //  detB0->set_int_param("ispipe", 0); //for future pipe implementation
          detB0->set_double_param("pipe_hole", pipe_hole);
          detB0->set_double_param("cable_hole", cable_hole);
          detB0->set_double_param("outer_radius", b0_radius);
          detB0->set_double_param("d_radius", d_radius);
          detB0->set_double_param("length", b0Si_zLen);
          detB0->set_string_param("material", "G4_Si");
          detB0->set_double_param("startAngle",start_angle);
          detB0->set_double_param("spanningAngle",spanning_angle);
          detB0->set_double_param("detid",i);
          detB0->set_double_param("pipe_x", pipe_x);
          detB0->set_double_param("pipe_y", 0);
          detB0->set_double_param("pipe_z", 0);
          detB0->set_double_param("pipe_hole_r", pipe_hole_r);
          detB0->set_double_param("cable_x", cable_x);
          detB0->set_double_param("cable_y", 0);
          detB0->set_double_param("cable_z", 0);
          detB0->set_double_param("place_z", b0tr_z);  // relative to B0 magnet
          detB0->SetActive(true);
          if (verbosity)
          detB0->Verbosity(verbosity);
          detB0->OverlapCheck(overlapCheck);
          detB0->SetMotherSubsystem(hFarFwdBeamLine::B0Magnet);
          g4Reco->registerSubsystem(detB0);
// For B0 Tracking Implementation
          if (Enable::B0TRACKING){ 
		  if (B0TRACKING::FastKalmanFilter)
   		  {
   	 		B0TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + Form("b0Truth_%d", i) ,           //      const std::string& phg4hitsNames,
        	                             B0TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             G4B0TRACKING::PositionResolution,           //      const float radres,
                                             G4B0TRACKING::PositionResolution,           //      const float phires,
                                             0,              //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
	   		B0TRACKING::FastKalmanFilter->add_zplane_state(Form("b0Truth_%d", i), b0Mag_zCent+b0tr_z);
   		 	B0TRACKING::FastKalmanFilterB0Track->add_phg4hits(string("G4HIT_") + Form("b0Truth_%d", i) ,           //      const std::string& phg4hitsNames,
        	                             B0TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             G4B0TRACKING::PositionResolution,           //      const float radres,
                                             G4B0TRACKING::PositionResolution,           //      const float phires,
                                             0,              //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
   			 B0TRACKING::FastKalmanFilterB0Track->add_zplane_state(Form("b0Truth_%d", i), b0Mag_zCent+b0tr_z);
	   		 B0TRACKING::B0ProjectionNames.insert(Form("b0Truth_%d", i));
 		 }
	  }

          auto *detB0e = new EICG4B0Subsystem(Form("b0Dead_%d", i), i);
          detB0e->SuperDetector("b0Dead");
      //  detB0e->set_int_param("ispipe", 0); //for future pipe implementation
          detB0e->set_double_param("pipe_hole", pipe_hole);
          detB0e->set_double_param("place_x", 0);
          detB0e->set_double_param("place_y", 0);
          detB0e->set_double_param("d_radius", d_radius);
          detB0e->set_double_param("pipe_x", pipe_x);
          detB0e->set_double_param("pipe_y", 0);
          detB0e->set_double_param("pipe_z", 0);
          detB0e->set_double_param("pipe_hole_r", pipe_hole_r);
          detB0e->set_double_param("cable_x", cable_x);
          detB0e->set_double_param("cable_y", 0);
          detB0e->set_double_param("cable_z", 0);
          detB0e->set_double_param("outer_radius", b0_radius);
          detB0e->set_double_param("length", b0Cu_zLen);
          detB0e->set_string_param("material", "G4_Cu");
          detB0e->set_double_param("detid",i);
          detB0e->set_double_param("startAngle",start_angle);
          detB0e->set_double_param("spanningAngle",spanning_angle);
          detB0e->set_double_param("place_z", b0tr_z +(b0Cu_zLen+b0Si_zLen)/2) ;  // relative to B0 magnet
          detB0e->SetActive(false);
          if (verbosity)
            detB0e->Verbosity(verbosity);
          detB0e->OverlapCheck(overlapCheck);
          detB0e->SetMotherSubsystem(hFarFwdBeamLine::B0Magnet);
          g4Reco->registerSubsystem(detB0e);
        }
      
  if (Enable::B0ECAL) {
	pipe_hole = b0Mag_zLen*cross_angle;
	pipe_x =  cross_angle*b0Mag_zCent - hFarFwdBeamLine::B0Magnet_x;
	if (Enable::B0_CIRCLE_PIPE_HOLE){
		pipe_hole = 0.1;
		pipe_hole_r = pipe_hole_r + b0Mag_zLen*cross_angle/2;
	}
	cout <<"Starting B0 ECAL "<<endl;
	cout <<"Pipe Hole: "<< pipe_hole<<"\t"<<pipe_x<<endl;
        if (Enable::B0ECALTOWERS){				//Use this option to do physics studies
//	pipe_x=-1.25;
//	pipe_hole=3.0;
	cout << hFarFwdBeamLine::B0Magnet_x<<endl;
	  	ostringstream mapping_b0ecal;
		mapping_b0ecal << getenv("CALIBRATIONROOT") << "/B0Ecal/mapping/B0ECAL_mapping_ip8_v1.txt"; // Specify the mapping file for B0 ECal Towers here
//		mapping_b0ecal << "B0ECAL_mapping_ip8_v1.txt"; // Specify the mapping file for B0 ECal Towers here
		//cout <<"Will use B0 mapping file "<< mapping_b0ecal.str()<<endl;  
	        auto *B0Ecal = new EICG4B0ECALSubsystem("B0ECAL");
		B0Ecal->SetTowerMappingFile(mapping_b0ecal.str());
	        B0Ecal->SuperDetector("B0ECAL");
	        B0Ecal->set_double_param("pipe_hole", pipe_hole);
       		B0Ecal->set_double_param("place_x", 0);
	        B0Ecal->set_double_param("place_y", 0);
        	B0Ecal->set_double_param("place_z", b0Ecal_z);
	        B0Ecal->set_double_param("pipe_x", pipe_x);
        	B0Ecal->set_double_param("pipe_y", 0);
	        B0Ecal->set_double_param("pipe_z", 0);
         	B0Ecal->set_double_param("pipe_hole_r", pipe_hole_r);
          	B0Ecal->set_double_param("cable_x", cable_x);
          	B0Ecal->set_double_param("cable_y", 0);
         	B0Ecal->set_double_param("cable_z", 0);
	        B0Ecal->set_double_param("length", b0Ecal_zLen);
	        B0Ecal->set_double_param("outer_radius", b0_radius);
	        B0Ecal->set_double_param("d_radius", d_radius);
	        B0Ecal->set_string_param("material", "G4_PbWO4");
	        B0Ecal->set_double_param("startAngle",start_angle);
	        B0Ecal->set_double_param("spanningAngle",spanning_angle);
	        B0Ecal->set_double_param("detid",0);
	        B0Ecal->set_double_param("global_x",hFarFwdBeamLine::B0Magnet_x);
	        B0Ecal->set_double_param("global_y",hFarFwdBeamLine::B0Magnet_y);
	        B0Ecal->set_double_param("global_z",hFarFwdBeamLine::B0Magnet_z);
		B0Ecal->set_int_param("lightyield",1); 		//Note additional parameter for storing Light Yield in B0 Ecal
		B0Ecal->SetActive(true);
	        if (verbosity)
        	  B0Ecal->Verbosity(verbosity);
	        B0Ecal->OverlapCheck(overlapCheck);
        	B0Ecal->SetMotherSubsystem(hFarFwdBeamLine::B0Magnet);
       		g4Reco->registerSubsystem(B0Ecal);
	}
	else {					//Use this option to have a circular packman-shape of the B0 ECal for plots.
	        auto *B0Ecal = new EICG4B0Subsystem(Form("b0Truth_%d", 2*b0DetNr), 2*b0DetNr);
        	B0Ecal->SuperDetector("b0Truth");
	        B0Ecal->set_double_param("pipe_hole", pipe_hole);
       		B0Ecal->set_double_param("place_x", 0);
	        B0Ecal->set_double_param("place_y", 0);
        	B0Ecal->set_double_param("place_z", b0Ecal_z);
	        B0Ecal->set_double_param("pipe_x", pipe_x);
        	B0Ecal->set_double_param("pipe_y", 0);
	        B0Ecal->set_double_param("pipe_z", 0);
         	B0Ecal->set_double_param("pipe_hole_r", pipe_hole_r);
          	B0Ecal->set_double_param("cable_x", cable_x);
          	B0Ecal->set_double_param("cable_y", 0);
         	B0Ecal->set_double_param("cable_z", 0);
	        B0Ecal->set_double_param("length", b0Ecal_zLen);
	        B0Ecal->set_double_param("outer_radius", b0_radius);
	        B0Ecal->set_double_param("d_radius", d_radius);
	        B0Ecal->set_string_param("material", "G4_PbWO4");
	        B0Ecal->set_double_param("startAngle",start_angle);
	        B0Ecal->set_double_param("spanningAngle",spanning_angle);
	        B0Ecal->set_double_param("detid",2*b0DetNr);
		B0Ecal->SetActive(true);
        	if (verbosity)
		  B0Ecal->Verbosity(verbosity);
	        B0Ecal->OverlapCheck(overlapCheck);
	        B0Ecal->SetMotherSubsystem(hFarFwdBeamLine::B0Magnet);
        	g4Reco->registerSubsystem(B0Ecal);
	}
      
        auto *B0Ecale = new EICG4B0Subsystem(Form("b0Dead_%d", b0DetNr), b0DetNr); //B0 ECal dead layer is the same subsystem as other four dead layers
        B0Ecale->SuperDetector("b0Dead");
      //  B0Ecale->set_int_param("ispipe", 0); //for future pipe implementation
        B0Ecale->set_double_param("pipe_hole", pipe_hole);
        B0Ecale->set_double_param("place_x", 0);
        B0Ecale->set_double_param("place_y", 0);
        B0Ecale->set_double_param("place_z", b0Ecal_z + (b0Ecal_zLen + b0Cu_zLen)/2);
        B0Ecale->set_double_param("pipe_x", pipe_x);
        B0Ecale->set_double_param("pipe_y", 0);
        B0Ecale->set_double_param("pipe_z", 0);
          B0Ecale->set_double_param("pipe_hole_r", pipe_hole_r);
          B0Ecale->set_double_param("cable_x", cable_x);
          B0Ecale->set_double_param("cable_y", 0);
          B0Ecale->set_double_param("cable_z", 0);
        B0Ecale->set_double_param("length", b0Cu_zLen);
        B0Ecale->set_double_param("d_radius", d_radius);
        B0Ecale->set_double_param("outer_radius", b0_radius);
        B0Ecale->set_string_param("material", "G4_Cu");
        B0Ecale->set_double_param("startAngle",start_angle);
        B0Ecale->set_double_param("spanningAngle",spanning_angle);
        B0Ecale->set_double_param("detid",b0DetNr+1);
        //B0Ecale->SetActive(true);
        B0Ecale->SetActive(false);
        if (verbosity)
          B0Ecale->Verbosity(verbosity);
        B0Ecale->OverlapCheck(overlapCheck);
        B0Ecale->SetMotherSubsystem(hFarFwdBeamLine::B0Magnet);
        g4Reco->registerSubsystem(B0Ecale);
	}
    } else {

       if (Enable::B0_FULLHITPLANE) {

	// Choice 2 circular hit planes
	cout << "Circular hit planes"<<endl;

       	    for (int i = 0; i < b0DetNr; i++)
       	    {
	      b0tr_z = b0tr[i] - b0Mag_zLen / 2;
       	      auto *detB0 = new PHG4CylinderSubsystem(Form("b0Truth_%d", i), i);
              detB0->SuperDetector("b0Truth");
              //detB0->SuperDetector(Form("b0Truth_%d", i));
       	      detB0->set_double_param("radius", 0);
       	      detB0->set_double_param("thickness", 20);
       	      detB0->set_double_param("length", 0.1);
       	      detB0->set_string_param("material", "G4_Si");
       	      detB0->set_double_param("place_z", b0tr_z);  // relative to B0 magnet
       	      detB0->SetActive(true);
       	      if (verbosity) detB0->Verbosity(verbosity);
       	      detB0->OverlapCheck(overlapCheck);
       	    
       	      detB0->SetMotherSubsystem(hFarFwdBeamLine::B0Magnet);
       	    
       	      g4Reco->registerSubsystem(detB0);
          if (Enable::B0TRACKING){ 
		  if (B0TRACKING::FastKalmanFilter)
   		  {
   	 		B0TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + Form("b0Truth_%d", i) ,           //      const std::string& phg4hitsNames,
        	                             B0TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             G4B0TRACKING::PositionResolution,           //      const float radres,
                                             G4B0TRACKING::PositionResolution,           //      const float phires,
                                             0,              //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
	   		B0TRACKING::FastKalmanFilter->add_zplane_state(Form("b0Truth_%d", i), b0Mag_zCent+b0tr_z);
   		 	B0TRACKING::FastKalmanFilterB0Track->add_phg4hits(string("G4HIT_") + Form("b0Truth_%d", i) ,           //      const std::string& phg4hitsNames,
        	                             B0TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             G4B0TRACKING::PositionResolution,           //      const float radres,
                                             G4B0TRACKING::PositionResolution,           //      const float phires,
                                             0,              //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
   			 B0TRACKING::FastKalmanFilterB0Track->add_zplane_state(Form("b0Truth_%d", i), b0Mag_zCent+b0tr_z);
	   		 B0TRACKING::B0ProjectionNames.insert(Form("b0Truth_%d", i));
 		 }
	  }
       	    
       	    }

	} else {

        /// Fun4All default B0 planes
	/// Choice 3 Hit planes with real detector geometry
	cout << "Realistic hit planes"<<endl;
	    
	    for (int i = 0; i < b0DetNr; i++) {
	if (Enable::B0_VAR_PIPE_HOLE){
		pipe_hole = b0tr[i]*cross_angle;
		pipe_x =  cross_angle*(b0Mag_zCent - b0Mag_zLen/2 + b0tr[i]/2) - hFarFwdBeamLine::B0Magnet_x;
	}
	else if (Enable::B0_CIRCLE_PIPE_HOLE){
		pipe_hole = 0.1;
		pipe_hole_r = pipe_hole_r + b0tr[b0DetNr-1]*cross_angle/2;
		pipe_x = - cross_angle*(b0Mag_zCent - b0Mag_zLen/2 + b0tr[b0DetNr-1]/2) - hFarFwdBeamLine::B0Magnet_x;
	}
	else {
		pipe_hole = b0tr[b0DetNr-1]*cross_angle;
		pipe_x =  cross_angle*(b0Mag_zCent - b0Mag_zLen/2 + b0tr[b0DetNr-1]/2) - hFarFwdBeamLine::B0Magnet_x;
	}
	cout <<"Starting B0 Tracker layer "<<i+1<<endl;
	cout <<"Pipe Hole: "<< pipe_hole<<"\t"<<pipe_x<<endl;
	  b0tr_z = b0tr[i] - b0Mag_zLen / 2;
          auto *detB0 = new EICG4B0Subsystem(Form("b0Truth_%d", i), i);
          detB0->SuperDetector(Form("b0Truth_%d", i));
          detB0->set_double_param("place_x", 0);
          detB0->set_double_param("place_y", 0);
      //  detB0->set_int_param("ispipe", 0); //for future pipe implementation
          detB0->set_double_param("pipe_hole", pipe_hole);
          detB0->set_double_param("cable_hole", cable_hole);
          detB0->set_double_param("outer_radius", b0_radius);
          detB0->set_double_param("d_radius", d_radius);
          detB0->set_double_param("length", b0Si_zLen);
          detB0->set_string_param("material", "G4_Si");
          detB0->set_double_param("startAngle",start_angle);
          detB0->set_double_param("spanningAngle",spanning_angle);
          detB0->set_double_param("detid",i);
          detB0->set_double_param("pipe_x", pipe_x);
          detB0->set_double_param("pipe_y", 0);
          detB0->set_double_param("pipe_z", 0);
          detB0->set_double_param("pipe_hole_r", pipe_hole_r);
          detB0->set_double_param("cable_x", cable_x);
          detB0->set_double_param("cable_y", 0);
          detB0->set_double_param("cable_z", 0);
          detB0->set_double_param("place_z", b0tr_z);  // relative to B0 magnet
	        detB0->SetActive(true);
	        if (verbosity)
	          detB0->Verbosity(verbosity);
	        detB0->OverlapCheck(overlapCheck);
	        detB0->SetMotherSubsystem(hFarFwdBeamLine::B0Magnet);
	        g4Reco->registerSubsystem(detB0);
          if (Enable::B0TRACKING){ 
		  if (B0TRACKING::FastKalmanFilter)
   		  {
   	 		B0TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + Form("b0Truth_%d", i) ,           //      const std::string& phg4hitsNames,
        	                             B0TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             G4B0TRACKING::PositionResolution,           //      const float radres,
                                             G4B0TRACKING::PositionResolution,           //      const float phires,
                                             0,              //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
	   		B0TRACKING::FastKalmanFilter->add_zplane_state(Form("b0Truth_%d", i), b0Mag_zCent+b0tr_z);
   		 	B0TRACKING::FastKalmanFilterB0Track->add_phg4hits(string("G4HIT_") + Form("b0Truth_%d", i) ,           //      const std::string& phg4hitsNames,
        	                             B0TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             G4B0TRACKING::PositionResolution,           //      const float radres,
                                             G4B0TRACKING::PositionResolution,           //      const float phires,
                                             0,              //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
   			 B0TRACKING::FastKalmanFilterB0Track->add_zplane_state(Form("b0Truth_%d", i), b0Mag_zCent+b0tr_z);
	   		 B0TRACKING::B0ProjectionNames.insert(Form("b0Truth_%d", i));
 		 }
	  }
	    }

	} 
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
