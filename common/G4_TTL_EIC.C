#ifndef MACRO_G4TTLEIC_C
#define MACRO_G4TTLEIC_C

#include "GlobalVariables.C"

#include <g4ttl/PHG4TTLSubsystem.h>
#include <g4detectors/PHG4CylinderSubsystem.h>

#include <g4main/PHG4Reco.h>

#include <string>

R__LOAD_LIBRARY(libg4detectors.so)

int make_forward_station(string name, PHG4Reco *g4Reco, double zpos, double Rmin,
                          double Rmax,double tSilicon, double xoffset=0);
int make_forward_station_basic(string name, PHG4Reco *g4Reco, double zpos, double Rmin,
                          double Rmax,double tSilicon);
int make_barrel_layer_basic(string name, PHG4Reco *g4Reco,
                      double radius, double halflength, double tSilicon, double zOffset);
int make_barrel_layer_LYSO_basic(string name, PHG4Reco *g4Reco,
                      double radius, double halflength, double tSilicon, double zOffset);
int make_barrel_layer(string name, PHG4Reco *g4Reco, 
                      double radius, double halflength, double tSilicon, double zOffset);

//-----------------------------------------------------------------------------------//
namespace Enable
{
  bool FTTL = false;
  bool ETTL = false;
  bool CTTL = false;
}

namespace G4TTL
{
  int layer[3]                  = { 2, 1, 2};
  double positionToVtx[3][3]    = { {-169., -172., -309.5}, {80., 114.7, 0. }, { 287., 289., 340.} };
  double minExtension[3][3]     = { {8, 8, 15.3}, {218, 180, 0 }, {11.62, 11.7, 13.8 } };
  double maxExtension[3][3]     = { {61., 61. , 200}, {-40, 0, 0 }, {170., 170., 250  } };
  double xoffsetFTTLIP6[3]      = { -6., -6., -6.};
  double xoffsetFTTLIP8[3]      = { 8.4, 8.4, 8.4};
  namespace SETTING
  {
    bool optionCEMC     = false;
    bool optionEEMCH    = true;
    bool optionBasicGeo = false;
    bool optionLYSO     = false;
    int optionDR        = 0;
    int optionGeo       = 7;
    int optionGran      = 1;
  }  // namespace SETTING



  // 2, LGAD based ToF by Wei Li:
  // Present the recent simulation studies and detector layout for the LGAD based ToF,
  // which can be used as an outer tracker. The current detector can reach around 30 micron spatial
  // resolution with 500 micron AC-LGAD sensor.
  // This detector can provide better coverage and provide further constraints in track
  // reconstruction within high pseudorapidity regions.

  const double PositionResolution(30e-4);

}  // namespace G4TTL


//-----------------------------------------------------------------------------------//
void TTL_Init()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 250.);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 350.);
  
  if (!G4TTL::SETTING::optionEEMCH){
    G4TTL::positionToVtx[0][0]  = -155.5;
    G4TTL::positionToVtx[0][1]  = -158.5;
  }
  
  if (G4TTL::SETTING::optionCEMC){
    G4TTL::maxExtension[0][0]   = 80.;
    G4TTL::maxExtension[0][1]   = 80.;
    G4TTL::positionToVtx[1][0]  = 92.;
    if(!G4TTL::SETTING::optionBasicGeo) G4TTL::positionToVtx[1][0] = 89.;
  } 
  
  if(G4TTL::SETTING::optionBasicGeo){
    G4TTL::minExtension[1][0]   = 180;
    G4TTL::maxExtension[1][0]   = -25;
  }

  if (G4TTL::SETTING::optionGeo == 1){
    cout << "TTL setup infront of ECals with 2 layers fwd/bwd & 1 layer barrel" << endl;  
  }
  if (G4TTL::SETTING::optionGeo == 2){
    cout << "TTL setup infront of ECals with 2 layers fwd/bwd & 1 layer barrel, lower barrel layer" << endl;  
    G4TTL::positionToVtx[1][0] = 50.;
    if(!G4TTL::SETTING::optionBasicGeo) G4TTL::positionToVtx[1][0] = 65.;
    G4TTL::minExtension[1][0]  = 100.;
    G4TTL::maxExtension[1][0]  = 0.;
  }
  if (G4TTL::SETTING::optionGeo == 3 || G4TTL::SETTING::optionGeo == 5 || G4TTL::SETTING::optionGeo == 6){
    cout << "TTL setup infront of ECals with 1 layers fwd/bwd & 1 layer barrel, lower barrel layer" << endl;  
    G4TTL::positionToVtx[1][0] = 50.;
    if(!G4TTL::SETTING::optionBasicGeo) G4TTL::positionToVtx[1][0] = 65.;
    G4TTL::minExtension[1][0]  = 100.;
    G4TTL::maxExtension[1][0]  = 0.;
    G4TTL::layer[0]            = 1;
    G4TTL::layer[2]            = 1;
    G4TTL::positionToVtx[0][0] = G4TTL::positionToVtx[0][1];
    G4TTL::positionToVtx[2][0] = G4TTL::positionToVtx[2][1];
    G4TTL::minExtension[0][0]  = G4TTL::minExtension[0][1];
    G4TTL::minExtension[2][0]  = G4TTL::minExtension[2][1];
    G4TTL::maxExtension[0][0]  = G4TTL::maxExtension[0][1];
    G4TTL::maxExtension[2][0]  = G4TTL::maxExtension[2][1];
  }
  if (G4TTL::SETTING::optionGeo == 4){
    cout << "TTL setup infront of ECals  with 2 layers fwd/bwd & 1 layer barrel, 1 layer before HCals everywhere" << endl;  
    G4TTL::layer[0]    = 3;
    G4TTL::layer[1]    = 2;
    if(!G4TTL::SETTING::optionBasicGeo) G4TTL::layer[1]    = 1;
    G4TTL::layer[2]    = 3;
  }
  if (G4TTL::SETTING::optionGeo == 5){
    // Option 5 is 1 layer fwd/bwd, with LYSO in the central barrel. We use the geometry for option 3 since it's the same.
    // However, we create a separate option to simplify setting up the configuration.
    cout << "TTL setup using LYSO in central barrel" << endl;
    G4TTL::SETTING::optionLYSO = true;
  }
  if(G4TTL::SETTING::optionGeo == 7){
    cout << "TTL one forward disk in front of dRICH and one backward disk in front of EEMC, barrel CTTL center at radius 64cm" << endl;
    // single disk in front of dRICH (full eta)
    G4TTL::layer[2]             = 1;
    G4TTL::minExtension[2][0]   = 7.0;
    G4TTL::maxExtension[2][0]   = 87;
    G4TTL::positionToVtx[2][0]  = 182.;
    G4TTL::xoffsetFTTLIP6[0]    = -2.7;
    G4TTL::xoffsetFTTLIP8[0]    = 3.0;

    // single disk in front of EEMC
    G4TTL::layer[0]             = 1;
    G4TTL::maxExtension[0][0]   = 64;

    // barrel layer at 64cm
    G4TTL::positionToVtx[1][0]  = 64.;
    G4TTL::minExtension[1][0]   = 140;
    G4TTL::maxExtension[1][0]   = 0;
  }
  if(G4TTL::SETTING::optionGeo == 8){
    G4TTL::layer[2]             = 3;
    cout << "TTL forward disk 1 in front of dRICH" << endl;
    G4TTL::minExtension[2][0]   = 7.0;
    G4TTL::maxExtension[2][0]   = 87;
    G4TTL::positionToVtx[2][0]  = 182.;
    cout << "additional two small TTL disks in front of FEMC" << endl;
    G4TTL::minExtension[2][1]   = 11.62;
    G4TTL::minExtension[2][2]   = 11.7;
    G4TTL::maxExtension[2][1]   = 60.;
    G4TTL::maxExtension[2][2]   = 60.;
    G4TTL::positionToVtx[2][1]  = 287.;
    G4TTL::positionToVtx[2][2]  = 289.;

    G4TTL::layer[0]             = 1;
    G4TTL::maxExtension[0][0]   = 64;

    G4TTL::positionToVtx[1][0]  = 64.;
    G4TTL::minExtension[1][0]   = 140;
    G4TTL::maxExtension[1][0]   = 0;
  }

  if (G4TTL::SETTING::optionDR == 2 && G4TTL::SETTING::optionGeo == 4 ){
     cout << "TTL setup infront of ECals  with 2 layers fwd/bwd & 1 layer barrel, 1 layer before HCals everywhere choosen!" << endl;  
     cout << "conflicting in forward region with DR calo, reducing by 1 layer!" << endl;  
     G4TTL::layer[2]    = 2;
  } else if (G4TTL::SETTING::optionDR == 1 && G4TTL::SETTING::optionGeo == 4 ){
     cout << "TTL setup infront of ECals  with 2 layers fwd/bwd & 1 layer barrel, 1 layer before HCals everywhere choosen!" << endl;  
     cout << "conflicting in forward region with DR calo, reducing adding larger cutout!" << endl;  
     G4TTL::minExtension[2][2]   = 2.5;
  }   
  
}
//-----------------------------------------------------------------------------------//
void FTTLSetup(PHG4Reco *g4Reco, TString fttloption = "")
{
  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = .1 * cm;
  const double um = 1e-3 * mm;

  for (Int_t i = 0; i < G4TTL::layer[2]; i++){
    cout << G4TTL::positionToVtx[2][i] << "\t" << G4TTL::minExtension[2][i] << "\t" << G4TTL::maxExtension[2][i] << endl;
    if(!G4TTL::SETTING::optionBasicGeo){
      make_forward_station(Form("FTTL_%d", i), g4Reco, G4TTL::positionToVtx[2][i],  G4TTL::minExtension[2][i], G4TTL::maxExtension[2][i], 85*um, Enable::IP8 ? G4TTL::xoffsetFTTLIP8[i] : G4TTL::xoffsetFTTLIP6[i]);
    } else {
      make_forward_station_basic(Form("FTTL_%d", i), g4Reco, G4TTL::positionToVtx[2][i],  G4TTL::minExtension[2][i], G4TTL::maxExtension[2][i], 85*um);
    }
  }
}


//-----------------------------------------------------------------------------------//
void ETTLSetup(PHG4Reco *g4Reco, TString ettloption = "")
{
  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = .1 * cm;
  const double um = 1e-3 * mm;
  for (Int_t i = 0; i < G4TTL::layer[0]; i++){
    cout << G4TTL::positionToVtx[0][i] << "\t" << G4TTL::minExtension[0][i] << "\t" << G4TTL::maxExtension[0][i] << endl;
    if(!G4TTL::SETTING::optionBasicGeo){
      make_forward_station(Form("ETTL_%d", i), g4Reco, G4TTL::positionToVtx[0][i],  G4TTL::minExtension[0][i], G4TTL::maxExtension[0][i], 85*um, 0);
    } else {
      make_forward_station_basic(Form("ETTL_%d", i), g4Reco, G4TTL::positionToVtx[0][i],  G4TTL::minExtension[0][i], G4TTL::maxExtension[0][i], 85*um);
    }
  }
}

//-----------------------------------------------------------------------------------//
void CTTLSetup(PHG4Reco *g4Reco, TString cttloption = "")
{
  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = .1 * cm;
  const double um = 1e-3 * mm;

  for (Int_t i = 0; i < G4TTL::layer[1]; i++){
    cout << "Radius: " << G4TTL::positionToVtx[1][i] << "\tLength: " << G4TTL::minExtension[1][i] << "\tz-Offset: " << G4TTL::maxExtension[1][i] << endl;
     if(G4TTL::SETTING::optionBasicGeo && G4TTL::SETTING::optionLYSO){
      make_barrel_layer_LYSO_basic(Form("CTTL_%d",i), g4Reco, G4TTL::positionToVtx[1][i],  G4TTL::minExtension[1][i], 85*um, G4TTL::maxExtension[1][i]);
    } else if(G4TTL::SETTING::optionBasicGeo){
      make_barrel_layer_basic(Form("CTTL_%d",i), g4Reco, G4TTL::positionToVtx[1][i],  G4TTL::minExtension[1][i], 85*um, G4TTL::maxExtension[1][i]);
    } else {
      make_barrel_layer(Form("CTTL_%d",i), g4Reco, G4TTL::positionToVtx[1][i],  G4TTL::minExtension[1][i], 85*um, G4TTL::maxExtension[1][i]);
    }
  }
}


//-----------------------------------------------------------------------------------//
int make_forward_station(string name, PHG4Reco *g4Reco,
        double zpos, double rMin, double rMax,
        double tSilicon, //silicon thickness
        double xoffset = 0 )
{
  cout << "r min: " << rMin << "\t r max: " << rMax << "\t z: " <<  zpos << endl;
  
  // always facing the interaction point
  double polar_angle = 0;
  double place_z(zpos);
  if (place_z < 0){
    place_z = -place_z;
    polar_angle = M_PI;
  }
  PHG4TTLSubsystem *ttl;
  ttl = new PHG4TTLSubsystem(name);
  ttl->SetDetailed(false);
  ttl->SuperDetector(name);
  ttl->set_double_param("polar_angle", polar_angle);                    //
  ttl->set_int_param("isForward", 1);                    //
  ttl->set_double_param("place_z", place_z * cm);                    //
  ttl->set_double_param("rMin", rMin * cm);                    //
  ttl->set_double_param("rMax", rMax * cm);                    //
  ttl->set_double_param("offset_x", xoffset * cm);                    //
  ttl->set_double_param("tSilicon", tSilicon);                    //
  // ttl->OverlapCheck(true);
  ttl->OverlapCheck(Enable::OVERLAPCHECK);
  
  g4Reco->registerSubsystem(ttl);


  if (TRACKING::FastKalmanFilter)
  {
    TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             G4TTL::PositionResolution,           //      const float radres,
                                             G4TTL::PositionResolution,           //      const float phires,
                                             tSilicon / sqrt(12.),              //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
    TRACKING::FastKalmanFilter->add_zplane_state(name, zpos);
    TRACKING::ProjectionNames.insert(name);
  }
  if (TRACKING::FastKalmanFilterInnerTrack and zpos>0)
  {
    TRACKING::FastKalmanFilterInnerTrack->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             G4TTL::PositionResolution,           //      const float radres,
                                             G4TTL::PositionResolution,           //      const float phires,
                                             tSilicon / sqrt(12.),              //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
  }

  return 0;
}


//-----------------------------------------------------------------------------------//
int make_forward_station_basic(string name, PHG4Reco *g4Reco,
        double zpos, double rmin, double rmax,
        double tSilicon) //silicon thickness
{
  //  cout
  //      << "make_GEM_station - GEM construction with PHG4SectorSubsystem - make_GEM_station_EdgeReadout  of "
  //      << name << endl;

  // always facing the interaction point
  double etamin = -TMath::Log(TMath::Tan(rmin/TMath::Abs(zpos)/2));
  double etamax = -TMath::Log(TMath::Tan(rmax/TMath::Abs(zpos)/2));

  double polar_angle = 0;
  if (zpos < 0){
    zpos = -zpos;
    polar_angle = M_PI;
    etamin = -etamin;
    etamax = -etamax;
  }
  if (etamax < etamin){
    double t = etamax;
    etamax = etamin;
    etamin = t;
  }

  PHG4SectorSubsystem *ttl;
  ttl = new PHG4SectorSubsystem(name);

  ttl->SuperDetector(name);

  ttl->get_geometry().set_normal_polar_angle(polar_angle);
  ttl->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm());
  ttl->get_geometry().set_min_polar_angle(PHG4Sector::Sector_Geometry::eta_to_polar_angle(etamax));
  ttl->get_geometry().set_max_polar_angle(PHG4Sector::Sector_Geometry::eta_to_polar_angle(etamin));
  ttl->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::ConeEdge());
  ttl->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::ConeEdge());
  ttl->get_geometry().set_N_Sector(1);
  ttl->get_geometry().set_material("G4_AIR");
  ttl->OverlapCheck(Enable::OVERLAPCHECK);
  
  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = .1 * cm;
  const double um = 1e-3 * mm;
  // build up layers

  ttl->get_geometry().AddLayer("SiliconSensor", "G4_Si", tSilicon, true, 100);
  ttl->get_geometry().AddLayer("Metalconnection", "G4_Al", 100 * um, false, 100);
  ttl->get_geometry().AddLayer("HDI", "G4_KAPTON", 20 * um, false, 100);
  ttl->get_geometry().AddLayer("Cooling", "G4_WATER", 100 * um, false, 100);
  ttl->get_geometry().AddLayer("Support", "G4_GRAPHITE", 50 * um, false, 100);
  ttl->get_geometry().AddLayer("Support_Gap", "G4_AIR", 1 * cm, false, 100);
  ttl->get_geometry().AddLayer("Support2", "G4_GRAPHITE", 50 * um, false, 100);

  g4Reco->registerSubsystem(ttl);


  if (TRACKING::FastKalmanFilter)
  {
    TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             G4TTL::PositionResolution,           //      const float radres,
                                             G4TTL::PositionResolution,           //      const float phires,
                                             tSilicon / sqrt(12.),              //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
    TRACKING::FastKalmanFilter->add_zplane_state(name, zpos);
    TRACKING::ProjectionNames.insert(name);
  }
  if (TRACKING::FastKalmanFilterInnerTrack and zpos>0)
  {
    TRACKING::FastKalmanFilterInnerTrack->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             G4TTL::PositionResolution,           //      const float radres,
                                             G4TTL::PositionResolution,           //      const float phires,
                                             tSilicon / sqrt(12.),              //      const float lonres, *ignored in plane detector*
                                             1,                                 //      const float eff,
                                             0);                                //      const float noise
  }

  return 0;
}

//-----------------------------------------------------------------------------------//
int make_barrel_layer(string name, PHG4Reco *g4Reco, 
                      double radius, double halflength, double tSilicon, double zOffset )
{
  // cout << "r min: " << rMin << "\t r max: " << rMax << "\t z: " <<  zpos << endl;
  if(G4TTL::SETTING::optionCEMC){
    cout << "The improved barrel TTL layer is placed at a fixed radius of rMin = 80 * cm and is meant to only be used with the BECAL!" << endl;
  }
  PHG4TTLSubsystem *ttl;
  ttl = new PHG4TTLSubsystem(name);
  ttl->SetDetailed(false);
  ttl->SuperDetector(name);
  ttl->set_int_param("isForward", 0);                    //
  ttl->set_double_param("place_z", zOffset * cm);                    //
  ttl->set_double_param("rMin", radius * cm);                    //
  ttl->set_double_param("length", 2.0 * halflength * cm);
  ttl->set_double_param("tSilicon", tSilicon);                    //
  ttl->OverlapCheck(Enable::OVERLAPCHECK);

  g4Reco->registerSubsystem(ttl);


  if (TRACKING::FastKalmanFilter)
  {
    TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + name,     //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
                                             tSilicon / sqrt(12.),        //      const float radres,
                                             G4TTL::PositionResolution,     //      const float phires,
                                             G4TTL::PositionResolution,     //      const float lonres,
                                             1,                           //      const float eff,
                                             0);                          //      const float noise
    TRACKING::FastKalmanFilter->add_cylinder_state(name, radius);

    TRACKING::ProjectionNames.insert(name);
  }
  if (TRACKING::FastKalmanFilterInnerTrack)
  {
    TRACKING::FastKalmanFilterInnerTrack->add_phg4hits(string("G4HIT_") + name,     //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
                                             tSilicon / sqrt(12.),        //      const float radres,
                                             G4TTL::PositionResolution,     //      const float phires,
                                             G4TTL::PositionResolution,     //      const float lonres,
                                             1,                           //      const float eff,
                                             0);                          //      const float noise

  }

  return 0;
}




//-----------------------------------------------------------------------------------//
int make_barrel_layer_basic(string name, PHG4Reco *g4Reco, 
                      double radius, double halflength, double tSilicon, double zOffset){

  //---------------------------------
  //build barrel layer
  //---------------------------------
  const int nSubLayer = 7;

  string layerName[nSubLayer] = {"SiliconSensor", "Metalconnection", "HDI", "Cooling",
                                 "Support1", "Support_Gap", "Support2"};
  string material[nSubLayer] = {"G4_Si", "G4_Al", "G4_KAPTON", "G4_WATER",
                                "G4_GRAPHITE", "G4_AIR", "G4_GRAPHITE"};
  double thickness[nSubLayer] = {tSilicon , 15 * um, 20 * um, 100 * um,
                                 50 * um, 1, 50 * um};

  double max_bh_radius = 0.;
  PHG4CylinderSubsystem* cyl;
//   cout << "started to create cylinder layer: " << name << endl;
  
  double currRadius = radius;
//   cout << currRadius << endl;
  for (int l = 0; l < nSubLayer; l++) {
//     cout << name <<"_"<< layerName[l] << endl;
    cyl = new PHG4CylinderSubsystem(name + "_" + layerName[l],l);
    cyl->SuperDetector(name);
    cyl->set_double_param("radius", currRadius);
    cyl->set_double_param("length", 2.0 * halflength);
    cyl->set_string_param("material", material[l]);
    cyl->set_double_param("thickness", thickness[l]);
    cyl->set_double_param("place_x", 0.);
    cyl->set_double_param("place_y", 0.);
    cyl->set_double_param("place_z", zOffset);

    if (l == 0) cyl->SetActive();  //only the Silicon Sensor is active
    cyl->OverlapCheck(true);
    g4Reco->registerSubsystem(cyl);
    currRadius = currRadius+thickness[l];
//     cout << currRadius << endl;
  }


  if (TRACKING::FastKalmanFilter)
  {
    TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + name,     //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
                                             tSilicon / sqrt(12.),        //      const float radres,
                                             G4TTL::PositionResolution,     //      const float phires,
                                             G4TTL::PositionResolution,     //      const float lonres,
                                             1,                           //      const float eff,
                                             0);                          //      const float noise
    TRACKING::FastKalmanFilter->add_cylinder_state(name, radius);

    TRACKING::ProjectionNames.insert(name);
  }
  if (TRACKING::FastKalmanFilterInnerTrack)
  {
    TRACKING::FastKalmanFilterInnerTrack->add_phg4hits(string("G4HIT_") + name,     //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
                                             tSilicon / sqrt(12.),        //      const float radres,
                                             G4TTL::PositionResolution,     //      const float phires,
                                             G4TTL::PositionResolution,     //      const float lonres,
                                             1,                           //      const float eff,
                                             0);                          //      const float noise

  }

  return 0;
}

//-----------------------------------------------------------------------------------//
int make_barrel_layer_LYSO_basic(string name, PHG4Reco *g4Reco,
                      double radius, double halflength, double tSilicon, double zOffset){

  //---------------------------------
  //build barrel layer
  //---------------------------------
  const int nSubLayer = 4;

  string layerName[nSubLayer] = {"Cooling", "Crystal","SiliconSensor", "Motherboard"};
  string material[nSubLayer] = {"G4_Al", "LSO", "G4_Si", "FR4"};
  double thickness[nSubLayer] = {0.031 * cm, 0.029 * cm, tSilicon, 0.033 * cm};

  double max_bh_radius = 0.;
  PHG4CylinderSubsystem* cyl;
//   cout << "started to create cylinder layer: " << name << endl;

  double currRadius = radius;
//   cout << currRadius << endl;
  for (int l = 0; l < nSubLayer; l++) {
//     cout << name <<"_"<< layerName[l] << endl;
    cyl = new PHG4CylinderSubsystem(name + "_" + layerName[l],l);
    cyl->SuperDetector(name);
    cyl->set_double_param("radius", currRadius);
    cyl->set_double_param("length", 2.0 * halflength);
    cyl->set_string_param("material", material[l]);
    cyl->set_double_param("thickness", thickness[l]);
    cyl->set_double_param("place_x", 0.);
    cyl->set_double_param("place_y", 0.);
    cyl->set_double_param("place_z", zOffset);

    if (l == 2) cyl->SetActive();  //only the Silicon Sensor is active
    cyl->OverlapCheck(Enable::OVERLAPCHECK);
    g4Reco->registerSubsystem(cyl);
    currRadius = currRadius+thickness[l];
//     cout << currRadius << endl;
  }

  if (TRACKING::FastKalmanFilter)
  {
    float resLGAD_barrel = G4TTL::PositionResolution;
    if(G4TTL::SETTING::optionLYSO){
      resLGAD_barrel = 35e-1; // https://cds.cern.ch/record/2667167/files/CMS-TDR-020.pdf page 33 bottom
    }
    TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + name,     //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
                                             999,                         //      const float radres,
                                             resLGAD_barrel,              //      const float phires,
                                             resLGAD_barrel,              //      const float lonres,
                                             1,                           //      const float eff, NOTE: Different from 0.95 in Modular
                                             0);                          //      const float noise
    TRACKING::FastKalmanFilter->add_cylinder_state(name, radius);

    TRACKING::ProjectionNames.insert(name);
  }

  return 0;
}

#endif

//-----------------------------------------------------------------------------------//
