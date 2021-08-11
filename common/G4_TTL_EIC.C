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
  double positionToVtx[3][3]    = { {-185.5, -188.5, -309.5}, {80., 114.7, 0. }, { 287., 289., 340.} };
  double minExtension[3][3]     = { {10, 10, 15.3}, {180, 180, 0 }, {11.62, 11.7, 13.8 } };
  double maxExtension[3][3]     = { {67., 67. , 200}, {-25, 0, 0 }, {170., 170., 250  } };
  namespace SETTING
  {
    bool optionCEMC  = true;
    bool optionEEMCH = true;
    int optionDR    = 0;
    int optionGeo   = 1;
    int optionGran  = 1;
  }  // namespace SETTING
}  // namespace G4FHCAL


//-----------------------------------------------------------------------------------//
void TTL_Init()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 200.);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 350.);
  
  if (!G4TTL::SETTING::optionEEMCH){
    G4TTL::positionToVtx[0][0]  = -155.5;
    G4TTL::positionToVtx[0][1]  = -158.5;
  }
  
  if (G4TTL::SETTING::optionCEMC){
    G4TTL::maxExtension[0][0]   = 80.;
    G4TTL::maxExtension[0][1]   = 80.;    
    G4TTL::positionToVtx[1][0]  = 92.;
  } 
    
  if (G4TTL::SETTING::optionGeo == 1){
    cout << "TTL setup infront of ECals with 2 layers fwd/bwd & 1 layer barrel" << endl;  
  } if (G4TTL::SETTING::optionGeo == 2){
    cout << "TTL setup infront of ECals with 2 layers fwd/bwd & 1 layer barrel, lower barrel layer" << endl;  
    G4TTL::positionToVtx[1][0] = 50.;
    G4TTL::minExtension[1][0]  = 100.;
    G4TTL::maxExtension[1][0]  = 0.;
  } if (G4TTL::SETTING::optionGeo == 3){
    cout << "TTL setup infront of ECals with 1 layers fwd/bwd & 1 layer barrel, lower barrel layer" << endl;  
    G4TTL::positionToVtx[1][0] = 50.;
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
  } if (G4TTL::SETTING::optionGeo == 4){
    cout << "TTL setup infront of ECals  with 2 layers fwd/bwd & 1 layer barrel, 1 layer before HCals everywhere" << endl;  
    G4TTL::layer[0]    = 3;
    G4TTL::layer[1]    = 2;
    G4TTL::layer[2]    = 3;
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
    make_forward_station(Form("FTTL_%d", i), g4Reco, G4TTL::positionToVtx[2][i],  G4TTL::minExtension[2][i], G4TTL::maxExtension[2][i], 85*um, 6);
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
    make_forward_station(Form("ETTL_%d", i), g4Reco, G4TTL::positionToVtx[0][i],  G4TTL::minExtension[0][i], G4TTL::maxExtension[0][i], 85*um);
  }
}

//-----------------------------------------------------------------------------------//
void CTTLSetup(PHG4Reco *g4Reco, TString cttloption = "")
{
  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = .1 * cm;
  const double um = 1e-3 * mm;
  
  for (Int_t i = 0; i < G4TTL::layer[1]; i++){
    cout << G4TTL::positionToVtx[1][i] << "\t" << G4TTL::minExtension[1][i] << "\t" << G4TTL::maxExtension[1][i] << endl;
    make_barrel_layer(Form("CTTL_%d",i), g4Reco, G4TTL::positionToVtx[1][i],  G4TTL::minExtension[1][i], 85*um, G4TTL::maxExtension[1][i]);     
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
  if (zpos < 0){
    zpos = -zpos;
    polar_angle = M_PI;
  }
  PHG4TTLSubsystem *ttl;
  ttl = new PHG4TTLSubsystem(name);
  ttl->SetDetailed(false);
  ttl->SuperDetector(name);
  ttl->set_double_param("polar_angle", polar_angle);                    //
  ttl->set_double_param("place_z", zpos * cm);                    //
  ttl->set_double_param("rMin", rMin * cm);                    //
  ttl->set_double_param("rMax", rMax * cm);                    //
  ttl->set_double_param("offset_x", xoffset * cm);                    //
  ttl->set_double_param("tSilicon", tSilicon);                    //
  ttl->OverlapCheck(false);

  g4Reco->registerSubsystem(ttl);
  return 0;
}




//-----------------------------------------------------------------------------------//
int make_barrel_layer(string name, PHG4Reco *g4Reco, 
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
    cyl->OverlapCheck(false);
    g4Reco->registerSubsystem(cyl);
    currRadius = currRadius+thickness[l];
//     cout << currRadius << endl;
  }

  return 0;
}

#endif

//-----------------------------------------------------------------------------------//
