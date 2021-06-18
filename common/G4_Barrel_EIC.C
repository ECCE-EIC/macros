#ifndef MACRO_G4BARREL_C
#define MACRO_G4BARREL_C

#include "GlobalVariables.C"
#include <g4lblvtx/AllSiliconTrackerSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4lblvtx/SimpleNtuple.h>
#include <g4main/PHG4Reco.h>

#include <string>

R__LOAD_LIBRARY(libg4detectors.so)

using namespace std;

namespace Enable
{
  bool BARREL = false;
}  // namespace Enable
//-----------------------------------------------------------------------------------//
void BarrelInit()
{
}
//-----------------------------------------------------------------------------------//
void Barrel(PHG4Reco *g4Reco,int det_ver=3)
{
  // Loading All-Si Tracker from dgml file
  AllSiliconTrackerSubsystem *allsili = new AllSiliconTrackerSubsystem();
  allsili->set_string_param("GDMPath",string(getenv("CALIBRATIONROOT"))+Form("/AllSiliconTracker/genfitGeom_AllSi_v%d.gdml",det_ver));
  allsili->AddAssemblyVolume("VST");      // Barrel
  
  // this is for plotting single logical volumes for debugging and geantino scanning they end up at the center, you can plot multiple
  // allsili->AddLogicalVolume("VstStave00");
  
  allsili->SuperDetector("Barrel");
  allsili->SetActive();          // this saves hits in the MimosaCore volumes
  allsili->SetAbsorberActive();  // this saves hits in all volumes (in the absorber node)
  g4Reco->registerSubsystem(allsili);

  float pitch=10e-4;
  int nBarrel = 6;

  if (TRACKING::FastKalmanFilter)
    {
      for (int i=10; i<10+nBarrel; i++)
	{
	  TRACKING::FastKalmanFilter->add_phg4hits(Form("G4HIT_LBLVTX_CENTRAL_%d", i),//      const std::string& phg4hitsNames,
						   PHG4TrackFastSim::Cylinder,        //      const DETECTOR_TYPE phg4dettype,
						   999,                               //      const float radres, not used
						   pitch / sqrt(12.),                 //      const float phires,
						   pitch / sqrt(12.),                 //      const float lonres, *ignored in plane detector*
						   1,                                 //      const float eff,
						   0);                                //      const float noise
	}
    }

}
//-----------------------------------------------------------------------------------//
void add_AllSi_hits(SimpleNtuple * hits, int det_ver=3)
{
  char nodename[100];
  int nBarrel = 6;
  int nDisks = 5;
  
  if(det_ver==3)
  {
    nBarrel = 7;
    nDisks = 7;
  }
  
  hits->AddNode("ABSORBER_LBLVTX",0); // hits in the passive volumes
  for (int i = 10; i < 10+nBarrel; i++)
  {
    sprintf(nodename, "LBLVTX_CENTRAL_%d", i);    
    hits->AddNode(nodename, i); // hits in the  MimosaCore volumes
  }
  for (int i = 20; i < 20+nDisks ; i++)
  {
    sprintf(nodename, "LBLVTX_FORWARD_%d", i);  
    hits->AddNode(nodename, i); // hits in the  MimosaCore volumes
  }
  for (int i = 30; i < 30+nDisks ; i++)
  {
    sprintf(nodename, "LBLVTX_BACKWARD_%d",i); 
    hits->AddNode(nodename, i); // hits in the  MimosaCore volumes
  }
}

#endif
