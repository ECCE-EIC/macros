#ifndef MACRO_GLOBALVARIABLES_C
#define MACRO_GLOBALVARIABLES_C

#include <g4decayer/EDecayType.hh>
#include <set>
#include <string>

double no_overlapp = 0.0001;

// These Input settings are needed in multiple Input selections
// Putting those here avoids include file ordering problems
namespace Input
{
  bool HEPMC = false;
  bool EMBED = false;
  bool READEIC = false;

  bool UPSILON = false;
  std::set<int> UPSILON_EmbedIds;
}  // namespace Input

namespace DstOut
{
  string OutputDir = ".";
  string OutputFile = "test.root";
}  // namespace DstOut

// Global settings affecting multiple subsystems
namespace Enable
{
  bool OVERLAPCHECK = false;
  bool ABSORBER = false;
  bool DSTOUT = false;
  bool DSTOUT_COMPRESS = false;
  int VERBOSITY = 0;

  // IP selection require explicit choice in the main macros
  bool IP6 = false;
  bool IP8 = false;

  float HFARFWD_ION_ENERGY = 0;
  float HFARBWD_E_ENERGY = 0;
  TString BEAM_COLLISION_SETTING;

}  // namespace Enable

// every G4 subsystem needs to implement this
// rather than forcing another include file,
// let's put this into the GlobalVariables.C
namespace BlackHoleGeometry
{
  double max_radius = 0.;  // this is needed for the overall dimension of the black hole
  double min_z = 0.;
  double max_z = 0.;
  double gap = no_overlapp;
  bool visible = false;
};  // namespace BlackHoleGeometry

namespace G4P6DECAYER
{
  EDecayType decayType = EDecayType::kAll;
}

// our various tracking macro
class PHG4TrackFastSim;
namespace TRACKING
{
  string TrackNodeName = "TrackMap";

  PHG4TrackFastSim * FastKalmanFilter(nullptr);

  PHG4TrackFastSim * FastKalmanFilterSiliconTrack(nullptr);

  PHG4TrackFastSim * FastKalmanFilterInnerTrack(nullptr);

  std::set<std::string> ProjectionNames;
}

//For B0 Tracking
class B0TrackFastSim;
namespace B0TRACKING
{
  string TrackNodeName = "TrackMap";

  B0TrackFastSim * FastKalmanFilter(nullptr);

  B0TrackFastSim * FastKalmanFilterB0Track(nullptr);

  std::set<std::string> B0ProjectionNames;
}

namespace G4MAGNET
{
  // initialize to garbage values - the override is done in the respective
  // MagnetInit() functions. If used standalone (without the G4_Magnet include)
  // like in the tracking - those need to be set in the Fun4All macro
  double magfield_rescale = NAN;
  string magfield;
}  // namespace G4MAGNET
#endif
