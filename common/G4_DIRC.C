#ifndef MACRO_G4DIRC_C
#define MACRO_G4DIRC_C

#include <GlobalVariables.C>

#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4detectors/PHG4SectorSubsystem.h>
#include <g4detectors/PHG4ConeSubsystem.h>

#include <g4main/PHG4Reco.h>

#include <cmath>

R__LOAD_LIBRARY(libg4detectors.so)

/*!
 * \file G4_DIRC.C
 * \brief Macro setting up the barrel DIRC
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.3 $
 * \date $Date: 2013/10/09 01:08:17 $
 */

namespace Enable
{
  bool DIRC = false;
  bool DIRC_OVERLAPCHECK = false;
}  // namespace Enable

namespace G4DIRC
{
  
  
  double z_prism  = 30;
  double z_start  = -275+z_prism;
  double z_end    = 125;
  double length   = 0.;
  double z_shift  = 0.;
  double outer_skin_radius = 89.25;
  double dRad     = 5.6;
  double dInSkin  = 7.54;
  namespace SETTING
  {
    bool USECEMCGeo = true;
  }
  
}  // namespace G4DIRC

void DIRCInit()
{
  if (!G4DIRC::SETTING::USECEMCGeo){
    G4DIRC::outer_skin_radius = 78.;
    G4DIRC::dRad              = 8.;
    G4DIRC::dInSkin           = 9.;
    G4DIRC::z_end             = -287 + G4DIRC::z_prism;
    G4DIRC::z_end             = +168;
  }
  G4DIRC::z_shift           = 0.5 * (G4DIRC::z_end + G4DIRC::z_start);
  G4DIRC::length            = G4DIRC::z_end - G4DIRC::z_start;
  
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4DIRC::outer_skin_radius);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, G4DIRC::z_end);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, G4DIRC::z_start- G4DIRC::z_prism);
}

//! Babar DIRC (Without most of support structure)
//! Ref: I. Adam et al. The DIRC particle identification system for the BaBar experiment.
//! Nucl. Instrum. Meth., A538:281-357, 2005. doi:10.1016/j.nima.2004.08.129.
double DIRCSetup(PHG4Reco *g4Reco)
{
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::DIRC_OVERLAPCHECK;

  double radiator_R = G4DIRC::outer_skin_radius-G4DIRC::dRad;
  
  PHG4SectorSubsystem *dirc;
  dirc = new PHG4SectorSubsystem("DIRC");
  dirc->get_geometry().set_normal_polar_angle(M_PI / 2);
  dirc->get_geometry().set_normal_start(radiator_R * PHG4Sector::Sector_Geometry::Unit_cm());
  dirc->get_geometry().set_min_polar_angle(atan2(radiator_R, G4DIRC::z_end));
  dirc->get_geometry().set_max_polar_angle(atan2(radiator_R, G4DIRC::z_start));
  dirc->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
  dirc->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
  dirc->get_geometry().set_material("Quartz");
  dirc->get_geometry().set_N_Sector(12);
  dirc->OverlapCheck(OverlapCheck);
  dirc->get_geometry().AddLayer("Radiator", "Quartz", 1.7 * PHG4Sector::Sector_Geometry::Unit_cm(), true);
  g4Reco->registerSubsystem(dirc);

  PHG4CylinderSubsystem *cyl;

  //  The cylinder skins provide most of the strength
  //  and stiffness of the CST. The thickness of the inner
  //  and outer skins is 1.27 and 0.76 mm, respectively

  double inner_R    = G4DIRC::outer_skin_radius-G4DIRC::dInSkin;
  // Inner skin:
  cyl = new PHG4CylinderSubsystem("DIRC_CST_Inner_Skin", 10);
  cyl->set_double_param("radius", inner_R);
  cyl->set_double_param("length", G4DIRC::length + G4DIRC::z_prism);
  cyl->set_string_param("material", "G4_Al");
  cyl->set_double_param("thickness", 0.127);
  cyl->set_double_param("place_x", 0.);
  cyl->set_double_param("place_y", 0.);
  cyl->set_double_param("place_z", G4DIRC::z_shift - G4DIRC::z_prism * 0.5);
  cyl->SetActive(0);
  cyl->SuperDetector("DIRC");
  cyl->OverlapCheck(OverlapCheck);

  g4Reco->registerSubsystem(cyl);

  // Outer skin:
  cyl = new PHG4CylinderSubsystem("DIRC_CST_Outer_Skin", 11);
  cyl->set_double_param("radius", G4DIRC::outer_skin_radius - 0.076);
  cyl->set_double_param("length", G4DIRC::length);
  cyl->set_string_param("material", "G4_Al");
  cyl->set_double_param("thickness", 0.076);
  cyl->set_double_param("place_x", 0.);
  cyl->set_double_param("place_y", 0.);
  cyl->set_double_param("place_z", G4DIRC::z_shift);
  cyl->SetActive(0);
  cyl->SuperDetector("DIRC");
  cyl->OverlapCheck(OverlapCheck);

  g4Reco->registerSubsystem(cyl);

  // simple approximation for DIRC prism
  PHG4ConeSubsystem *cone = new PHG4ConeSubsystem("DIRC_Prism");
  cone->set_color(0, 1, 0);
  cone->SetR1(radiator_R, radiator_R + 20);
  cone->SetR2(radiator_R, radiator_R + 2);
  cone->SetZlength(0.5 * G4DIRC::z_prism);
  cone->SetPlaceZ(G4DIRC::z_start - 0.5 * G4DIRC::z_prism);
  cone->SetMaterial("Quartz");
  cone->SetActive(0);
  cone->SuperDetector("DIRC");
  cone->OverlapCheck(OverlapCheck);
  g4Reco->registerSubsystem(cone);
  
  // Done
  return G4DIRC::outer_skin_radius;
}
#endif
