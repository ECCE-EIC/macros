#ifndef MACRO_FUN4ALLG4EICDETECTOR_C
#define MACRO_FUN4ALLG4EICDETECTOR_C

#include <GlobalVariables.C>

#include <DisplayOn.C>
#include <G4Setup_EICDetector.C>
#include <G4_DSTReader_EICDetector.C>
#include <G4_EventEvaluator.C>
#include <G4_FwdJets.C>
#include <G4_Global.C>
#include <G4_Input.C>
#include <G4_Production.C>
#include <G4_User.C>

#include <TROOT.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

#include <RooUnblindPrecision.h>

R__LOAD_LIBRARY(libfun4all.so)

int Fun4All_G4_EICDetector(
    const int nEvents = 1,
    const string &inputFile = "https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/sPHENIX_G4Hits_sHijing_9-11fm_00000_00010.root",
    const string &outputFile = "G4EICDetector.root",
    const string &embed_input_file = "https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/sPHENIX_G4Hits_sHijing_9-11fm_00000_00010.root",
    const int skip = 0,
    const string &outdir = ".")
{
  //---------------
  // Fun4All server
  //---------------
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);
  //Opt to print all random seed used for debugging reproducibility. Comment out to reduce stdout prints.
  //PHRandomSeed::Verbosity(1);

  // just if we set some flags somewhere in this macro
  recoConsts *rc = recoConsts::instance();
  // By default every random number generator uses
  // PHRandomSeed() which reads /dev/urandom to get its seed
  // if the RANDOMSEED flag is set its value is taken as initial seed
  // which will produce identical results so you can debug your code
  // rc->set_IntFlag("RANDOMSEED", 12345);

  bool generate_seed = false;

  if (generate_seed)
  {
    size_t findSlash = inputFile.find_last_of("/");
    string inputFileName = inputFile.substr(findSlash + 1, inputFile.size());

    RooRealVar dummyVal("dummy", "", 0);
    RooUnblindPrecision blindVal("blindVal", "blindVal", inputFileName.c_str(), nEvents, skip + 1, dummyVal, kFALSE);
    rc->set_IntFlag("RANDOMSEED", abs(ceil(blindVal.getVal() * 1e2)));
  }

  //===============
  // Input options
  //===============

  // switching IPs by comment/uncommenting the following lines
  // used for both beamline setting and for the event generator crossing boost
  Enable::IP6 = true;
  //Enable::IP8 = true;

   
  //===============
  // The following Ion energy and electron energy setting needs to be speficied
  // The setting options for e-p high divergence setting (p energy x e energy):
  // Option: 275x18, 275x10, 100x10, 100x5, 41x5
  //
  // The setting options for e-p high divergence setting (p energy x e energy):
  // Option: 275x18, 275x10, 100x10, 100x5, 41x5
  //
  // The setting options for e-p high divergence setting (A energy x e energy):
  // Option: 110x18, 110x10, 110x5, 41x5

  // Setting proton beam pipe energy. If you don't know what to set here, leave it at 275
  Enable::HFARFWD_ION_ENERGY = 275;

  // Setting electron beam pipe energy. If you don't know what to set here, leave it at 18
  Enable::HFARBWD_E_ENERGY = 18;

  // Beam Scattering configuration setting specified by CDR
  //
  // Option 1: ep-high-acceptance
  // Option 2: ep-high-divergence
  // Option 3: eA
  //
  // Enable::BEAM_COLLISION_SETTING = "ep-high-divergence";
  // If you don't know what to put here, set it to ep-high-divergence   
  //
  // Enable::BEAM_COLLISION_SETTING = "eA";
  Enable::BEAM_COLLISION_SETTING = "ep-high-divergence";

  // Either:
  // read previously generated g4-hits files, in this case it opens a DST and skips
  // the simulations step completely. The G4Setup macro is only loaded to get information
  // about the number of layers used for the cell reco code
  //
  //Input::READHITS = true;
  INPUTREADHITS::filename[0] = inputFile;
  // if you use a filelist
  // INPUTREADHITS::listfile[0] = inputFile;

  // Or:
  // Use one or more particle generators
  // It is run if Input::<generator> is set to true
  // all other options only play a role if it is active
  // In case embedding into a production output, please double check your G4Setup_EICDetector.C and G4_*.C consistent with those in the production macro folder
  //  Input::EMBED = true;
  INPUTEMBED::filename[0] = embed_input_file;
  // if you use a filelist
  //INPUTEMBED::listfile[0] = embed_input_file;

  // Use Pythia 8
  // Input::PYTHIA8 = true;

  // Use Pythia 6
  // Input::PYTHIA6 = true;

  // Use Sartre
  //   Input::SARTRE = true;

  // Simple multi particle generator in eta/phi/pt ranges
  Input::SIMPLE = true;
  // Input::SIMPLE_NUMBER = 2; // if you need 2 of them
  // Input::SIMPLE_VERBOSITY = 1;

  // Particle gun (same particles in always the same direction)
  // Input::GUN = true;
  // Input::GUN_NUMBER = 3; // if you need 3 of them
  // Input::GUN_VERBOSITY = 0;

  // Particle ion gun
  // Input::IONGUN = true; 

  // Upsilon generator
  // Input::UPSILON = true;
  // Input::UPSILON_NUMBER = 3; // if you need 3 of them
  // Input::UPSILON_VERBOSITY = 0;

  // And/Or read generated particles from file

  // eic-smear output
  // Input::READEIC = true;
  INPUTREADEIC::filename = inputFile;

  // HepMC2 files
  //  Input::HEPMC = true;
  Input::VERBOSITY = 0;
  INPUTHEPMC::filename = inputFile;

  //-----------------
  // Initialize the selected Input/Event generation
  //-----------------
  InputInit();
  //--------------
  // Set generator specific options
  //--------------
  // can only be set after InputInit() is called

  // Simple Input generator:
  // if you run more than one of these Input::SIMPLE_NUMBER > 1
  // add the settings for other with [1], next with [2]...
  if (Input::SIMPLE)
  {
    INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("pi-", 5);
    if (Input::HEPMC || Input::EMBED)
    {
      INPUTGENERATOR::SimpleEventGenerator[0]->set_reuse_existing_vertex(true);
      INPUTGENERATOR::SimpleEventGenerator[0]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
    }
    else
    {
      INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
                                                                                PHG4SimpleEventGenerator::Uniform,
                                                                                PHG4SimpleEventGenerator::Uniform);
      INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., 0.);
      INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0., 0., 5.);
    }
    INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(-3, 3);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(-M_PI, M_PI);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_pt_range(0.1, 20.);
  }
  // Upsilons
  // if you run more than one of these Input::UPSILON_NUMBER > 1
  // add the settings for other with [1], next with [2]...
  if (Input::UPSILON)
  {
    INPUTGENERATOR::VectorMesonGenerator[0]->add_decay_particles("mu", 0);
    INPUTGENERATOR::VectorMesonGenerator[0]->set_rapidity_range(-1, 1);
    INPUTGENERATOR::VectorMesonGenerator[0]->set_pt_range(0., 10.);
    // Y species - select only one, last one wins
    INPUTGENERATOR::VectorMesonGenerator[0]->set_upsilon_1s();
    if (Input::HEPMC || Input::EMBED)
    {
      INPUTGENERATOR::VectorMesonGenerator[0]->set_reuse_existing_vertex(true);
      INPUTGENERATOR::VectorMesonGenerator[0]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
    }
  }
  // particle gun
  // if you run more than one of these Input::GUN_NUMBER > 1
  // add the settings for other with [1], next with [2]...
  if (Input::GUN)
  {
    INPUTGENERATOR::Gun[0]->AddParticle("pi-", 0, 1, 0);
    INPUTGENERATOR::Gun[0]->set_vtx(0, 0, 0);
  }

  if (Input::IONGUN)
   {
     float theta = -25e-3;
 
     INPUTGENERATOR::IonGun[0]->SetA(197);
     INPUTGENERATOR::IonGun[0]->SetZ(79);
     INPUTGENERATOR::IonGun[0]->SetCharge(79);
     INPUTGENERATOR::IonGun[0]->SetMom(sin(theta)*110*197, 0,cos(theta)*110*197); // -25 mrad                        

     INPUTGENERATOR::IonGun[0]->Print();

   }

  // pythia6
  if (Input::PYTHIA6)
  {
    INPUTGENERATOR::Pythia6->set_config_file(string(getenv("CALIBRATIONROOT")) + "/Generators/phpythia6_ep.cfg");
    //! apply EIC beam parameter following EIC CDR
    Input::ApplyEICBeamParameter(INPUTGENERATOR::Pythia6);
  }
  // pythia8
  if (Input::PYTHIA8)
  {
    //! apply EIC beam parameter following EIC CDR
    Input::ApplyEICBeamParameter(INPUTGENERATOR::Pythia8);
  }
  // Sartre
  if (Input::SARTRE)
  {
    //! apply EIC beam parameter following EIC CDR
    Input::ApplyEICBeamParameter(INPUTGENERATOR::Sartre);
  }

  //--------------
  // Set Input Manager specific options
  //--------------
  // can only be set after InputInit() is called

  if (Input::HEPMC)
  {
    //! apply EIC beam parameter following EIC CDR
    Input::ApplyEICBeamParameter(INPUTMANAGER::HepMCInputManager);
    // optional overriding beam parameters
    //INPUTMANAGER::HepMCInputManager->set_vertex_distribution_width(100e-4, 100e-4, 30, 0);  //optional collision smear in space, time
    //    INPUTMANAGER::HepMCInputManager->set_vertex_distribution_mean(0,0,0,0);//optional collision central position shift in space, time
    // //optional choice of vertex distribution function in space, time
    // INPUTMANAGER::HepMCInputManager->set_vertex_distribution_function(PHHepMCGenHelper::Gaus, PHHepMCGenHelper::Gaus, PHHepMCGenHelper::Gaus, PHHepMCGenHelper::Gaus);
    //! embedding ID for the event
    //! positive ID is the embedded event of interest, e.g. jetty event from pythia
    //! negative IDs are backgrounds, .e.g out of time pile up collisions
    //! Usually, ID = 0 means the primary Au+Au collision background
    //INPUTMANAGER::HepMCInputManager->set_embedding_id(2);
  }

  // register all input generators with Fun4All
  InputRegister();

  // Reads event generators in EIC smear files, which is registered in InputRegister
  if (Input::READEIC)
  {
    //! apply EIC beam parameter following EIC CDR
    INPUTGENERATOR::EICFileReader->SetFirstEntry(skip);
    Input::ApplyEICBeamParameter(INPUTGENERATOR::EICFileReader);
  }

  // set up production relatedstuff
  //   Enable::PRODUCTION = true;

  //======================
  // Write the DST
  //======================

  // Enable::DSTOUT = true;
  DstOut::OutputDir = outdir;
  DstOut::OutputFile = outputFile;
  Enable::DSTOUT_COMPRESS = true;  // Compress DST files

  //Option to convert DST to human command readable TTree for quick poke around the outputs
  // Enable::DSTREADER = true;

  // turn the display on (default off)
  // Enable::DISPLAY = true;

  //======================
  // What to run
  //======================
  // Global options (enabled for all subsystems - if implemented)
  //  Enable::ABSORBER = true;
  //  Enable::OVERLAPCHECK = true;
  //  Enable::VERBOSITY = 1;

  // whether to simulate the Be section of the beam pipe
  Enable::PIPE = true;
  // If need to disable EIC beam pipe extension beyond the Be-section:
  G4PIPE::use_forward_pipes = true;
  //EIC hadron far forward magnets and detectors. IP6 and IP8 are incompatible (pick either or);
  Enable::HFARFWD_MAGNETS = true;
  Enable::HFARFWD_VIRTUAL_DETECTORS = true;

  Enable::HFARBWD_MAGNETS = true;
  Enable::HFARBWD_VIRTUAL_DETECTORS = true;

  // gems
  Enable::EGEM = false;
  Enable::FGEM = false; // deactivated as it's replaced by a FTTL layer
  // Enable::BGEM = true; // not yet defined in this model
  Enable::RWELL = true;
  // barrel tracker
  Enable::TrackingService = true;
  // Enable::TrackingService_VERBOSITY = INT_MAX - 10;
  Enable::BARREL = true;
  // fst
  Enable::FST = true;

  //AC-LGAD  TOFs
  Enable::FTTL = true;
  Enable::ETTL = true;
  Enable::CTTL = true;

  //mRPC TOFs
  Enable::BTOF = false;
  Enable::ETOF = false;
  Enable::HTOF = false;
  Enable::ETOF_GAS = Enable::ETOF && true;
  Enable::HTOF_GAS = Enable::HTOF && true;

  Enable::TRACKING = true;
  Enable::TRACKING_EVAL = Enable::TRACKING && true;
  G4TRACKING::DISPLACED_VERTEX = true;  // this option exclude vertex in the track fitting and use RAVE to reconstruct primary and 2ndary vertexes
                                        // projections to calorimeters
  G4TRACKING::PROJECTION_EEMC = true;
  G4TRACKING::PROJECTION_BECAL = true;
  G4TRACKING::PROJECTION_EHCAL = true;
  G4TRACKING::PROJECTION_CEMC = true;
  G4TRACKING::PROJECTION_HCALIN = true;
  G4TRACKING::PROJECTION_HCALOUT = true;
  G4TRACKING::PROJECTION_FEMC = true;
  G4TRACKING::PROJECTION_FHCAL = true;
  G4TRACKING::PROJECTION_LFHCAL = true;

  // enable barrel calos & magnet
  Enable::BECAL   = true;
  Enable::HCALIN  = true;
  Enable::MAGNET  = true;
  Enable::HCALOUT = true;

  // EICDetector geometry - barrel
  Enable::DIRC = true;
  Enable::DIRC_RECO = Enable::DIRC && true;

  Enable::BMMG = false;
  // Enable::DIRC_VERBOSITY = 2;

  // EICDetector geometry - 'hadron' direction
  Enable::RICH = true;
  Enable::RICH_RECO = Enable::RICH && true;

  Enable::TRD = false;
  Enable::TRD_GAS = false;
  // Enable::RICH_VERBOSITY = 2;


  // enable forward calos
  Enable::FEMC    = true;
  Enable::DRCALO  = false;
  Enable::LFHCAL  = true;

  // EICDetector geometry - 'electron' direction
  Enable::mRICH = true;
  Enable::mRICH_RECO = Enable::mRICH && true;
  // Enable::mRICH_VERBOSITY = 2;
  
  // EICDetector geometry - 'electron' direction
  Enable::EEMCH   = true;
  G4EEMCH::SETTING::USECUSTOMMAPUPDATED = true; // enable proper carbon structure
  G4TTL::SETTING::optionEEMCH           = Enable::EEMCH && true;
  Enable::EHCAL   = false;

  Enable::FFR_EVAL = Enable::HFARFWD_MAGNETS && Enable::HFARFWD_VIRTUAL_DETECTORS && true;

  Enable::PLUGDOOR = false;

  // Other options
  Enable::GLOBAL_RECO = G4TRACKING::DISPLACED_VERTEX;  // use reco vertex for global event vertex
  Enable::GLOBAL_FASTSIM = true;

  // jet reconstruction
  Enable::FWDJETS = true;
  Enable::FWDJETS_EVAL = Enable::FWDJETS && false;

  // new settings using Enable namespace in GlobalVariables.C
  Enable::BLACKHOLE = true;
  //Enable::BLACKHOLE_SAVEHITS = false; // turn off saving of bh hits
  //BlackHoleGeometry::visible = true;

  // ZDC
  // Enable::ZDC = true;
  // Enable::ZDC_DISABLE_BLACKHOLE = true;

  // B0
  // B0 shape
  // Enable::B0_DISABLE_HITPLANE = true;
  // Enable::B0_FULLHITPLANE = true;
  // Enable::B0_VAR_PIPE_HOLE = true;
  // Enable::B0_CIRCLE_PIPE_HOLE = false;
  
  // B0 Tracking
  // Enable::B0TRACKING = false; // For B0 Tracking
  // Enable::B0TRACKING_EVAL = Enable::B0TRACKING && true; //For B0 Tracking

  // B0 ECAL
  // Enable::B0ECALTOWERS = true;  //To Construct Towers of B0ECal instead of one single volume
  // Enable::B0ECAL = Enable::B0_DISABLE_HITPLANE && true;
  // Enable::B0ECAL_CELL = Enable::B0ECAL && true;
  // Enable::B0ECAL_TOWER = Enable::B0ECAL_CELL && true;
  // Enable::B0ECAL_CLUSTER = Enable::B0ECAL_TOWER && true;
  // Enable::B0ECAL_EVAL = Enable::B0ECAL_CLUSTER && true;
    
  // RP
  // Enable::RP_DISABLE_HITPLANE = true;
   
   //Far Backward detectors
//  Enable::BWD = true;
//  Enable::BWDN[0]=true; // 1st Q2 tagger
//  Enable::BWDN[1]=true; // 2nd Q2 tagger
//  Enable::BWDN[2]=true; // 1st Lumi
//  Enable::BWDN[3]=true; // 2nd Lumi (+)
//  Enable::BWDN[4]=true; // 3rd Lumi (-)
//  Enable::BWD_CELL = Enable::BWD && true;
//  Enable::BWD_TOWER = Enable::BWD_CELL && true;
//  Enable::BWD_CLUSTER = Enable::BWD_TOWER && true;
//  Enable::BWD_EVAL = Enable::BWD_CLUSTER && true;

  //************************************************************
  // details for calos: cells, towers, clusters
  //************************************************************
  Enable::BECAL_CELL    = Enable::BECAL && true;
  Enable::BECAL_TOWER   = Enable::BECAL_CELL && true;
  Enable::BECAL_CLUSTER = Enable::BECAL_TOWER && true;
  Enable::BECAL_EVAL    = Enable::BECAL_CLUSTER && true;

  //  Enable::HCALIN_ABSORBER = true;
  Enable::HCALIN_CELL     = Enable::HCALIN && true;
  Enable::HCALIN_TOWER    = Enable::HCALIN_CELL && true;
  Enable::HCALIN_CLUSTER  = Enable::HCALIN_TOWER && true;
  Enable::HCALIN_EVAL     = Enable::HCALIN_CLUSTER && true;

  //  Enable::HCALOUT_ABSORBER = true;
  Enable::HCALOUT_CELL    = Enable::HCALOUT && true;
  Enable::HCALOUT_TOWER   = Enable::HCALOUT_CELL && true;
  Enable::HCALOUT_CLUSTER = Enable::HCALOUT_TOWER && true;
  Enable::HCALOUT_EVAL    = Enable::HCALOUT_CLUSTER && true;
  
  //  Enable::FEMC_ABSORBER = true;
  Enable::FEMC_TOWER      = Enable::FEMC && true;
  Enable::FEMC_CLUSTER    = Enable::FEMC_TOWER && true;
  Enable::FEMC_EVAL       = Enable::FEMC_CLUSTER && true;
  
  Enable::DRCALO_CELL     = Enable::DRCALO && true;
  Enable::DRCALO_TOWER    = Enable::DRCALO_CELL && true;
  Enable::DRCALO_CLUSTER  = Enable::DRCALO_TOWER && true;
  Enable::DRCALO_EVAL     = Enable::DRCALO_CLUSTER && false;

  Enable::LFHCAL_ABSORBER = false;
  Enable::LFHCAL_CELL     = Enable::LFHCAL && true;
  Enable::LFHCAL_TOWER    = Enable::LFHCAL_CELL && true;
  Enable::LFHCAL_CLUSTER  = Enable::LFHCAL_TOWER && true;
  Enable::LFHCAL_EVAL     = Enable::LFHCAL_CLUSTER && true;

  Enable::EEMCH_TOWER     = Enable::EEMCH && true;
  Enable::EEMCH_CLUSTER   = Enable::EEMCH_TOWER && true;
  Enable::EEMCH_EVAL      = Enable::EEMCH_CLUSTER && true;

  Enable::EHCAL_CELL      = Enable::EHCAL && true;
  Enable::EHCAL_TOWER     = Enable::EHCAL_CELL && true;
  Enable::EHCAL_CLUSTER   = Enable::EHCAL_TOWER && true;
  Enable::EHCAL_EVAL      = Enable::EHCAL_CLUSTER && true;
  
  // Enabling the event evaluator?
  Enable::EVENT_EVAL            = true;
  Enable::EVENT_EVAL_DO_HITS    = true;
  Enable::EVENT_EVAL_DO_HEPMC   = Input::PYTHIA6 or Input::PYTHIA8 or Input::SARTRE or Input::HEPMC or Input::READEIC;
  Enable::EVENT_EVAL_DO_EVT_LVL = Input::PYTHIA6 or Input::PYTHIA8 or Input::READEIC;

  //Enable::USER = true;

  //---------------
  // World Settings
  //---------------
  //  G4WORLD::PhysicsList = "FTFP_BERT"; //FTFP_BERT_HP best for calo
  //  G4WORLD::WorldMaterial = "G4_AIR"; // set to G4_GALACTIC for material scans
  //  G4WORLD::WorldMaterial = "G4_Galactic"; // set to G4_GALACTIC for material scans

  //---------------
  // Magnet Settings
  //---------------

  //  const string magfield = "1.5"; // alternatively to specify a constant magnetic field, give a float number, which will be translated to solenoidal field in T, if string use as fieldmap name (including path)
  //  G4MAGNET::magfield = string(getenv("CALIBRATIONROOT")) + string("/Field/Map/sPHENIX.2d.root");  // default map from the calibration database
  G4MAGNET::magfield_rescale = -1.4 / 1.5;  // make consistent with expected Babar field strength of 1.4T

  //---------------
  // Pythia Decayer
  //---------------
  // list of decay types in
  // $OFFLINE_MAIN/include/g4decayer/EDecayType.hh
  // default is All:
  // G4P6DECAYER::decayType = EDecayType::kAll;

  // Initialize the selected subsystems
  G4Init();

  //---------------------
  // GEANT4 Detector description
  //---------------------

  // If "readhepMC" is also set, the Upsilons will be embedded in Hijing events, if 'particles" is set, the Upsilons will be embedded in whatever particles are thrown
  if (!Input::READHITS)
  {
    G4Setup();
  }

  //------------------
  // Detector Division
  //------------------
  if (Enable::CEMC_CELL) CEMC_Cells();

  if (Enable::HCALIN_CELL) HCALInner_Cells();

  if (Enable::HCALOUT_CELL) HCALOuter_Cells();

  //-----------------------------
  // CEMC towering and clustering
  //-----------------------------

  if (Enable::CEMC_TOWER) CEMC_Towers();
  if (Enable::CEMC_CLUSTER) CEMC_Clusters();

  //-----------------------------
  // HCAL towering and clustering
  //-----------------------------

  if (Enable::HCALIN_TOWER) HCALInner_Towers();
  if (Enable::HCALIN_CLUSTER) HCALInner_Clusters();

  if (Enable::HCALOUT_TOWER) HCALOuter_Towers();
  if (Enable::HCALOUT_CLUSTER) HCALOuter_Clusters();

  //-----------------------------
  // e, h direction Calorimeter  towering and clustering
  //-----------------------------

  if (Enable::FEMC_TOWER) FEMC_Towers();
  if (Enable::FEMC_CLUSTER) FEMC_Clusters();

  if (Enable::FHCAL_TOWER) FHCAL_Towers();
  if (Enable::FHCAL_CLUSTER) FHCAL_Clusters();

  if (Enable::DRCALO_TOWER) DRCALO_Towers();
  if (Enable::DRCALO_CLUSTER) DRCALO_Clusters();

  if (Enable::LFHCAL_TOWER) LFHCAL_Towers();
  if (Enable::LFHCAL_CLUSTER) LFHCAL_Clusters();

  if (Enable::EEMC_TOWER) EEMC_Towers();
  if (Enable::EEMC_CLUSTER) EEMC_Clusters();

  if (Enable::EEMCH_TOWER) EEMCH_Towers();
  if (Enable::EEMCH_CLUSTER) EEMCH_Clusters();

  if (Enable::EHCAL_TOWER) EHCAL_Towers();
  if (Enable::EHCAL_CLUSTER) EHCAL_Clusters();

  if (Enable::BECAL_TOWER) BECAL_Towers();
  if (Enable::BECAL_CLUSTER) BECAL_Clusters();
    
  if (Enable::B0ECAL_TOWER) B0ECAL_Towers(); // For B0Ecal
  if (Enable::B0ECAL_CLUSTER) B0ECAL_Clusters(); //For B0Ecal
    
  if (Enable::BWD_TOWER) BWD_Towers(); // For Bwd
  if (Enable::BWD_CLUSTER) BWD_Clusters(); //For Bwd

  if (Enable::DSTOUT_COMPRESS) ShowerCompress();

  //--------------
  // Tracking and PID
  //--------------

  if (Enable::TRACKING) Tracking_Reco();

  if (Enable::DIRC_RECO) DIRCReco();

  if (Enable::mRICH_RECO ) mRICHReco();

  if (Enable::RICH_RECO) RICHReco();

  //-----------------
  // Global Vertexing
  //-----------------

  if (Enable::GLOBAL_RECO)
  {
    Global_Reco();
  }
  else if (Enable::GLOBAL_FASTSIM)
  {
    Global_FastSim();
  }

  //---------
  // Jet reco
  //---------

  if (Enable::FWDJETS) Jet_FwdReco();

  string outputroot = outputFile;
  string remove_this = ".root";
  size_t pos = outputroot.find(remove_this);
  if (pos != string::npos)
  {
    outputroot.erase(pos, remove_this.length());
  }

  if (Enable::DSTREADER) G4DSTreader_EICDetector(outputroot + "_DSTReader.root");

  //----------------------
  // Simulation evaluation
  //----------------------

  if (Enable::EVENT_EVAL) Event_Eval(outputroot + "_eventtree.root");

  if (Enable::TRACKING_EVAL) Tracking_Eval(outputroot + "_g4tracking_eval.root");

  if (Enable::CEMC_EVAL) CEMC_Eval(outputroot + "_g4cemc_eval.root");

  if (Enable::HCALIN_EVAL) HCALInner_Eval(outputroot + "_g4hcalin_eval.root");

  if (Enable::HCALOUT_EVAL) HCALOuter_Eval(outputroot + "_g4hcalout_eval.root");

  if (Enable::FEMC_EVAL) FEMC_Eval(outputroot + "_g4femc_eval.root");

  if (Enable::FHCAL_EVAL) FHCAL_Eval(outputroot + "_g4fhcal_eval.root");

  if (Enable::EEMC_EVAL) EEMC_Eval(outputroot + "_g4eemc_eval.root");

  if (Enable::FFR_EVAL) FFR_Eval(outputroot + "_g4ffr_eval.root");

  if (Enable::B0ECAL_EVAL) B0ECAL_Eval(outputroot + "_g4b0ecal_eval_test.root"); // For B0Ecal
    
  if (Enable::BWD_EVAL) BWD_Eval(outputroot + "_g4bwd_eval_e0100_debug"); // For Bwd
    
  if (Enable::FWDJETS_EVAL) Jet_FwdEval();

  if (Enable::USER) UserAnalysisInit();

  //--------------
  // Set up Input Managers
  //--------------

  InputManagers();

  //--------------
  // Set up Output Manager
  //--------------
  if (Enable::PRODUCTION)
  {
    Production_CreateOutputDir();
  }

  if (Enable::DSTOUT)
  {
    string FullOutFile = DstOut::OutputDir + "/" + DstOut::OutputFile;
    Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", FullOutFile);
    if (Enable::DSTOUT_COMPRESS) DstCompress(out);
    se->registerOutputManager(out);
  }

  //-----------------
  // Event processing
  //-----------------
  if (Enable::DISPLAY)
  {
    DisplayOn();

    gROOT->ProcessLine("Fun4AllServer *se = Fun4AllServer::instance();");
    gROOT->ProcessLine("PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");");

    cout << "-------------------------------------------------" << endl;
    cout << "You are in event display mode. Run one event with" << endl;
    cout << "se->run(1)" << endl;
    cout << "Run Geant4 command with following examples" << endl;
    gROOT->ProcessLine("displaycmd()");

    return 0;
  }
  // if we use a negative number of events we go back to the command line here
  if (nEvents < 0)
  {
    return 0;
  }
  // if we run any of the particle generators and use 0 it'll run forever
  if (nEvents == 0 && !Input::READHITS && !Input::HEPMC && !Input::READEIC)
  {
    cout << "using 0 for number of events is a bad idea when using particle generators" << endl;
    cout << "it will run forever, so I just return without running anything" << endl;
    return 0;
  }

  se->skip(skip);
  se->run(nEvents);

  //-----
  // Exit
  //-----

  se->End();
  std::cout << "All done" << std::endl;
  delete se;
  if (Enable::PRODUCTION)
  {
    Production_MoveOutput();
  }
  gSystem->Exit(0);
  return 0;
}
#endif
