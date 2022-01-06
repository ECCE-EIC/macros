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

int Fun4All_G4_ECCEModular(
    const int nEvents               = 1,
    const double particlemomMin     = -1,
    const double particlemomMax     = -1,
    const string detectorSettings   = "TTLGEO_5",
    const TString generatorSettings = "PYTHIA8",
    const string &inputFile         = "https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/sPHENIX_G4Hits_sHijing_9-11fm_00000_00010.root",
    const string &outputFile        = "G4EICDetector.root",
    const string &embed_input_file  = "https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/sPHENIX_G4Hits_sHijing_9-11fm_00000_00010.root",
    const int skip                  = 0,
    const string &outdir            = ".")
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
  // Enable::IP8 = true;

  // Setting proton beam pipe energy. If you don't know what to set here, leave it at 275
  Enable::HFARFWD_ION_ENERGY = 275;

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

  if(particlemomMin==-1 && particlemomMax==-1){
    if (generatorSettings.Contains("PYTHIA6")) {
      Input::PYTHIA6  = true;
    } else if (generatorSettings.Contains("PYTHIA8")) {
      Input::PYTHIA8  = true;
    } else if (generatorSettings.Contains("SATRE")) {
      Input::SARTRE   = true;
    } else if (generatorSettings.Contains("READEIC")) {
      Input::READEIC  = true;
    } else if (generatorSettings.Contains("HEPMCINPUT")) {
      Input::HEPMC    = true;
      INPUTHEPMC::filename = inputFile;
    }
  }
  // Simple multi particle generator in eta/phi/pt ranges
  Input::SIMPLE = false;
  if (particlemomMin>-1 && particlemomMax>-1){
    Input::SIMPLE = true;
    Input::SIMPLE_VERBOSITY = 0;
    if (generatorSettings.Contains("Multi"))
      Input::SIMPLE_NUMBER = 3; // if you need 2 of them
  }

  Input::VERBOSITY = 0;
  
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
  if (Input::SIMPLE){
      if (generatorSettings.Contains("Multi")){
        for(int igen=0;igen<Input::SIMPLE_NUMBER;igen++){
          if (generatorSettings.Contains("PiPrEl")){
            if(igen==0)INPUTGENERATOR::SimpleEventGenerator[igen]->add_particles("pi-", 1);
            else if(igen==1)INPUTGENERATOR::SimpleEventGenerator[igen]->add_particles("e-", 1);
            else if(igen==2)INPUTGENERATOR::SimpleEventGenerator[igen]->add_particles("proton", 1);
          }else if (generatorSettings.Contains("Pion"))
            INPUTGENERATOR::SimpleEventGenerator[igen]->add_particles("pi-", 1);
          else if (generatorSettings.Contains("Kaon"))
            INPUTGENERATOR::SimpleEventGenerator[igen]->add_particles("kaon-", 1);
          else if (generatorSettings.Contains("Proton"))
            INPUTGENERATOR::SimpleEventGenerator[igen]->add_particles("proton", 1);
          else if (generatorSettings.Contains("Muon"))
            INPUTGENERATOR::SimpleEventGenerator[igen]->add_particles("mu-", 1);
          else if (generatorSettings.Contains("Photon"))
            INPUTGENERATOR::SimpleEventGenerator[igen]->add_particles("gamma", 1);
          else if (generatorSettings.Contains("Neutron"))
            INPUTGENERATOR::SimpleEventGenerator[igen]->add_particles("neutron", 1);
          else if (generatorSettings.Contains("Lambda"))
            INPUTGENERATOR::SimpleEventGenerator[igen]->add_particles("lambda", 1);
          else if (generatorSettings.Contains("K0S"))
            INPUTGENERATOR::SimpleEventGenerator[igen]->add_particles("kaon0S", 1);
          else if (generatorSettings.Contains("Electron"))
            INPUTGENERATOR::SimpleEventGenerator[igen]->add_particles("e-", 1);
          else if (generatorSettings.Contains("PiZero"))
            INPUTGENERATOR::SimpleEventGenerator[igen]->add_particles("pi0", 1);
          else {
            std::cout << "You didn't specify which particle you wanted to generate, exiting" << std::endl;
            return 0;
          }
          INPUTGENERATOR::SimpleEventGenerator[igen]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
                                                                                    PHG4SimpleEventGenerator::Uniform,
                                                                                    PHG4SimpleEventGenerator::Uniform);
          INPUTGENERATOR::SimpleEventGenerator[igen]->set_vertex_distribution_mean(0., 0., 0.);
          INPUTGENERATOR::SimpleEventGenerator[igen]->set_vertex_distribution_width(0., 0., 0.);

          bool strictrange = false;
          if(generatorSettings.Contains("strict")) strictrange = true;
          if (igen==0){
            if(strictrange)
              INPUTGENERATOR::SimpleEventGenerator[igen]->set_eta_range(-0.7, 0.7);
            else
              INPUTGENERATOR::SimpleEventGenerator[igen]->set_eta_range(-1.7, 1.2);
          } else if (igen==1) {
            if(strictrange)
              INPUTGENERATOR::SimpleEventGenerator[igen]->set_eta_range(-3.0, -2.5);
            else
              INPUTGENERATOR::SimpleEventGenerator[igen]->set_eta_range(-4, -1.7);
          } else if (igen==2) {
            if(strictrange)
              INPUTGENERATOR::SimpleEventGenerator[igen]->set_eta_range(2.5, 3.0);
            else
              INPUTGENERATOR::SimpleEventGenerator[igen]->set_eta_range(1.2, 4.0);
          } else {
            INPUTGENERATOR::SimpleEventGenerator[igen]->set_eta_range(-4.0, 4.0);
          }
          INPUTGENERATOR::SimpleEventGenerator[igen]->set_phi_range(-M_PI, M_PI);
          INPUTGENERATOR::SimpleEventGenerator[igen]->set_p_range(particlemomMin, particlemomMax);
        }
    } else {
      if (generatorSettings.Contains("SimplePion"))
        INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("pi-", 1);
      else if (generatorSettings.Contains("SimpleKaon"))
        INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("kaon-", 1);
      else if (generatorSettings.Contains("SimpleProton"))
        INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("proton", 1);
      else if (generatorSettings.Contains("SimplePhoton"))
        INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("gamma", 1);
      else if (generatorSettings.Contains("SimpleNeutron"))
        INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("neutron", 1);
      else if (generatorSettings.Contains("SimpleMuon"))
        INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("mu-", 1);
      else if (generatorSettings.Contains("SimpleLambda"))
        INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("lambda", 1);
      else if (generatorSettings.Contains("SimpleK0S"))
        INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("kaon0S", 1);
      else if (generatorSettings.Contains("SimpleElectron"))
        INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("e-", 1);
      else if (generatorSettings.Contains("SimplePiZero"))
        INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("pi0", 1);
      else {
        std::cout << "You didn't specify which particle you wanted to generate, exiting" << std::endl;
        return 0;
      }
      INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
                                                                                PHG4SimpleEventGenerator::Uniform,
                                                                                PHG4SimpleEventGenerator::Uniform);
      INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., 0.);
      INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0., 0., 0.);
      if (generatorSettings.Contains("central"))
        INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(-1.8, 1.2);
      else if (generatorSettings.Contains("bck"))
        INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(-4, -1.7);
      else if (generatorSettings.Contains("fwd"))
        INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(1.2, 4.0);
      else
        INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(-4.0, 4.0);
      INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(-M_PI, M_PI);
      INPUTGENERATOR::SimpleEventGenerator[0]->set_p_range(particlemomMin, particlemomMax);
    }
  }
  if(particlemomMin>-1 && particlemomMax == -1){
    PHG4ParticleGenerator *gen = new PHG4ParticleGenerator("PGENERATOR");
    gen->set_name("pi-");
    // gen->set_name("pi0");
    gen->set_vtx(0, 0, 0);
    gen->set_eta_range(-4.0, 4.0);            // around midrapidity
    if(particlemomMin > -1) {
      gen->set_mom_range(particlemomMin, particlemomMin);                   // fixed 4 GeV/c
    }
    else {
      gen->set_mom_range(1, 60);                   // fixed 4 GeV/c
    }
    gen->set_phi_range(0., 2* M_PI);  // 0-90 deg
    // gen->Verbosity(1);  // 0-90 deg
    se->registerSubsystem(gen);
  }
  // pythia6
  if (Input::PYTHIA6)
  {
    //INPUTGENERATOR::Pythia6->set_config_file(string(getenv("CALIBRATIONROOT")) + "/Generators/phpythia6_ep.cfg");
    INPUTGENERATOR::Pythia6->set_config_file(inputFile);
    //! apply EIC beam parameter following EIC CDR
    Input::ApplyEICBeamParameter(INPUTGENERATOR::Pythia6);
  }
  // pythia8
  if (Input::PYTHIA8)
  {
    // Configuration file
    PYTHIA8::config_file = inputFile;
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
    INPUTREADEIC::filename = inputFile;
    //! apply EIC beam parameter following EIC CDR
    INPUTGENERATOR::EICFileReader->SetFirstEntry(skip);
    Input::ApplyEICBeamParameter(INPUTGENERATOR::EICFileReader);
  }

  // set up production relatedstuff
  //   Enable::PRODUCTION = true;

  //======================
  // Write the DST
  //======================

  Enable::DSTOUT = false;
  DstOut::OutputDir = outdir;
  DstOut::OutputFile = outputFile;
  Enable::DSTOUT_COMPRESS = false;  // Compress DST files

  //Option to convert DST to human command readable TTree for quick poke around the outputs
  // Enable::DSTREADER = true;

  // turn the display on (default off)
  if (detectorSettings.find("display") != std::string::npos) {
    Enable::DISPLAY = true;
  }
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
  // Enable::HFARFWD_MAGNETS = true;
  // Enable::HFARFWD_VIRTUAL_DETECTORS = true;

  // Enable::HFARBWD_MAGNETS = true;
  // Enable::HFARBWD_VIRTUAL_DETECTORS = true;

  //***********************************************
  // barrel trackers
  //***********************************************
  Enable::RWELL = true;
  // barrel tracker
  Enable::TrackingService = true;
  // Enable::TrackingService_VERBOSITY = INT_MAX - 10;
  Enable::BARREL = true;
  // fst
  Enable::FST = true;

  //***********************************************
  // TOFs
  //***********************************************
  Enable::FTTL          = true;
  Enable::ETTL          = true;
  Enable::CTTL          = true;
  G4TTL::SETTING::optionCEMC = false;
  std::string ttlSettingToFind = "TTLGEO_";
  if (detectorSettings.find(ttlSettingToFind) != std::string::npos) {
    auto pos = detectorSettings.find(ttlSettingToFind);
    G4TTL::SETTING::optionGeo = std::stoi(detectorSettings.substr(pos + ttlSettingToFind.size(), pos + ttlSettingToFind.size() + 1));
  }  else {
    G4TTL::SETTING::optionGeo = 1;
  }

  //***********************************************
  // gems fwd & bwd
  //***********************************************
  Enable::EGEM = false;
  Enable::FGEM = false;

  //***********************************************
  // tracking macro settings
  //***********************************************
  Enable::TRACKING              = true;
  Enable::TRACKING_EVAL         = Enable::TRACKING && false;
  G4TRACKING::DISPLACED_VERTEX  = true;  // this option exclude vertex in the track fitting and use RAVE to reconstruct primary and 2ndary vertexes
                                        // projections to calorimeters
  G4TRACKING::PROJECTION_EEMC     = true;
  G4TRACKING::PROJECTION_BECAL    = true;
  G4TRACKING::PROJECTION_EHCAL    = true;
  G4TRACKING::PROJECTION_CEMC     = true;
  G4TRACKING::PROJECTION_HCALIN   = true;
  G4TRACKING::PROJECTION_HCALOUT  = true;
  G4TRACKING::PROJECTION_FEMC     = true;
  G4TRACKING::PROJECTION_FHCAL    = true;
  G4TRACKING::PROJECTION_LFHCAL   = true;

  //***********************************************
  // barrel calos & magnet
  //***********************************************
  Enable::BECAL   = true;
  Enable::HCALIN  = true;
  Enable::MAGNET  = true;
  Enable::HCALOUT = true;

  //***********************************************
  // cherenkov's
  //***********************************************
  Enable::DIRC    = true;
  Enable::RICH    = true;
  Enable::mRICH   = true;

  //***********************************************
  // fwd calos
  //***********************************************
  Enable::FEMC    = true;
  Enable::DRCALO  = false;
  G4TTL::SETTING::optionDR = 1;
  Enable::LFHCAL  = true;

  //***********************************************
  // bwd calos
  //***********************************************
  Enable::EEMCH = true;
  if (detectorSettings.find("EEMAPNC") != std::string::npos) {
    G4EEMCH::SETTING::USECUSTOMMAPNOCARBON = true;
  }
  if (detectorSettings.find("EEMAP30CM") != std::string::npos) {
    G4EEMCH::SETTING::USECUSTOMMAP30CM = true;
  }
  if (detectorSettings.find("EEMAPCARBON") != std::string::npos) {
    G4EEMCH::SETTING::USECUSTOMMAPCARBON = true;
  }
  if (detectorSettings.find("EEMAPUPDATE") != std::string::npos) {
    G4EEMCH::SETTING::USECUSTOMMAPUPDATED = true;
  }
  Enable::EHCAL     = false;
  Enable::PLUGDOOR  = true;

  // Other options
  Enable::GLOBAL_RECO = G4TRACKING::DISPLACED_VERTEX;  // use reco vertex for global event vertex
  Enable::GLOBAL_FASTSIM = true;

  // jet reconstruction
  Enable::FWDJETS       = false;
  Enable::FWDJETS_EVAL  = Enable::FWDJETS && false;

  // new settings using Enable namespace in GlobalVariables.C
  Enable::BLACKHOLE = true;
  bool BLACKHOLE_SAVEHITS = false;
  if(detectorSettings.find("BHH")!= std::string::npos ){
    Enable::BLACKHOLE_SAVEHITS = true; // turn off saving of bh hits
    Enable::EVENT_EVAL_DO_HITS_BLACKHOLE = true; // turn off saving of bh hits
  }
  // BlackHoleGeometry::visible = true;
  
  // ZDC
  // Enable::ZDC = true;
  // Enable::ZDC_DISABLE_BLACKHOLE = true;

  // B0
  // Enable::B0_DISABLE_HITPLANE = true;
  // Enable::B0_FULLHITPLANE = true;

  // Enable::B0ECALTOWERS = true;  //To Construct Towers of B0ECal instead of one single volume
  // Enable::B0ECAL = Enable::B0_DISABLE_HITPLANE && true;
  // Enable::B0ECAL_CELL = Enable::B0ECAL && true;
  // Enable::B0ECAL_TOWER = Enable::B0ECAL_CELL && true;
  // Enable::B0ECAL_CLUSTER = Enable::B0ECAL_TOWER && true;
  // Enable::B0ECAL_EVAL = Enable::B0ECAL_CLUSTER && true;
    
  // RP
  // Enable::RP_DISABLE_HITPLANE = true;
  // Enable::RP_FULLHITPLANE = true;

  // RP after 2nd focus for IP8 only
  // Enable::RP2nd_DISABLE_HITPLANE = true;
  // Enable::RP2nd_FULLHITPLANE = true;

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // special settings for Calo standalone studies
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // deactivate all respective detector systems for standalone studies
  if(detectorSettings.find("STANDALONE")!= std::string::npos){
    Enable::PIPE = false;
    G4PIPE::use_forward_pipes = false;
    Enable::HFARFWD_MAGNETS = false;
    Enable::HFARFWD_VIRTUAL_DETECTORS = false;
    // Enable::TPC_ENDCAP = false;
    // G4TRACKING::PROJECTION_CEMC   = false;
    // G4TRACKING::PROJECTION_FEMC   = false;
    // G4TRACKING::PROJECTION_FHCAL  = false;
    // G4TRACKING::PROJECTION_EHCAL  = false;
    // G4TRACKING::PROJECTION_DRCALO = false;
    // G4TRACKING::PROJECTION_EEMC   = false;
    Enable::MAGNET = false;
    Enable::DIRC = false;
    Enable::RICH = false;
    Enable::mRICH = false;
    // Enable::AEROGEL = false;
    Enable::CEMC = false;
    Enable::HCALOUT = false;
    Enable::HCALIN = false;
    Enable::EHCAL = false;
    Enable::EEMC = false;
    Enable::EEMCH = false;
    Enable::FEMC = false;
    Enable::FHCAL = false;
    Enable::LFHCAL = false;
    Enable::BECAL = false;
    Enable::FTTL = false;
    Enable::CTTL = false;
    Enable::ETTL = false;
    Enable::EEMCH = false;
    Enable::RWELL = false;
    Enable::TrackingService = false;
    Enable::BARREL = false;
    Enable::FST = false;
    Enable::EGEM = false;
    Enable::FGEM = false;
    if(detectorSettings.find("PIPE")!= std::string::npos ){
      Enable::PIPE = true;
      G4PIPE::use_forward_pipes = true;
    }
    if(detectorSettings.find("Magnet")!= std::string::npos )
      Enable::MAGNET = true;
    if(detectorSettings.find("dRICH")!= std::string::npos )
      Enable::RICH = true;
    // if(detectorSettings.find("ALLSILICON")!= std::string::npos )
    //   Enable::ALLSILICON = true;
    if(detectorSettings.find("CEMC")!= std::string::npos )
      Enable::CEMC = true;
    if(detectorSettings.find("HCALOUT")!= std::string::npos ){
      Enable::HCALOUT = true;
    }
    if(detectorSettings.find("HCALIN")!= std::string::npos ){
      Enable::HCALIN = true;
    }
    if(detectorSettings.find("HCALINOUT")!= std::string::npos ){
      Enable::HCALOUT = true;
      Enable::HCALIN = true;
    }
    if(detectorSettings.find("DR")!= std::string::npos )
      Enable::DRCALO = true;
    if(detectorSettings.find("FEMC")!= std::string::npos )
      Enable::FEMC = true;
    if(detectorSettings.find("FGEM")!= std::string::npos )
      Enable::FGEM = true;
    if((detectorSettings.find("FHCAL")!= std::string::npos) && !(detectorSettings.find("LFHCAL")!= std::string::npos) )
      Enable::FHCAL = true;
    if(detectorSettings.find("LFHCAL")!= std::string::npos )
      Enable::LFHCAL = true;
    if(detectorSettings.find("BECAL")!= std::string::npos )
      Enable::BECAL = true;
    if(detectorSettings.find("EHCAL")!= std::string::npos )
      Enable::EHCAL = true;
    if(detectorSettings.find("EEMCH")!= std::string::npos )
      Enable::EEMCH = true;
    if(detectorSettings.find("RWELL")!= std::string::npos )
      Enable::RWELL = true;
    if(detectorSettings.find("CHCAL")!= std::string::npos ){
      Enable::HCALIN   = true;
      Enable::HCALOUT  = true;
    }
    if(detectorSettings.find("DIRC")!= std::string::npos )
      Enable::DIRC = true;
    if(detectorSettings.find("SUPPORT")!= std::string::npos ){
      Enable::TrackingService = true;
      Enable::TrackingService_OVERLAPCHECK = true;
    }
    if(detectorSettings.find("FWDCALO")!= std::string::npos ){
      Enable::FEMC    = true;
      Enable::FHCAL   = true;
    }
    if(detectorSettings.find("FWDLCALO")!= std::string::npos ){
      Enable::FEMC    = true;
      Enable::LFHCAL  = true;
    }
    if(detectorSettings.find("BARCALO")!= std::string::npos ){
      Enable::BECAL    = true;
      Enable::HCALIN   = true;
      Enable::HCALOUT  = true;
      Enable::MAGNET = true;
    }
    if(detectorSettings.find("BCKCALO")!= std::string::npos ){
      Enable::EHCAL    = true;
      Enable::EEMCH    = true;
    }
    if(detectorSettings.find("TTL")!= std::string::npos ){
      // Enable::PIPE = true;
      // G4PIPE::use_forward_pipes = true;
      // LGAD layers
      if(detectorSettings.find("FTTL")!= std::string::npos )
        Enable::FTTL = true;
      if(detectorSettings.find("ETTL")!= std::string::npos ){
        Enable::ETTL = true;
        // G4DIRC::SETTING::USECEMCGeo   = false;
        G4TTL::SETTING::optionCEMC    = false;
      }
      if(detectorSettings.find("CTTL")!= std::string::npos ){
        Enable::CTTL = true;
        // Enable::DIRC = true;
        // Enable::CEMC = true;
        // Enable::BECAL = true;
        // Enable::ALLSILICON = true;
        // G4DIRC::SETTING::USECEMCGeo   = false;
        G4TTL::SETTING::optionCEMC    = false;
      }
    }
  }

  //************************************************************
  // details for calos: cells, towers, clusters
  //************************************************************
  Enable::BECAL_CELL      = Enable::BECAL && true;
  Enable::BECAL_TOWER     = Enable::BECAL_CELL && true;
  Enable::BECAL_CLUSTER   = Enable::BECAL_TOWER && true;
  Enable::BECAL_EVAL      = Enable::BECAL_CLUSTER && false;

  Enable::HCALIN_CELL     = Enable::HCALIN && true;
  Enable::HCALIN_TOWER    = Enable::HCALIN_CELL && true;
  Enable::HCALIN_CLUSTER  = Enable::HCALIN_TOWER && true;
  Enable::HCALIN_EVAL     = Enable::HCALIN_CLUSTER && false;

  Enable::HCALOUT_CELL    = Enable::HCALOUT && true;
  Enable::HCALOUT_TOWER   = Enable::HCALOUT_CELL && true;
  Enable::HCALOUT_CLUSTER = Enable::HCALOUT_TOWER && true;
  Enable::HCALOUT_EVAL    = Enable::HCALOUT_CLUSTER && false;

  Enable::DIRC_RECO       = Enable::DIRC && true;
  Enable::RICH_RECO       = Enable::DIRC && true;
  Enable::mRICH_RECO      = Enable::DIRC && true;

  Enable::FEMC_TOWER      = Enable::FEMC && true;
  Enable::FEMC_CLUSTER    = Enable::FEMC_TOWER && true;
  Enable::FEMC_EVAL       = Enable::FEMC_CLUSTER && false;

  Enable::DRCALO_CELL     = Enable::DRCALO && true;
  Enable::DRCALO_TOWER    = Enable::DRCALO_CELL && true;
  Enable::DRCALO_CLUSTER  = Enable::DRCALO_TOWER && true;
  Enable::DRCALO_EVAL     = Enable::DRCALO_CLUSTER && false;

  Enable::LFHCAL_ABSORBER = false;
  Enable::LFHCAL_CELL     = Enable::LFHCAL && true;
  Enable::LFHCAL_TOWER    = Enable::LFHCAL_CELL && true;
  Enable::LFHCAL_CLUSTER  = Enable::LFHCAL_TOWER && true;
  Enable::LFHCAL_EVAL     = Enable::LFHCAL_CLUSTER && false;

  Enable::EEMCH_TOWER     = Enable::EEMCH && true;
  Enable::EEMCH_CLUSTER   = Enable::EEMCH_TOWER && true;
  Enable::EEMCH_EVAL      = Enable::EEMCH_CLUSTER && false;
  G4TTL::SETTING::optionEEMCH = Enable::EEMCH && true;

  Enable::EHCAL_CELL      = Enable::EHCAL && true;
  Enable::EHCAL_TOWER     = Enable::EHCAL_CELL && true;
  Enable::EHCAL_CLUSTER   = Enable::EHCAL_TOWER && true;
  Enable::EHCAL_EVAL      = Enable::EHCAL_CLUSTER && false;

  Enable::FFR_EVAL        = Enable::HFARFWD_MAGNETS && Enable::HFARFWD_VIRTUAL_DETECTORS && false;

  // Enabling the event evaluator?
  Enable::EVENT_EVAL = true;
  if (detectorSettings.find("HITS") != std::string::npos) {
    Enable::EVENT_EVAL_DO_HITS = true;
    if (detectorSettings.find("HITSABS") != std::string::npos) {
      Enable::EVENT_EVAL_DO_HITS_ABSORBER = true;
    }
    if (detectorSettings.find("HITSC") != std::string::npos) {
      Enable::EVENT_EVAL_DO_HITS_CALO = true;
    }
  }
  Enable::EVENT_EVAL_DO_HEPMC   = Input::PYTHIA6 or Input::PYTHIA8 or Input::SARTRE or Input::HEPMC or Input::READEIC;
  Enable::EVENT_EVAL_DO_EVT_LVL = Input::PYTHIA6 or Input::PYTHIA8 or Input::READEIC;
  // EVENT_EVALUATOR::Verbosity = 1;
  
  //Enable::USER = true;

  //---------------
  // World Settings
  //---------------
  //  G4WORLD::PhysicsList = "FTFP_BERT"; //FTFP_BERT_HP best for calo
  //  G4WORLD::WorldMaterial = "G4_AIR"; // set to G4_GALACTIC for material scans
    //  G4WORLD::WorldMaterial = "G4_Galactic"; // set to G4_GALACTIC for material scans
  // ---------------
  // Magnet Settings
  //---------------
  if (detectorSettings.find("NOFIELD") != std::string::npos) {
    const string magfield = "0.0"; // alternatively to specify a constant magnetic field, give a float number, which will be translated to solenoidal field in T, if string use as fieldmap name (including path)
    G4MAGNET::magfield = magfield;
    G4WORLD::WorldMaterial = "G4_Galactic"; // set to G4_GALACTIC for material scans
  } else {
    // const string magfield = "1.5"; // alternatively to specify a constant magnetic field, give a float number, which will be translated to solenoidal field in T, if string use as fieldmap name (including path)
    //  G4MAGNET::magfield = string(getenv("CALIBRATIONROOT")) + string("/Field/Map/sPHENIX.2d.root");  // default map from the calibration database
    G4MAGNET::magfield_rescale = -1.4 / 1.5;  // make consistent with expected Babar field strength of 1.4T
  }
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

  string outputroot = outdir + "/" + outputFile;
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

  if (Enable::FWDJETS_EVAL) Jet_FwdEval(outputroot);

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
    if (detectorSettings.find("viewer") != std::string::npos){
      gROOT->ProcessLine("PHG4Reco *g4 = QTGui();"); // alternative to DisplayOn
    } else {
      DisplayOn();
      gROOT->ProcessLine("Fun4AllServer *se = Fun4AllServer::instance();");
      gROOT->ProcessLine("PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");");
    }
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
