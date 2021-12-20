#ifndef MACRO_FUN4ALLG4RUNEVALUATORS_C
#define MACRO_FUN4ALLG4RUNEVALUATORS_C

#include <dirent.h>
#include <fstream>
#include <map>
#include <stdlib.h>

#include <GlobalVariables.C>

#include <G4Setup_EICDetector.C>
#include <G4_DSTReader_EICDetector.C>
#include <G4_EventEvaluator.C>
#include <G4_FwdJets.C>
#include <G4_Global.C>
#include <G4_Input.C>
#include <G4_Production.C>
#include <G4_User.C>

#include <fun4all/Fun4AllServer.h>

using namespace std;

R__LOAD_LIBRARY(libfun4all.so)

bool checkForDir(string name)
{
  DIR* dir = opendir(name.c_str());

  return dir == NULL ? 0 : 1;
}

bool checkForFile(string dir, string base, map<string, string> extension)
{
  bool fileExists = false;
  string fileName;
  for (auto const& key : extension)
  {
    fileName = dir + "/" + base + key.second;
    ifstream file(fileName.c_str());
    if (file.good()) fileExists = true;
  }

  return fileExists;
}

int Fun4All_runEvaluators(
    const int nEvents = 1,
    const string &inputFile = "myInputFile.root",
    const string &inputDir = ".",
    const int skip = 0,
    const string &outdir = ".")
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  Input::READHITS = true;
  INPUTREADHITS::filename[0] = inputDir + "/" + inputFile;

  InputInit();

  //-----
  // What to run
  //-----
  
  Enable::DSTREADER = false;
  Enable::USER = false; // Enable this to run custom code from G4_User.C

  Enable::EVENT_EVAL = true;
  Enable::TRACKING_EVAL = true;
  Enable::BECAL_EVAL = true;
  Enable::HCALIN_EVAL = true;
  Enable::HCALOUT_EVAL = true;
  Enable::FEMC_EVAL = true;
  Enable::LFHCAL_EVAL = true;
  Enable::EEMCH_EVAL = true;
  Enable::EHCAL_EVAL = true;
  Enable::FWDJETS_EVAL = false;
  Enable::FFR_EVAL = true;

  map<string, string> evaluatorNames;
  evaluatorNames["event"] = "_g4event_eval.root";
  evaluatorNames["tracking"] = "_g4tracking_eval.root";
  evaluatorNames["becal"] = "_g4becal_eval.root";
  evaluatorNames["hcalin"] = "_g4hcalin_eval.root";
  evaluatorNames["hcalout"] = "_g4hcalout_eval.root";
  evaluatorNames["femc"] = "_g4femc_eval.root";
  evaluatorNames["fhcal"] = "_g4fhcal_eval.root";
  evaluatorNames["eemch"] = "_g4eemch_eval.root";
  evaluatorNames["ehcal"] = "_g4ehcal_eval.root";
  evaluatorNames["farfwd"] = "_g4farfwd_eval.root";

  /*
   * Enable extra features in the evaluators
   * WARNING!!! Ensure your DST has the right nodes
   */
  Enable::TRACKING = true;
  Enable::BECAL = true;
  Enable::HCALIN = true;
  Enable::HCALOUT = true;
  Enable::FEMC = true;
  Enable::LFHCAL = true;
  Enable::EEMCH = true;
  Enable::EHCAL = true;

  TRACKING::ProjectionNames = {"BECAL"
                              ,"CTTL_0"
                              ,"EEMC"
                              ,"EHCAL"
                              ,"ETTL_0"
                              ,"ETTL_1" 
                              ,"FEMC"
                              ,"FTTL_0" 
                              ,"FTTL_1"
                              ,"HCALOUT"
                              ,"LFHCAL"
                              ,"RICH"
                              ,"hpDIRC"
                              ,"mRICH"
                              };

  Enable::DIRC_RECO = true;
  Enable::mRICH_RECO = true;
  Enable::RICH_RECO = true;

  //-----
  // Output file headers and path
  //-----

  //Get base file name
  string baseFile = inputFile;
  string remove_this = ".root";
  size_t pos = baseFile.find(remove_this);
  if (pos != string::npos)
  {
    baseFile.erase(pos, remove_this.length());
  }
 
  string evalDir = outdir;
  string outdirLastChar = outdir.substr(outdir.size() - 1, 1);
  if (outdirLastChar != "/") evalDir += "/";

  unsigned int revisionWidth = 5;
  unsigned int revisionNumber = 0;
  ostringstream evalRevision;
  evalRevision << setfill('0') << setw(revisionWidth) << to_string(revisionNumber);
  evalDir += "eval_" + evalRevision.str();

  while (checkForDir(evalDir))
  {
    bool evalFileStatus = checkForFile(evalDir, baseFile, evaluatorNames);
    if (!evalFileStatus) break;
    evalDir = evalDir.substr(0, evalDir.size() - revisionWidth);
    revisionNumber++;
    evalRevision.str("");
    evalRevision.clear();
    evalRevision << setfill('0') << setw(revisionWidth) << to_string(revisionNumber);
    evalDir += evalRevision.str(); 
  }

  string makeDirectory = "mkdir -p " + evalDir;
  system(makeDirectory.c_str());

  string outputroot = evalDir + "/" + baseFile;

  //-----
  // Reader and User analysis
  //-----

  if (Enable::DSTREADER) G4DSTreader_EICDetector(outputroot + "_DSTReader.root");
  if (Enable::USER) UserAnalysisInit();

  //-----
  // Evaluators
  //-----

  if (Enable::EVENT_EVAL) Event_Eval(outputroot + evaluatorNames.find("event")->second);
  if (Enable::TRACKING_EVAL) Tracking_Eval(outputroot + evaluatorNames.find("tracking")->second);
  if (Enable::BECAL_EVAL) BECAL_Eval(outputroot + evaluatorNames.find("becal")->second);
  if (Enable::HCALIN_EVAL) HCALInner_Eval(outputroot + evaluatorNames.find("hcalin")->second);
  if (Enable::HCALOUT_EVAL) HCALOuter_Eval(outputroot + evaluatorNames.find("hcalout")->second);
  if (Enable::FEMC_EVAL) FEMC_Eval(outputroot + evaluatorNames.find("femc")->second);
  if (Enable::FHCAL_EVAL) FHCAL_Eval(outputroot + evaluatorNames.find("fhcal")->second);
  if (Enable::EEMCH_EVAL) EEMCH_Eval(outputroot + evaluatorNames.find("eemch")->second);
  if (Enable::EHCAL_EVAL) EHCAL_Eval(outputroot + evaluatorNames.find("ehcal")->second);
  if (Enable::FWDJETS_EVAL) Jet_FwdEval(outputroot);
  if (Enable::FFR_EVAL) FFR_Eval(outputroot + evaluatorNames.find("farfwd")->second);

  //--------------
  // Set up Input Managers
  //--------------

  InputManagers();

  //-----
  // Run
  //-----

  se->skip(skip);
  se->run(nEvents);

  //-----
  // Exit
  //-----

  se->End();
  cout << "All done" << endl;
  delete se;

  gSystem->Exit(0);
  return 0;
}

#endif
