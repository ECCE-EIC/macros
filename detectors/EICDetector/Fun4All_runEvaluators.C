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

  Enable::EVENT_EVAL = false;
  Enable::TRACKING_EVAL = true;
  Enable::CEMC_EVAL = true;
  Enable::HCALIN_EVAL = true;
  Enable::HCALOUT_EVAL = true;
  Enable::FEMC_EVAL = true;
  Enable::FHCAL_EVAL = true;
  Enable::EEMC_EVAL = true;
  Enable::FWDJETS_EVAL = false;

  map<string, string> evaluatorNames;
  evaluatorNames["event"] = "_g4event_eval.root";
  evaluatorNames["tracking"] = "_g4tracking_eval.root";
  evaluatorNames["cemc"] = "_g4cemc_eval.root";
  evaluatorNames["hcalin"] = "_g4hcalin_eval.root";
  evaluatorNames["hcalout"] = "_g4hcalout_eval.root";
  evaluatorNames["femc"] = "_g4femc_eval.root";
  evaluatorNames["fhcal"] = "_g4fhcal_eval.root";
  evaluatorNames["eemc"] = "_g4eemc_eval.root";

  Enable::TRACKING = true;
  Enable::CEMC_CLUSTER = true;
  Enable::EEMC_CLUSTER = true;
  Enable::FEMC_CLUSTER = true;
  Enable::HCALIN_CLUSTER = true;
  Enable::HCALOUT_CLUSTER = true;
  Enable::FHCAL_CLUSTER = true;

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
  if (Enable::CEMC_EVAL) CEMC_Eval(outputroot + evaluatorNames.find("cemc")->second);
  if (Enable::HCALIN_EVAL) HCALInner_Eval(outputroot + evaluatorNames.find("hcalin")->second);
  if (Enable::HCALOUT_EVAL) HCALOuter_Eval(outputroot + evaluatorNames.find("hcalout")->second);
  if (Enable::FEMC_EVAL) FEMC_Eval(outputroot + evaluatorNames.find("femc")->second);
  if (Enable::FHCAL_EVAL) FHCAL_Eval(outputroot + evaluatorNames.find("fhcal")->second);
  if (Enable::EEMC_EVAL) EEMC_Eval(outputroot + evaluatorNames.find("eemc")->second);
  //if (Enable::FWDJETS_EVAL) Jet_FwdEval();

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
