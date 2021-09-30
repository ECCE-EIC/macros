#ifndef MACRO_EventEvaluator_C
#define MACRO_EventEvaluator_C

#include <fun4all/Fun4AllServer.h>
#include <eiceval/EventEvaluatorEIC.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libeiceval.so)

namespace Enable
{
  // use Enable::EVENT_EVAL = true; in your macro
  bool EVENT_EVAL = false;
}  // namespace Enable

namespace EVENT_EVALUATOR
{
  int Verbosity = 0;
  float EnergyThreshold = 0.05;
}  // namespace EVENT_EVALUATOR

void Event_Eval(const std::string &filename)
{
  Fun4AllServer *se = Fun4AllServer::instance();

  EventEvaluatorEIC *eval = new EventEvaluatorEIC("EVENTEVALUATOR", filename);
  eval->set_reco_tracing_energy_thresholdMC(EVENT_EVALUATOR::EnergyThreshold);
  eval->Verbosity(EVENT_EVALUATOR::Verbosity);

  if (Enable::TRACKING)
  {
    eval->set_do_TRACKS(true);
    //eval->set_do_HITS(true);
    eval->set_do_PROJECTIONS(true);
    if (G4TRACKING::DISPLACED_VERTEX)
      eval->set_do_VERTEX(true);
    if (Enable::DIRC_RECO or Enable::mRICH_RECO or Enable::RICH_RECO)
      eval->set_do_PID_LogLikelihood(true);
  }
  if (Enable::CEMC) eval->set_do_CEMC(true);
  if (Enable::EEMC || Enable::EEMCH) eval->set_do_EEMC(true);
  if (Enable::EEMCH && G4EEMCH::SETTING::USEHYBRID) eval->set_do_EEMCG(true);
  if (Enable::FEMC) eval->set_do_FEMC(true);
  if (Enable::EHCAL) eval->set_do_EHCAL(true);
  if (Enable::HCALIN) eval->set_do_HCALIN(true);
  if (Enable::HCALOUT) eval->set_do_HCALOUT(true);
  if (Enable::FHCAL) eval->set_do_FHCAL(true);
  if (Enable::FHCAL_CLUSTER || Enable::FEMC_CLUSTER || Enable::EEMC_CLUSTER) eval->set_do_CLUSTERS(true);
  if (Enable::DRCALO) eval->set_do_DRCALO(true);
  if (Enable::LFHCAL) eval->set_do_LFHCAL(true);
  if (Enable::BECAL) eval->set_do_BECAL(true);

  eval->set_do_MCPARTICLES(true);
  eval->set_do_HEPMC(Input::PYTHIA6 or Input::PYTHIA8 or Input::SARTRE or Input::HEPMC or Input::READEIC);
  eval->set_do_store_event_level_info(Input::PYTHIA6 or Input::PYTHIA8 or Input::READEIC);
  se->registerSubsystem(eval);

  return;
}

#endif
