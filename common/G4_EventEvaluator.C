#ifndef MACRO_EventEvaluator_C
#define MACRO_EventEvaluator_C

#include <fun4all/Fun4AllServer.h>
#include <g4eval/EventEvaluator.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4eval.so)

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

  EventEvaluator *eval = new EventEvaluator("EVENTEVALUATOR", filename);
  eval->set_reco_tracing_energy_threshold(EVENT_EVALUATOR::EnergyThreshold);
  eval->Verbosity(EVENT_EVALUATOR::Verbosity);

  if (Enable::TRACKING)
  {
    eval->set_do_TRACKS(true);
    //eval->set_do_HITS(true);
    eval->set_do_PROJECTIONS(true);
    if (G4TRACKING::DISPLACED_VERTEX)
      eval->set_do_VERTEX(true);
  }
  if (Enable::CEMC_CLUSTER) eval->set_do_CEMC(true);
  if (Enable::EEMC_CLUSTER) eval->set_do_EEMC(true);
  if (Enable::FEMC_CLUSTER) eval->set_do_FEMC(true);
  if (Enable::HCALIN_CLUSTER) eval->set_do_HCALIN(true);
  if (Enable::HCALOUT_CLUSTER) eval->set_do_HCALOUT(true);
  if (Enable::FHCAL_CLUSTER) eval->set_do_FHCAL(true);
  if (Enable::FHCAL_CLUSTER || Enable::FEMC_CLUSTER || Enable::EEMC_CLUSTER)
    eval->set_do_CLUSTERS(true);

  eval->set_do_MCPARTICLES(true);
  se->registerSubsystem(eval);

  return;
}

#endif
