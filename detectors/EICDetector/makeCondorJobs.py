import math
import sys, os
from os import environ
from ROOT import TFile, TObjString

simulationsTopDir = '/sphenix/user/cdean/ECCE/DST_files'

nArgs = len(sys.argv)
if nArgs != 5:
    print("Usage: python makeCondorJob <nEventsPerJob> <physics WG> <generator> <collision>")
    sys.exit()

myShell = str(environ['SHELL'])
goodShells = ['/bin/bash', '/bin/tcsh']
if myShell not in goodShells:
    print("Your shell {} was not recognised".format(myShell))
    sys.exit()


nEventsPerJob = int(sys.argv[1])

thisWorkingGroup = sys.argv[2]
ecceWorkingGroup = ['SIDIS', 'HFandJets', 'ExclusiveReactions']
if thisWorkingGroup not in ecceWorkingGroup:
    print("Physics WG: {} was not recognised".format(thisWorkingGroup))
    sys.exit()
else:
    print("Physics WG: {}".format(thisWorkingGroup))


thisGenerator = sys.argv[3]
ecceGenerator = ['pythia6', 'pythia8']
if thisGenerator not in ecceGenerator:
    print("Generator: {} was not recognised".format(thisGenerator))
    sys.exit()
else:
    print("Generator: {}".format(thisGenerator))


thisCollision = sys.argv[4]
ecceCollision = ['ep_10x100', 'ep_18x100']
if thisCollision not in ecceCollision:
    print("Collision: {} was not recognised".format(thisCollision))
    sys.exit()
else:
    print("Collision: {}".format(thisCollision))


def getNumEvtsInFile(theFile):
    inputFile = TFile(theFile)
    nEvents = inputFile.Get("nEvents")
    return nEvents.GetString().Atoi()


def makeCondorJob(PWG):
    print("Creating condor submission files for {} production".format(PWG))
    #Find and open the PWG list of input event files
    inputFileList = "eic-smear_{}_{}_{}.list".format(PWG, thisGenerator, thisCollision)
    infile = open(inputFileList, "r")
    line = infile.readline()
    #Get the current working directory to write submissions and logs to
    myOutputPath = os.getcwd()
    condorDir = "{}/condorJobs".format(myOutputPath)
    os.makedirs("{}/log".format(condorDir), exist_ok=True)
    submitScriptName = "{}/submitJobs.sh".format(condorDir)
    submitScript = open("{}".format(submitScriptName), "w")
    submitScript.write("#!/bin/bash\n")
    #Now make output directory (plus eval folder)
    outputPath = "{}/{}/{}/{}".format(simulationsTopDir, thisWorkingGroup, thisGenerator, thisCollision)
    outputEvalPath = outputPath + "/eval"
    os.makedirs(outputPath, exist_ok=True)
    os.makedirs(outputEvalPath, exist_ok=True)
    #Print input/output info
    print("Input file list: {}".format(inputFileList))
    print("Output directory: {}".format(outputPath))
    #Now loop over all input trees and make a submission script that fits the request
    nJobs = 0
    while line:
       inputFile = line.replace("\n", "")
       nEventsInFile = getNumEvtsInFile(inputFile)
       nJobsFromFile = math.ceil(nEventsInFile/nEventsPerJob)
       for i in range(nJobsFromFile):

            jobNumber = nJobs
            skip = i*nEventsPerJob

            condorOutputInfo = "{0}/log/condor-{1}_{2}_{3}-{4:05d}".format(condorDir, PWG, thisGenerator, thisCollision, jobNumber)

            condorFileName = "condorJob_{0}_{1}_{2}-{3:05d}.job".format(PWG, thisGenerator, thisCollision, jobNumber)
            condorFile = open("{0}/{1}".format(condorDir, condorFileName), "w")
            condorFile.write("Universe        = vanilla\n")
            if myShell == '/bin/bash': condorFile.write("Executable      = {}/run_EIC_production.sh\n".format(myOutputPath))
            if myShell == '/bin/tcsh': condorFile.write("Executable      = {}/run_EIC_production.csh\n".format(myOutputPath))
            outputFile = "DST_{}_{}_{}-{:05d}.root".format(PWG, thisGenerator, thisCollision, jobNumber)
            argument = "{} {} {} {} {}".format(nEventsPerJob, inputFile, outputFile, skip, outputPath)
            condorFile.write("Arguments       = \"{}\"\n".format(argument))
            condorFile.write("Output          = {}.out\n".format(condorOutputInfo))
            condorFile.write("Error           = {}.err\n".format(condorOutputInfo))
            condorFile.write("Log             = {}.log\n".format(condorOutputInfo))
            condorFile.write("Initialdir      = {}\n".format(myOutputPath))
            condorFile.write("PeriodicHold    = (NumJobStarts>=1 && JobStatus == 1)\n")
            condorFile.write("request_memory  = 4GB\n")
            condorFile.write("Priority        = 20\n")
            condorFile.write("job_lease_duration = 3600\n")
            condorFile.write("Queue 1\n")

            submitScript.write("condor_submit {}\n".format(condorFileName))

            nJobs += 1
       line = infile.readline()

    print("Condor submission files have been written to:\n{}/condorJobs".format(myOutputPath))
    print("This setup will submit {} jobs".format(nJobs))
    print("You can submit your jobs with the script:\n{}".format(submitScriptName))


makeCondorJob(thisWorkingGroup)
