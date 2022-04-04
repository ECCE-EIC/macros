#!/bin/bash

source /cvmfs/eic.opensciencegrid.org/ecce/gcc-8.3/opt/fun4all/core/bin/ecce_setup.sh -n $6

export ROOT_INCLUDE_PATH=$(pwd)/../../common:$ROOT_INCLUDE_PATH

metaDataFile=${5}/${3}.txt
tmpLogFile=${3}.out

d=`date +%Y/%m/%d`
t=`date +%H:%M`

# Print production details to screen and to metadata file simultaneously
cat << EOF | tee ${metaDataFile}
====== Your production details ======
Production started: ${d} ${t}
Production site: $9
Production Host: ${HOSTNAME}
ECCE build: $6
ECCE macros branch: ${12}
ECCE macros hash: $8
PWG: $7
Generator: ${10} 
Collision type: ${11}
Input file: $2
Output file: $3
Output dir: $5
Number of events: $1
Skip: $4
=====================================
EOF

inputGenerator=
particleType="e-"
if [ "${12}" == "production_AI_Optimization" ]
then
  if [ "${10}" == "particleGun" ]
  then
    inputGenerator="${10}"
    if [ "${11}" == "singlePion" ] 
    then
      particleType="pi-"
    elif [ "${11}" == "singlePionPlus" ] 
    then
      particleType="pi+"
    elif [ "${11}" == "singleMuon" ] 
    then
      particleType="mu-"
    elif [ "${11}" == "singleMuonPlus" ] 
    then
      particleType="mu+"
    elif [ "${11}" == "singleElectron" ] 
    then
      particleType="e-"
    elif [ "${11}" == "singlePositron" ] 
    then
      particleType="e+"
    else
      echo "Your particle type, ${11}, was not recognised"
      exit 1
    fi
  elif [ "${10}" == "pythia8" ]
  then
    inputGenerator="pythia8"
  else
    inputGenerator="EIC-smear"
  fi
fi

# Run Fun4all. Send output to stdout but also capture to temporary local file
if [ "${12}" == "production_AI_Optimization" ]
then
  echo running root.exe -q -b Fun4All_G4_EICDetector.C\($1,\"$2\",\"$3\",\"\",$4,\"$5\"\,\"${inputGenerator}\"\,\"${particleType}\"\)
  root.exe -q -b Fun4All_G4_EICDetector.C\($1,\"$2\",\"$3\",\"\",$4,\"$5\",\"${inputGenerator}\",\"${particleType}\"\) | tee ${tmpLogFile}
else
  echo running root.exe -q -b Fun4All_G4_EICDetector.C\($1,\"$2\",\"$3\",\"\",$4,\"$5\"\)
  root.exe -q -b Fun4All_G4_EICDetector.C\($1,\"$2\",\"$3\",\"\",$4,\"$5\"\) | tee ${tmpLogFile}
fi
rc_dst=$?
echo " rc for dst: $rc_dst"

# Do some basic error handling here: is this failed we need to abort! Continuing might cause broken files on all levels.
if [ ".$rc_dst" != ".0" ] || ! [ -e "$outputPath/$outputFile" ] 
then
  echo " DST production failed. EXIT here, no file copy will be initiated!"
  ls -lhrt $outputPath
  exit $rc_dst
fi

# Scan stdout of Fun4all for random number seeds and add to metadata file
echo production script finished, writing metadata
echo "" >> ${metaDataFile}
echo Seeds: >> ${metaDataFile}
grep 'PHRandomSeed::GetSeed()' ${tmpLogFile} | awk '{print $3}' >> ${metaDataFile}
rm ${tmpLogFile}

echo "" >> ${metaDataFile}
echo md5sum: >> ${metaDataFile}
md5sum ${5}/${3} | awk '{print $1}' >> ${metaDataFile}

echo "DST has been created"
echo "Now producing evaluators"

root.exe -q -b Fun4All_runEvaluators.C\(0,\"$3\",\"$5\",0,\"$5\"\)

rc_eval=$?
echo " rc for eval: $rc_eval"
# Do some more error handling here.
if [ ".$rc_eval" != ".0" ]
then
  echo " EVAL production failed. Delete the potentially broken or incomplete EVAL files."
  echo " --> but keeping the DST and continue to copy."
  evalFileHeader=$(echo $3 | sed 's/.root//')
  rm -rf ${outputPath}/eval_*/${evalFileHeader}*
fi

echo "script done"
