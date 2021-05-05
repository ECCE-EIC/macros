#!/bin/bash

export HOME=/sphenix/u/${LOGNAME}
source /cvmfs/eic.opensciencegrid.org/ecce/gcc-8.3/opt/fun4all/core/bin/ecce_setup.sh -n new.5

export ECCE=/sphenix/user/cdean/ECCE
export MYINSTALL=$ECCE/install
export LD_LIBRARY_PATH=$MYINSTALL/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$MYINSTALL/include:$ROOT_INCLUDE_PATH

source /cvmfs/eic.opensciencegrid.org/ecce/gcc-8.3/opt/fun4all/core/bin/setup_local.sh $MYINSTALL

export ROOT_INCLUDE_PATH=${ECCE}/macros/common:${ROOT_INCLUDE_PATH}

echo 'here comes your environment'
#printenv
echo arg1 \(nEvents\) : $1
echo arg2 \(input file\): $2
echo arg3 \(output file\): $3
echo arg4 \(skip\): $4
echo arg5 \(output dir\): $5
echo running root.exe -q -b Fun4All_G4_EICDetector.C\($1,\"$2\",\"$3\",\"\",$4,\"$5\"\)
root.exe -q -b Fun4All_G4_EICDetector.C\($1,\"$2\",\"$3\",\"\",$4,\"$5\"\)
echo "script done"
