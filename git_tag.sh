#!/bin/sh

IP="IP6"
#IP="IP8"

git checkout Diff_Tagg_$IP

#grep -n '^  Enable::HFARFWD_ION_ENERGY' detectors/EICDetector/Fun4All_G4_EICDetector.C

linenumber=$(grep -n '^  Enable::HFARFWD_ION_ENERGY' detectors/EICDetector/Fun4All_G4_EICDetector.C | cut -d : -f 1)

echo $linenumber


for VAR in 41 62 82 100 135 165 200 220 249 275
do
	echo "Setting up for $IP $VAR GeV beamline"
	sed -i.bak -e "$linenumber"' d' detectors/EICDetector/Fun4All_G4_EICDetector.C
	sed -i.bak -e "$linenumber"' i\' -e "  Enable::HFARFWD_ION_ENERGY = $VAR;" detectors/EICDetector/Fun4All_G4_EICDetector.C
	git add -u .
	git commit -m "$VAR GeV"
	git tag Diff_Tagg_"$IP"_"$VAR"GeV_r6 -f
#	exit
done



