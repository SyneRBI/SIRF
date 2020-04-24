#!/bin/sh

# C
in1="\nThis is software[.\s\S]+?\)\."
out1="\nThis is software developed originally for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance Imaging
(http://www.ccppetmr.ac.uk/) and now as part of CCP-SyneRBI."

# Matlab
in2="% This is software[.\s\S]+?\)\."
out2="% This is software developed originally for the Collaborative Computational
% Project in Positron Emission Tomography and Magnetic Resonance Imaging
% (http://www.ccppetmr.ac.uk/) and now as part of CCP-SyneRBI."

# Python
in3="# This is software[.\s\S]+?\)\."
out3="# This is software developed originally for the Collaborative Computational
# Project in Positron Emission Tomography and Magnetic Resonance Imaging
# (http://www.ccppetmr.ac.uk/) and now as part of CCP-SyneRBI."

# Stash unstaged changes
git stash -q --keep-index

readonly FILES=$(git diff --diff-filter=ACMRTUXB --cached --name-only)

ThereWereErrors=0

for f in $FILES; do
	perl -0777 -i -pe "s&${in1}&${out1}&g" $f &&
	perl -0777 -i -pe "s&${in2}&${out2}&g" $f &&
	perl -0777 -i -pe "s&${in3}&${out3}&g" $f
	if [ $? -ne 0 ]; then ThereWereErrors=1; fi
done

# Stage updated files
git add -u

# Re-apply original unstaged changes
git stash pop -q

exit $ThereWereErrors