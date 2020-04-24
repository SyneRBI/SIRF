#!/bin/sh

in1="^.*?This is software[.\s\S]+?\)\."
out1="This is software developed originally for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance Imaging
(http://www.ccppetmr.ac.uk/) and now as part of CCP-SyneRBI."

# Stash unstaged changes
git stash -q --keep-index

readonly FILES=$(git diff --diff-filter=ACMRTUXB --cached --name-only)

for f in $FILES; do
	perl -0777 -i -pe "s#${in1}#${out1}#g" $f
done

# Stage updated files
git add -u

# Re-apply original unstaged changes
git stash pop -q

exit 0