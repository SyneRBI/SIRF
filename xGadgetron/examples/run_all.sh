#!/bin/csh

#find . -name *.py -print | xargs python

foreach file (*.py)
	echo "==========="
	echo python $file
	python $file
end


