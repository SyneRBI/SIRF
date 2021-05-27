This directory contains a few demo files that are intended
to be run interactively. For instance, for Python using spyder.

You should run them cell-by-cell, and see what commands are used
and what plots are made.

In spyder (and other front-ends?) and if you use ipython,
you can choose if plots generate a separate figure or appear
"inline" as follows

     # separate figures
     %matplotlib

     # inline with other output
     %matplotlib inline

The files can be executed in one go, but the plots won't make sense
and depending on your settings, might even not appear (or immediately
disappear).

Input data sits in the `SIRF/data/examples/PET/*` directories. 
When you run the demos, output will be in a sub-directory called 
`working_folder`. You can safely delete this afterwards.

Recommended order
1. `01-display_and_projection.py`
2. `02-image_creation_and_projection.py` (optionally)
3. `03-basic_reconstruction.py`
4. `reconstruction.py`
