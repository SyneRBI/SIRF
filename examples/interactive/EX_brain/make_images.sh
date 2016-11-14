#! /bin/sh -e
# Script to convert BrainWeb data to useful input for the demo
# Relies on STIR utilities

phantom=phantom_1.0mm_normal_crisp
# decompress BrainWeb data
if [ ! -r ${phantom}.rawb.bz2 ]; then
  bunzip2 -k ${phantom}.rawb.bz2
fi

# convert BrainWeb segmented data to PET input
# Do this using Linux utilities (as we don't want to rely on matlab or so)
# This is of course a bit crazy, but it works for here.

# BrainWEB labels:
#0=Background, 1=CSF, 2=Grey Matter, 3=White Matter, 4=Fat, 5=Muscle/Skin, 6=Skin, 7=Skull, 8=Glial Matter, 9=Connective

# First create the emission image
# We will use "tr" to replace the labels with other integers 
# (we have to use octal notation)
tr '\001\002\003\004\005\006\007\010\011' '\000\002\007\001\001\001\000\001\001' <${phantom}.rawb  > preemission.rawb
# copy Interfile header
sed -e s/${phantom}/preemission/ < ${phantom}.hv > preemission.hv
# zoom z-spacing to what we need for the forward projection
zoom_image emission.hv preemission.hv 211 1 0 0 15 0.148148148 -42.75

# Now create the attenuation image
# We will use "tr" to replace the labels with other integers 
# We have to use unsigned bytes, so we will set "tissue" to 2, and "skull" to 3,
# which is about the correct ratio for 511 keV
tr '\001\002\003\004\005\006\007\010\011' '\000\002\002\002\002\002\003\002\002' <${phantom}.rawb  > preattenuation.rawb
# copy Interfile header
sed -e s/${phantom}/preattenuation/ < ${phantom}.hv > preattenuation.hv

# zoom z-spacing to what we need for the forward projection
zoom_image preattenuation_zoomed.hv preattenuation.hv 211 1 0 0 15 0.148148148 -42.75

# now convert to cm^-1 ("tissue" will become 0.096)
# we have to take the zoom-factor into account for this
stir_math --including-first --times-scalar 0.048 --divide-scalar 6.75 attenuation.hv preattenuation_zoomed.hv


rm pre*v pre*rawb
