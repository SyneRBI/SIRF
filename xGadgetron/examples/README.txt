coil_sensitivity_maps.py
Python only. Demonstrates methods for calculating coil sensitivity maps.

fully_samp_basic.m
Demo showing reconstruction of fully sampled MR data. 

gen_us_data.m  (requires add_noise.m)
Matlab function to simulate the ISMRMRD HDF5 data from a GRAPPA undersampled acquisition.

grappa_basic.m
Demo that shows basic reconstruction of data acquired with GRAPPA undersampling.

grappa_detail.m
Similar to grappa_basic.m but with more annotations in code and output of g-factors.




Previous Demos (now renamed)
----------------------------
demo1:

Lower level interface demo that illustrates creating and running a chain
of gadgets.

demo2:

Lower level interface demo that illustrates creating and running a chain
of gadgets and gadget sets.

demo3:

Lower level interface demo that illustrates creating and running gadget chains
of 3 types:
- acquisition processing chain
- reconstruction chain
- image processing chain

demo4:

Upper level interface demo that illustrates pre-processing of acquisitions,
reconstructing images and post-processing them.

demo5:

Upper level demo that illustrates the computation of coil sensitivity maps,
applying projection from the image space into acquisition space and back
defined by the aquisition model, and images and acquisitions algebra.

demo6:

Lower level 3-chain GRAPPA reconstruction demo.

demo7:

Upper level 3-steps GRAPPA reconstruction demo.

