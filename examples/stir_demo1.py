import numpy
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import sys
sys.path.append('../pSTIR')
import stir
import time

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    # direct all printing to a file
    printer = stir.printerTo('stir_demo1.txt')

    # create OSMAPOSL reconstructor
    recon = stir.OSMAPOSLReconstruction('OSMAPOSL_test_PM_QP2.par')

    # check/redefine some parameters
    f = stir.CylindricFilter(recon.get_inter_iteration_filter())
    f.set_strictly_less_than_radius(True)
    obj = recon.get_objective_function()
    prior = obj.get_prior()
    print('prior penalisation factor:', prior.get_penalisation_factor())
    prior.set_penalisation_factor(0.001)
    print('prior penalisation factor:', prior.get_penalisation_factor())
    am = stir.PoissonLogLh_LinModMean_AcqModData(obj).get_acquisition_model()
    print('tangential_LORs:', am.get_matrix().get_num_tangential_LORs())
    am.get_matrix().set_num_tangential_LORs(2)

    # read an initial estimate for the reconstructed image from a file
    image = stir.Image('my_image0.hv')

    if HAVE_PYLAB:
        # plot the initial image
        data = image.as_array()
        pylab.figure(1)
        pylab.imshow(data[20,:,:])
        pylab.show()

    # set up the reconstructor
    recon.set_up(image)

    # run reconstruction
    start_time = time.time()
    recon.reconstruct(image)
    elapsed_time = time.time() - start_time
    print('elapsed time:', elapsed_time)

    # compare the reconstructed image to the expected image
    expectedImage = stir.Image('expected_image.hv')
    diff = expectedImage.diff_from(image)
    print('difference from expected image:', diff)

    # compare the reconstructed image to the exact image
    exactImage = stir.Image('my_image.hv')

    if HAVE_PYLAB:
        # let the user inspect any z-crossections of the image they want to
        data = image.as_array()
        x_data = exactImage.as_array()
        nz = data.shape[0]
        while True:
            s = str(input('enter z-coordinate: '))
            if len(s) < 1:
                break
            z = int(s)
            if z < 0 or z >= nz:
                break
            pylab.figure(z)
            pylab.imshow(data[z,:,:])
            pylab.figure(100000 + z)
            pylab.imshow(x_data[z,:,:])
            pylab.show()

except stir.error as err:
    # display error information
    print('STIR exception occured:\n', err.value)
