import numpy
import pylab
import stir
import time

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    # direct all diagnostic printing to a file
    printer = stir.printerTo('stir_demo.txt')

    # create OSMAPOSL reconstructor
    recon = stir.OSMAPOSLReconstruction('OSMAPOSL_test_PM_QP.par')

    # check/redefine some parameters
    f = stir.CylindricFilter(recon.get_inter_iteration_filter())
    f.set_strictly_less_than_radius(True)
    obj = recon.get_objective_function()
    prior = obj.get_prior()
    print('prior penalisation factor:', prior.get_penalisation_factor())
    prior.set_penalisation_factor(0.5)
    am = stir.PLL_LMM_AMD(obj).get_acquisition_model()
    print('tangential_LORs:', am.get_matrix().get_num_tangential_LORs())
    am.get_matrix().set_num_tangential_LORs(2)

    # read an initial estimate for the reconstructed image from a file
    image = stir.Image('my_uniform_image_circular.hv')

    # set up the reconstructor
    recon.set_up(image)

    # run reconstruction
    start_time = time.time()
    recon.reconstruct(image)
    elapsed_time = time.time() - start_time
    print('elapsed time:', elapsed_time)

    # compare the reconstructed image to the expected image
    expectedImage = stir.Image('test_image_PM_QP_6.hv')
    diff = expectedImage.diff_from(image)
    print('difference from expected image:', diff)

    # let the user inspect any z-crossections of the image they want to
    data = image.density()
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
        pylab.show()

except stir.error as err:
    # display error information
    print('STIR exception occured:\n', err.value)
