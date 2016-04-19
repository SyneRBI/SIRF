import math
import pylab
import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../pGadgetron')

from pGadgetron import *
from pGadgets import *

def edge_weight(u):
    shape = u.shape
    nx = shape[0]
    ny = shape[1]
    mx = 2*nx - 1
    my = 2*ny - 1
    w = numpy.ndarray((mx, my), dtype = numpy.float32)
    max_grad = 0
    for iy in range(my):
        jy = iy/2
        for ix in range(mx):
            jx = ix/2
            if ix%2 == 0 and iy%2 == 0:
                continue
            if ix%2 == 1 and iy%2 == 1:
                w[ix, iy] = abs(u[jx + 1, jy + 1] - u[jx, jy])
            elif ix%2 == 1:
                w[ix, iy] = abs(u[jx + 1, jy] - u[jx, jy])
            else:
                w[ix, iy] = abs(u[jx, jy + 1] - u[jx, jy])
            max_grad = max(max_grad, w[ix, iy])
    for jy in range(ny):
        iy = 2*jy
        for jx in range(nx):
            ix = 2*jx
            s = 0
            n = 0
            for ky in range(-1, 2):
                ly = iy + ky
                if ly not in range(my):
                    continue
                for kx in range(-1, 2):
                    lx = ix + kx
                    if lx not in range(mx):
                        continue
                    if kx != ky:
                        s += w[lx, ly]
                        n += 1
            w[ix, iy] = s/n
    ng = mx*my
    nh = ng//10
    hg = numpy.ndarray((nh), dtype = numpy.int16)
    for i in range(nh):
        hg[i] = 0
    for iy in range(my):
        for ix in range(mx):
            j = int((nh - 1.0)*w[ix, iy]/max_grad)
            hg[j] += 1
##    for i in range(10):
##        print(hg[i])
    cutoff = 0
    step = max_grad/nh;
    for i in range(nh):
        if i > 10 and hg[i] < 1:
##            print(i)
            break
        cutoff += step
    print(cutoff, max_grad)
    for iy in range(my):
        for ix in range(mx):
            t = w[ix, iy]/cutoff
            if t > 1:
                w[ix, iy] = 0
            else:
                w[ix, iy] = 1 # (1 - t**4)**2
    for jy in range(ny):
        iy = 2*jy
        for jx in range(nx):
            ix = 2*jx
##            s = 0
##            n = 0
##            for ky in range(-1, 2):
##                ly = iy + ky
##                if ly not in range(my):
##                    continue
##                for kx in range(-1, 2):
##                    lx = ix + kx
##                    if lx not in range(mx):
##                        continue
##                    if kx != 0 or ky != 0:
##                        s += w[lx, ly]
##                        n += 1
##            w[ix, iy] = s/n
            w[ix, iy] = 1
    return w

def smoothen(u, w):
    shape = u.shape
    nx = shape[0]
    ny = shape[1]
    v = numpy.ndarray((nx, ny), dtype = numpy.float32)
    for iy in range(ny):
        jy = 2*iy
        for ix in range(nx):
            jx = 2*ix
            n = 0
            r = 0
            s = 0
            t = 0
            for ky in range(-1, 2):
                ly = iy + ky
                if ly not in range(ny):
                    continue
                for kx in range(-1, 2):
                    lx = ix + kx
                    if lx not in range(nx):
                        continue
                    if kx != 0 and ky != 0:
                        continue
                    if kx != 0 or ky != 0:
                        n += 1
                        r += u[lx, ly]
                        s += u[lx, ly]*w[jx + kx, jy + ky]
                        t += w[jx + kx, jy + ky]
            if t > 0:
                s /= t
                v[ix, iy] = (u[ix, iy] + s)/2
            else:
                v[ix, iy] = (u[ix, iy] + r/n)/2
    return v

try:
    # acquisitions will be read from this HDF file
    #file = str(input('raw data file: '))
    file = 'testdata.h5'
    input_data = MR_Acquisitions(file)

    processed_data = MR_remove_x_oversampling(input_data)

    # perform reconstruction
    recon = MR_BasicReconstruction()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()

    # extract real images from complex
    images = MR_extract_real_images(complex_images)

    # plot obtained images
    ni = images.number()
    for i in range(ni):
        data = images.image_as_array(i)
        pylab.figure(i)
        pylab.imshow(data[0,0,:,:])
        u = data[0,0,:,:]
        w = edge_weight(u)
        pylab.figure(i + ni)
        pylab.imshow(w)
        u = smoothen(u, w)
        pylab.figure(i + 2*ni)
        pylab.imshow(u)
        w = edge_weight(u)
        pylab.figure(i + 3*ni)
        pylab.imshow(w)
        u = smoothen(u, w)
        pylab.figure(i + 4*ni)
        pylab.imshow(u)
        w = edge_weight(u)
        pylab.figure(i + 5*ni)
        pylab.imshow(w)
        u = smoothen(u, w)
        pylab.figure(i + 6*ni)
        pylab.imshow(u)
        w = edge_weight(u)
        pylab.figure(i + 7*ni)
        pylab.imshow(w)
        u = smoothen(u, w)
        pylab.figure(i + 8*ni)
        pylab.imshow(u)
        pylab.show()

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)

