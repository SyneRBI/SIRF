#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 11:08:29 2018

@author: sirfuser
"""

import numpy

class Function(object):
    def __init__(self):
        pass
    
    def __call__(self, x):
        raise NotImplementedError

    def grad(self, x):
        raise NotImplementedError

    def prox(self, x, tau, out=None):
        raise NotImplementedError
    
    @property
    def convex_conj():
        raise NotImplementedError
        

class ZeroFun(Function):
    
    def __init__(self):
        super(ZeroFun, self).__init__()
    
    def __call__(self, x):
        return 0
    
    def prox(self, x, tau):
        return x

    
class KullbackLeibler(Function):
    def __init__(self, data, background):
        self.data = data
        self.background = background
        self.__offset = None
        
    def __call__(self, x):
        """Return the KL-diveregnce in the point ``x``.

        If any components of ``x`` is non-positive, the value is positive
        infinity.

        Needs one extra array of memory of the size of `prior`.
        """

        # define short variable names
        y = self.data
        r = self.background

        # Compute
        #   sum(x + r - y + y * log(y / (x + r)))
        # = sum(x - y * log(x + r)) + self.offset
        # Assume that
        #   x + r > 0

        # sum the result up
        obj = numpy.sum(x - y * numpy.log(x + r)) + self.offset()

        if numpy.isnan(obj):
            # In this case, some element was less than or equal to zero
            return numpy.inf
        else:
            return obj
    
    @property
    def convex_conj(self):
        """The convex conjugate functional of the KL-functional."""
        return KullbackLeiblerConvexConjugate(self.data, self.background)
    
    def offset(self):
        """The offset which is independent of the unknown."""
    
        if self.__offset is None:
            tmp = self.domain.element()
    
            # define short variable names
            y = self.data
            r = self.background
    
            tmp = self.domain.element(numpy.maximum(y, 1))
            tmp = r - y + y * numpy.log(tmp)
    
            # sum the result up
            self.__offset = numpy.sum(tmp)
    
        return self.__offset

#     def __repr__(self):
#         """to be added???"""
#         """Return ``repr(self)``."""
        # return '{}({!r}, {!r}, {!r})'.format(self.__class__.__name__,
          ##                                    self.domain, self.data,
           #                                   self.background)

    
class KullbackLeiblerConvexConjugate(Function):
    """The convex conjugate of Kullback-Leibler divergence functional.

    Notes
    -----
    The functional :math:`F^*` with prior :math:`g>0` is given by:

    .. math::
        F^*(x)
        =
        \\begin{cases}
            \\sum_{i} \left( -g_i \ln(1 - x_i) \\right)
            & \\text{if } x_i < 1 \\forall i
            \\\\
            +\\infty & \\text{else}
        \\end{cases}

    See Also
    --------
    KullbackLeibler : convex conjugate functional
    """

    def __init__(self, data, background):
        self.data = data
        self.background = background

    def __call__(self, x):
        y = self.data
        r = self.background

        tmp = numpy.sum(- x * r - y * numpy.log(1 - x))

        if numpy.isnan(tmp):
            # In this case, some element was larger than or equal to one
            return numpy.inf
        else:
            return tmp


    def prox(self, x, tau, out=None):
        # Let y = data, r = background, z = x + tau * r
        # Compute 0.5 * (z + 1 - sqrt((z - 1)**2 + 4 * tau * y))
        # Currently it needs 3 extra copies of memory.

        if out is None:
            out = x.clone()

        # define short variable names
        y = self.data.as_array()
        r = self.background.as_array()
        x = x.as_array()

        try:
            taua = tau.as_array()
        except:
            taua = tau

        z = x + tau * r
                    
        out.fill(0.5 * (z + 1 - numpy.sqrt((z - 1) ** 2 + 4 * taua * y)))
        
        return out 
   
    @property
    def convex_conj(self):
        return KullbackLeibler(self.data, self.background)


def mult(x, y):
    try:
        xa = x.as_array()
    except:
        xa = x
        
    out = y.clone()
    out.fill(xa * y.as_array())

    return out

def PowerMethodNonsquare(op, numiters, x0=None):
    # Initialise random
    # Jakob's
    #inputsize = op.size()[1]
    #x0 = ImageContainer(numpy.random.randn(*inputsize)
    # Edo's
    #vg = ImageGeometry(voxel_num_x=inputsize[0],
    #                   voxel_num_y=inputsize[1], 
    #                   voxel_num_z=inputsize[2])
    #
    #x0 = ImageData(geometry = vg, dimension_labels=['vertical','horizontal_y','horizontal_x'])
    #print (x0)
    #x0.fill(numpy.random.randn(*x0.shape))
    
    if x0 is None:
        x0 = op.create_image_data()
    
    #s = numpy.zeros(numiters)
    # Loop
    for it in numpy.arange(numiters):
        x1 = op.adjoint(op.direct(x0))
        x1norm = numpy.sqrt(x1.dot(x1))
        #print ("x0 **********" ,x0)
        #print ("x1 **********" ,x1)
        #s[it] = x1.dot(x0) / x0.dot(x0)
        x0 = (1.0/x1norm)*x1
    return numpy.sqrt(x1norm)
