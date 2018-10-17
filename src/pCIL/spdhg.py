class spdhg():
    """Computes a saddle point with a stochastic PDHG.

    This means, a solution (x*, y*), y* = (y*_1, ..., y*_n) such that

    (x*, y*) in arg min_x max_y sum_i=1^n <y_i, A_i> - f*[i](y_i) + g(x)

    where g : X -> IR_infty and f[i] : Y[i] -> IR_infty are convex, l.s.c. and
    proper functionals. For this algorithm, they all may be non-smooth and no
    strong convexity is assumed.

    Parameters
    ----------
    f : list of functions
        Functionals Y[i] -> IR_infty that all have a convex conjugate with a
        proximal operator, i.e.
        f[i].convex_conj.prox(sigma[i]) : Y[i] -> Y[i].
    g : function
        Functional X -> IR_infty that has a proximal operator, i.e.
        g.prox(tau) : X -> X.
    A : list of functions
        Operators A[i] : X -> Y[i] that possess adjoints: A[i].adjoint
    x : primal variable, optional
        By default equals 0.
    y : dual variable, optional
        Part of a product space. By default equals 0.
    z : variable, optional
        Adjoint of dual variable, z = A^* y. By default equals 0 if y = 0.
    tau : scalar / vector / matrix, optional
        Step size for primal variable. Note that the proximal operator of g
        has to be well-defined for this input.
    sigma : scalar, optional
        Scalar / vector / matrix used as step size for dual variable. Note that
        the proximal operator related to f (see above) has to be well-defined
        for this input.
    prob : list of scalars, optional
        Probabilities prob[i] that a subset i is selected in each iteration.
        If fun_select is not given, then the sum of all probabilities must
        equal 1.
    A_norms : list of scalars, optional
        Norms of the operators in A. Can be used to determine the step sizes
        tau and sigma and the probabilities prob.
    fun_select : function, optional
        Function that selects blocks at every iteration IN -> {1,...,n}. By
        default this is serial sampling, fun_select(k) selects an index
        i \in {1,...,n} with probability prob[i].

    References
    ----------
    [CERS2018] A. Chambolle, M. J. Ehrhardt, P. Richtarik and C.-B. Schoenlieb,
    *Stochastic Primal-Dual Hybrid Gradient Algorithm with Arbitrary Sampling
    and Imaging Applications*. SIAM Journal on Optimization, 28(4), 2783-2808
    (2018) http://doi.org/10.1007/s10851-010-0251-1 

    [E+2017] M. J. Ehrhardt, P. J. Markiewicz, P. Richtarik, J. Schott,
    A. Chambolle and C.-B. Schoenlieb, *Faster PET reconstruction with a
    stochastic primal-dual hybrid gradient method*. Wavelets and Sparsity XVII,
    58 (2017) http://doi.org/10.1117/12.2272946.
    
    [EMS2018] M. J. Ehrhardt, P. J. Markiewicz and C.-B. Schoenlieb, *Faster 
    PET Reconstruction with Non-Smooth Priors by Randomization and 
    Preconditioning*. (2018) ArXiv: http://arxiv.org/abs/1808.07150 
    """
    
    def __init__(self, f, g, A, x=None, y=None, z=None, tau=None, sigma=None, 
                 prob=None, A_norms=None, fun_select=None):
        # fun_select is optional and by default performs serial sampling
                
        if x is None:
            x = A[0].img_templ.copy()
            x.fill(0)
            
        if y is None:
            if z is not None:
                raise ValueError('y and z have to be defaulted together')
             
            import pSTIR as pet
            
            y = [Ai.sirf2sub(0 * pet.AcquisitionData(Ai.__op__.acq_templ)) 
                 for Ai in A]
            z = 0 * x.copy()
            
        else:
            if z is None:
                raise ValueError('y and z have to be defaulted together')
                            
        if A_norms is not None:
            if tau is not None or sigma is not None or prob is not None:
                raise ValueError('Either A_norms or (tau, sigma, prob) must '
                                 'be given')
                
            tau = 1 / sum(A_norms)
            sigma = [1 / nA for nA in A_norms]
            prob = [nA / sum(A_norms) for nA in A_norms]

            #uniform prob, needs different sigma and tau
            #n = len(A)
            #prob = [1./n] * n
                        
        if fun_select is None:
            if prob is None:
                raise ValueError('prob was not determined')
                
            import numpy
            
            def fun_select(k):
                return [int(numpy.random.choice(len(A), 1, p=prob))]

        self.iter = 0
        self.x = x
        
        self.y = y
        self.z = z
        
        self.f = f
        self.g = g
        self.A = A
        self.tau = tau
        self.sigma = sigma
        self.prob = prob
        self.fun_select = fun_select

        # Initialize variables
        self.z_relax = z.copy()
        self.tmp = self.x.copy()
        
    def update(self):
        # select block
        selected = self.fun_select(self.iter)

        # update primal variable
        #tmp = (self.x - self.tau * self.z_relax).as_array()
        #self.x.fill(self.g.prox(tmp, self.tau))
        self.tmp = - self.tau * self.z_relax
        self.tmp += self.x
        self.x = self.g.prox(self.tmp, self.tau)

        # update dual variable and z, z_relax
        self.z_relax = self.z.copy()
        for i in selected:
            # save old yi
            y_old = self.y[i].copy()

            # y[i]= prox(tmp)
            tmp = y_old + self.sigma[i] * self.A[i].direct(self.x)
            self.y[i] = self.f[i].convex_conj.prox(tmp, self.sigma[i])

            # update adjoint of dual variable
            dz = self.A[i].adjoint(self.y[i] - y_old)
            self.z += dz

            # compute extrapolation
            self.z_relax += (1 + 1 / self.prob[i]) * dz

        self.iter += 1            
