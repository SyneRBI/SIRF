class spdhg():
    """Computes a saddle point with a stochastic PDHG.

    This means, a solution (x*, y*), y* = (y*_1, ..., y*_n) such that

    (x*, y*) in arg min_x max_y sum_i=1^n <y_i, A_i> - f*[i](y_i) + g(x)

    where g : X -> IR_infty and f[i] : Y[i] -> IR_infty are convex, l.s.c. and
    proper functionals. For this algorithm, they all may be non-smooth and no
    strong convexity is assumed.

    Parameters
    ----------
    x : primal variable
        This vag.proximalriable is both input and output of the method.
    y : dual variable, optional
        Dual variable is part of a product space. By default equals 0.
    z : variable, optional
        Adjoint of dual variable, z = A^* y. By default equals 0 if y = 0.
    f : functions
        Functionals Y[i] -> IR_infty that all have a convex conjugate with a
        proximal operator, i.e.
        f[i].convex_conj.proximal(sigma[i]) : Y[i] -> Y[i].
    g : function
        Functional X -> IR_infty that has a proximal operator, i.e.
        g.proximal(tau) : X -> X.
    A : functions
        Operators A[i] : X -> Y[i] that possess adjoints: A[i].adjoint
    tau : scalar / vector / matrix
        Step size for primal variable. Note that the proximal operator of g
        has to be well-defined for this input.
    sigma : scalar
        Scalar / vector / matrix used as step size for dual variable. Note that
        the proximal operator related to f (see above) has to be well-defined
        for this input.
    fun_select : function
        Function that selects blocks at every iteration IN -> {1,...,n}. By
        default this is serial uniform sampling, fun_select(k) selects an index
        i \in {1,...,n} with probability 1/n.

    References
    ----------
    [CERS2017] A. Chambolle, M. J. Ehrhardt, P. Richtarik and C.-B. Schoenlieb,
    *Stochastic Primal-Dual Hybrid Gradient Algorithm with Arbitrary Sampling
    and Imaging Applications*. ArXiv: http://arxiv.org/abs/1706.04957 (2017).

    [E+2017] M. J. Ehrhardt, P. J. Markiewicz, P. Richtarik, J. Schott,
    A. Chambolle and C.-B. Schoenlieb, *Faster PET reconstruction with a
    stochastic primal-dual hybrid gradient method*. Wavelets and Sparsity XVII,
    58 (2017) http://doi.org/10.1117/12.2272946.
    """
    
    def __init__(self, x, y, z, f, g, A, tau, sigma, fun_select):
        self.iter = 0
        self.x = x
        
        # I would make y and z optional and default them with 0
        self.y = y
        self.z = z
        
        self.f = f
        self.g = g
        self.A = A
        self.tau = tau
        self.sigma = sigma
        self.fun_select = fun_select

        # Initialize variables
        self.z_relax = z.copy()

    def update(self):
        # select block
        selected = self.fun_select(self.iter)

        # update primal variable
        tmp = (self.x - self.tau * self.z_relax).as_array()
        self.x.fill(self.g.prox(tmp, self.tau))

        # update dual variable and z, z_relax
        self.z_relax = self.z.clone()
        for i in selected:
            # save old yi
            y_old = self.y[i].clone()

            # y[i]= prox(tmp)
            tmp = y_old + self.sigma[i] * self.A[i].direct(self.x)
            self.y[i] = self.f[i].convex_conj.prox(tmp, self.sigma[i])

            # update adjoint of dual variable
            dz = self.A[i].adjoint(self.y[i] - y_old)
            self.z += dz

            # compute extrapolation
            self.z_relax += 2 * dz

        self.iter += 1            
