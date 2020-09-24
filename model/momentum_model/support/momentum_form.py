import numpy as np
from dolfin import Constant, dx, FacetNormal, conditional, lt
from model.support.expressions import *

class VerticalBasis(object):
    def __init__(self,u,coef,dcoef):
        self.u = u
        self.coef = coef
        self.dcoef = dcoef

    def __call__(self,s):
        return sum([u*c(s) for u,c in zip(self.u,self.coef)])

    def ds(self,s):
        return sum([u*c(s) for u,c in zip(self.u,self.dcoef)])

    def dx(self,s,x):
        return sum([u.dx(x)*c(s) for u,c in zip(self.u,self.coef)])

class VerticalIntegrator(object):
    
    def __init__(self,points,weights):
        self.points = points
        self.weights = weights
        
    def integral_term(self,f,s,w):
        return w*f(s)
    
    def intz(self,f):
        return sum([self.integral_term(f,s,w) for s,w in zip(self.points,self.weights)])

class MomentumForm(object):

    def __init__(self, model):

        ### Ice model fields and constants
        ###############################################################

        physical_constants = model.model_wrapper.model_params['physical_constants']
        momentum_params = model.model_wrapper.model_params['momentum_params']
        
        # Load physical constants
        n = momentum_params['n']
        rho = physical_constants['rho']
        rho_w = physical_constants['rho_w']
        g = physical_constants['g']
        b = Constant(momentum_params['b'])
        eps_reg = Constant(1e-5)

        
        # Thickness
        H = model.H
        # Ice base
        B = model.Bhat
        # Basal traction
        beta2 = model.beta2
        # Surface
        S = model.S
        # Glacier width
        width = model.width
        # Effective pressure
        N = model.N
        # Facet normal vector
        nhat = FacetNormal(model.mesh)
        # Velocity components
        ubar = model.ubar
        udef = model.udef
        # Test functions
        phibar = model.phibar
        phidef = model.phidef
        # Floating indicator
        floating = model.floating


        ### Vertical basis
        ##############################################################
        
        # Vertical test functions
        coef  = [lambda s:1.0, lambda s:1./4.*(5*s**4-1.)] 
        dcoef = [lambda s:0.0, lambda s:5*s**3]     

        # This is a quadrature rule for vertical integration
        points  = np.array([0.0,0.4688,0.8302,1.0])
        weights = np.array([0.4876/2.,0.4317,0.2768,0.0476])

        # Vertical integrator
        vi = VerticalIntegrator(points, weights)
        
        u_ = [ubar, udef]
        phi_ = [phibar, phidef]

        u = VerticalBasis(u_, coef, dcoef)
        phi = VerticalBasis(phi_, coef, dcoef)

        self.u = u
        
        
        ### Momentum residual
        ##############################################################

        def dsdx(s):
            return 1./H*(S.dx(0) - s*H.dx(0))

        def dsdz(s):
            return -1./H

        def eta_v(s):
            return (b/2.)*((u.dx(s,0) + u.ds(s)*dsdx(s))**2 + \
                        0.25*((u.ds(s)*dsdz(s))**2) + eps_reg)**((1.-n)/(2*n))

        def tau_xx(s):
            return (phi.dx(s,0) + phi.ds(s)*dsdx(s))*H*eta_v(s)\
            *(4*(u.dx(s,0) + u.ds(s)*dsdx(s))) + (phi.dx(s,0) \
            + phi.ds(s)*dsdx(s))*H*eta_v(s)*\
            (2*u(s)/width*width.dx(0))     

        def tau_xz(s):
            return dsdz(s)**2*phi.ds(s)*H*eta_v(s)*u.ds(s)

        def tau_d(s):
            return rho*g*H*S.dx(0)*phi(s)

        tau_xy = 2.*H*b/width*((n+2)/(width))**(1./n)\
            *(ubar**2 + Constant(0.01))**((1./n - 1)/2.)*ubar

        tau_b = floating*beta2*N*(abs(u(1)) + Constant(1e-10))

        R_momentum = (-vi.intz(tau_xx) - vi.intz(tau_xz) - phi(1)*tau_b\
                      - vi.intz(tau_d) - phibar*tau_xy)*dx

        R_momentum += Constant(0.5)*(P_0*H - P_w*D)*phibar*dx

        self.R_momentum = R_momentum
