from dolfin import MixedElement, FunctionSpace, FunctionAssigner, Function, \
    TrialFunction, TestFunction, split, Constant, NonlinearVariationalProblem, \
    ds, parameters, derivative, project, NonlinearVariationalSolver, assemble, \
    DirichletBC, SubDomain, near, solve
from model.momentum_model.support.momentum_form import *
from model.support.expressions import *
import matplotlib.pyplot as plt

parameters['form_compiler']['cpp_optimize'] = True
parameters["form_compiler"]["representation"] = "uflacs"
parameters['form_compiler']['quadrature_degree'] = 4
parameters['allow_extrapolation'] = True

class MomentumModel(object):

    def __init__(self, model_wrapper, params = {}):

        # Model inputs object
        self.model_wrapper = model_wrapper
        # Mesh
        self.mesh = model_wrapper.mesh
        

        #### Function spaces
        ########################################################################

        E_cg = self.model_wrapper.E_cg
        E_dg = self.model_wrapper.E_dg

        V_cg = self.model_wrapper.V_cg
        V_dg = self.model_wrapper.V_dg

        self.V_cg = V_cg
        self.V_dg = V_dg


        ### Mixed function spaces
        ########################################################################

        # Mixed element
        E_V = MixedElement(E_cg, E_cg)
        # Mixed space
        V = FunctionSpace(self.mesh, E_V)
        # For moving data between vector functions and scalar functions
        self.assigner_inv = FunctionAssigner([V_cg, V_cg], V)
        self.assigner     = FunctionAssigner(V, [V_cg, V_cg])
        self.V = V


        ### Model unknowns + trial and test functions
        ########################################################################

        # U contains both velocity components, the DG thickness, the CG-projected thickness,
        # and the length
        U = Function(V)
        # Trial Function
        dU = TrialFunction(V)
        # Test Function
        Phi = TestFunction(V)
        # Split vector functions into scalar components
        ubar, udef = split(U)
        phibar, phidef = split(Phi)

        self.ubar = ubar
        self.udef = udef
        self.phibar = phibar
        self.phidef = phidef
        self.U = U
        self.Phi = Phi
        self.zero_guess = Function(V_cg)
        self.ubar0 = Function(V_cg)
        self.udef0 = Function(V_cg)

        # Initialize guesses for unknowns
        self.assigner.assign(U, [self.zero_guess, self.zero_guess])

        
        ### Input functions
        ########################################################################

        # Bed elevation
        self.B = model_wrapper.cg_funcs['B']
        # Thickness
        self.H = model_wrapper.cg_funcs['H']
        # Basal traction
        self.beta2 = model_wrapper.cg_funcs['beta2']
        # Ice stream width
        self.width = model_wrapper.cg_funcs['width']
        # Ice base
        self.Bhat = model_wrapper.cg_funcs['Bhat']
        # Floating indicator
        self.floating = model_wrapper.cg_funcs['floating']
        # Ocean depth
        self.depth = model_wrapper.cg_funcs['depth']
        # Effective pressure
        self.N = model_wrapper.cg_funcs['N']
        # Water pressure
        self.P_w = model_wrapper.cg_funcs['P_w']
        # Overburden pressure
        self.P_0 = model_wrapper.cg_funcs['P_0']
        # Surface
        self.S = self.Bhat + self.H
        # Terminus marker measure
        self.ds_t = model_wrapper.ds_t
        

        ### Variational forms
        ########################################################################

        # Momentum balance residual
        momentum_form = MomentumForm(self)
        R_momentum = momentum_form.R_momentum
        # Total residual
        R = R_momentum
        # Derivative
        J = derivative(R, U, dU)
        self.R = R
        self.J = J


        ### Problem setup
        ########################################################################

        # Define variational problem subject to no Dirichlet BCs, but with a
        # thickness bound, plus form compiler parameters for efficiency.
        ffc_options = {"optimize": True}

        # SNES parameters for fixed domain problem
        self.snes_params = {'nonlinear_solver': 'newton',
                      'newton_solver': {
                       'relative_tolerance' : 5e-14,
                       'absolute_tolerance' : 8e-5,
                       'linear_solver': 'gmres',
                       'preconditioner': 'ilu',
                       'maximum_iterations': 35,
                       'report' : True
                       }}

        def left_boundary(x):
            return x[0] < 1e-10

        def right_boundary(x):
            return abs(x[0] - model_wrapper.x.max()) < 1e-10
        
        bcl = DirichletBC(V.sub(0), Constant(0.), left_boundary)
        bcr = DirichletBC(V.sub(0), Constant(0.), right_boundary)
        self.bcs = [bcl, bcr]
        
        
    # Step the model forward by one time step
    def solve(self):

        ### Solve for velocity
        ##############################################################
        
        self.assigner.assign(self.U, [self.zero_guess, self.zero_guess])
        # Solve
        solve(self.R == 0, self.U, self.bcs, solver_parameters = self.snes_params)
        # Update previous solutions
        self.assigner_inv.assign([self.ubar0, self.udef0], self.U)
        # Update the model wrapper
        self.model_wrapper.cg_funcs['u'].assign(self.ubar0)
        
