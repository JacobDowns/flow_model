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
        # Fields that need to be loaded
        self.fields = ['B', 'H', 'width', 'beta2']
        # Load model fields
        model_wrapper.load_fields(self.fields)
        

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


        ### Input functions
        ########################################################################

        # Bed elevation
        B = Function(V_cg)
        B.assign(model_wrapper.input_functions['B'])
        # Thickness
        H = Function(V_cg)
        H.assign(model_wrapper.input_functions['H'])
        # Basal traction
        beta2 = Function(V_cg)
        beta2.assign(model_wrapper.input_functions['beta2'])
        # Ice stream width
        width = Function(V_cg)
        width.assign(model_wrapper.input_functions['width'])
        # Ice base
        Bhat = Function(V_cg)
        # Floating indicator
        floating = Function(V_cg)
        # Ocean depth
        depth = Function(V_cg)
        # Effective pressure
        N = Function(V_cg)
        # Terminus (for applying hydrostatic bc)
        terminus = Function(V_cg)
        
        self.H = H
        self.B = B
        self.beta2 = beta2
        self.width = width
        self.S = Bhat + H
        self.Bhat = Bhat
        self.floating = floating
        self.depth = depth
        self.N = N
        self.terminus = terminus

        
        ### Function initialization
        ########################################################################

        # Initialize guesses for unknowns
        self.assigner.assign(U, [self.zero_guess, self.zero_guess])


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


    def update(self):

        ### Update ice base
        ##############################################################
        
        B = self.B.vector().get_local()
        H = self.H.vector().get_local()
        Bhat = np.maximum(B,-(917./1029.)*H)
        self.Bhat.vector().set_local(Bhat)

        ### Update basal traction on floating parts of domain
        ##############################################################

        depth = Bhat - B
        self.depth.vector().set_local(depth)
        floating = np.ones_like(depth)
        floating[depth > 1.] = 1e-12
        self.floating.vector().set_local(floating)


        ### Hydrostatic boundary condition at terminus
        ##############################################################
        
        min_thickness = self.model_wrapper.model_params['mass_params']['min_thickness']
        index = np.where(H == min_thickness)[0].max()
        terminus_vec = np.zeros_like(H)
        terminus_vec[index] = 1.
        self.terminus.vector().set_local(terminus_vec)

        
    # Step the model forward by one time step
    def solve(self):
         
        self.update()

        ### Solve for velocity
        ##############################################################
        
        self.assigner.assign(self.U, [self.zero_guess, self.zero_guess])
        # Solve
        solve(self.R == 0, self.U, self.bcs, solver_parameters = self.snes_params)
        # Update previous solutions
        self.assigner_inv.assign([self.ubar0, self.udef0], self.U)
