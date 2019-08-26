from dolfin import *
from support.ice_constants import *
from support.momentum_form_marine import *
from support.mass_form import *
from support.length_form_test import *
from ..support.expressions import *

parameters['form_compiler']['cpp_optimize'] = True
parameters["form_compiler"]["representation"] = "uflacs"
parameters['form_compiler']['quadrature_degree'] = 4
parameters['allow_extrapolation'] = True

class IceModel(object):

    def __init__(self, model_wrapper):

        # Model inputs object
        self.model_wrapper = model_wrapper
        # Mesh
        self.mesh = model_wrapper.mesh
        # Physical constants / parameters
        self.ice_constants = ice_constants
        # Model time
        self.t = float(self.model_wrapper.input_functions['t0'])
        

        #### Function spaces
        ########################################################################

        # Define finite element function spaces.  Here we use CG1 for
        # velocity computations, DG0 (aka finite volume) for mass cons,
        # and "Real" (aka constant) elements for the length

        E_cg = self.model_wrapper.E_cg
        E_dg = self.model_wrapper.E_dg
        E_r =  self.model_wrapper.E_r

        V_cg = self.model_wrapper.V_cg
        V_dg = self.model_wrapper.V_dg
        V_r =  self.model_wrapper.V_r

        self.V_cg = V_cg
        self.V_dg = V_dg
        self.V_r = V_r


        ### Mixed function spaces
        ########################################################################

        # Mixed element
        E_V = MixedElement(E_cg, E_cg, E_cg, E_dg, E_r)
        # Mixed space
        V = FunctionSpace(self.mesh, E_V)
        # For moving data between vector functions and scalar functions
        self.assigner_inv = FunctionAssigner([V_cg, V_cg, V_cg, V_dg, V_r], V)
        self.assigner     = FunctionAssigner(V, [V_cg, V_cg, V_cg, V_dg, V_r])
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
        ubar, udef, H_c, H, L = split(U)
        phibar, phidef, xsi_c, xsi, chi = split(Phi)

        # Values of model variables at previous time step
        un = Function(V_cg)
        u2n = Function(V_cg)
        H0_c = Function(V_cg)
        H0 = Function(V_dg)
        L0 = Function(V_r)

        self.ubar = ubar
        self.udef = udef
        self.H_c = H_c
        self.H = H
        self.L = L
        self.phibar = phibar
        self.phidef = phidef
        self.xsi_c = xsi_c
        self.xsi = xsi
        self.chi = chi
        self.U = U
        self.Phi = Phi
        self.un = un
        self.u2n = u2n
        self.H0_c = H0_c
        self.H0 = H0
        self.L0 = L0
        # Time step
        dt = Constant(1.)
        self.dt = dt
        # 0 function used as an initial velocity guess if velocity solve fails
        self.zero_guess = Function(V_cg)


        ### Input functions
        ########################################################################

        # Bed elevation
        B = Function(V_cg)
        # Basal traction
        beta2 = Function(V_cg)
        # SMB
        adot = Function(V_cg)
        # Ice stream width
        width = Function(V_cg)

        self.B = B
        self.beta2 = beta2
        self.adot = adot
        self.width = width
        # Facet function marking divide and margin boundaries
        self.boundaries = model_wrapper.boundaries
        # Boundary measure
        self.ds1 = dolfin.ds(subdomain_data = self.boundaries)


        ### Function initialization
        ########################################################################

        # Assign initial ice sheet length from data
        L0.vector()[:] = model_wrapper.L_init
        # Initialize initial thickness
        H0.assign(model_wrapper.input_functions['H0'])
        H0_c.assign(model_wrapper.original_cg_functions['H0_c'])
        # Initialize guesses for unknowns
        self.assigner.assign(U, [self.zero_guess, self.zero_guess, H0_c, H0, L0])


        ### Derived expressions + Parameters
        ########################################################################

        # Sea level
        sea_level = Constant(ice_constants['sea_level'])
        # Fraction above flotation where calving begins
        q = ice_constants['q']
        # Water surface, or the greater of bedrock topography or zero
        rho = ice_constants['rho']
        rho_w = ice_constants['rho_w']
        # Ice base
        Bhat = softplus(B,-rho/rho_w*H_c,alpha=0.2) 
        # Water depth
        D = softplus(-(Bhat - sea_level), Constant(0.))
        # Greater of bedrock elevation or water surface
        l = softplus(sea_level, B)
        S = Bhat + H_c
        # Ice surface as DG function
        S_dg = Bhat + H
        # Time derivatives
        dLdt = (L - L0) / dt
        dHdt = (H - H0) / dt
        # Effective pressure
        N = Function(self.V_cg)
        # CG ice thickness at last time step
        self.S0_c = Function(self.V_cg)

        self.S = S
        self.dLdt = dLdt
        self.dHdt = dHdt
        self.dt = dt
        self.N = N
        self.sea_level = sea_level
        self.l = l
        self.Bhat = Bhat
        self.D = D

        ### Temporary variables that store variable values before a step is accepted
        ########################################################################

        self.un_temp = Function(V_cg)
        self.u2n_temp = Function(V_cg)
        self.H0_c_temp = Function(V_cg)
        self.H0_temp = Function(V_dg)
        self.L0_temp = Function(V_r)


        ### Variational forms
        ########################################################################

        # Momentum balance residual
        momentum_form = MomentumForm(self)
        self.momentum_form = momentum_form
        R_momentum = momentum_form.R_momentum

        # Continuous thickness residual
        R_thickness = (H_c - H)*xsi_c*dx

        # Mass balance residual
        mass_form = MassForm(self)
        R_mass = mass_form.R_mass

        # Length residual
        length_form = LengthForm(self)
        self.length_form = length_form
        R_length = length_form.R_length

        # Total residual
        R = R_momentum + R_thickness + R_mass + R_length
        J = derivative(R, U, dU)


        ### Problem setup
        ########################################################################

        # Define variational problem subject to no Dirichlet BCs, but with a
        # thickness bound, plus form compiler parameters for efficiency.
        ffc_options = {"optimize": True}

        # SNES parameters for fixed domain problem
        self.snes_params = {'nonlinear_solver': 'newton',
                      'newton_solver': {
                       'relative_tolerance' : 5e-14,
                       'absolute_tolerance' : 9e-5,
                       'linear_solver': 'gmres',
                       'preconditioner': 'ilu',
                       'maximum_iterations': 35,
                       'report' : False
                       }}

        # Variational problem
        self.problem = NonlinearVariationalProblem(R, U, bcs=[], J=J, form_compiler_parameters = ffc_options)
        # Time step
        self.dt.assign(self.model_wrapper.dt)


    # Step the model forward by one time step
    def step(self, accept = True):
        
        ### Solve
        ####################################################################
        
        try:
            self.assigner.assign(self.U, [self.un, self.u2n, self.H0_c, self.H0, self.L0])
            solver = NonlinearVariationalSolver(self.problem)
            solver.parameters.update(self.snes_params)
            solver.solve()
        except:
            self.assigner.assign(self.U, [self.zero_guess, self.zero_guess, self.H0_c, self.H0, self.L0])
            solver = NonlinearVariationalSolver(self.problem)
            solver.parameters.update(self.snes_params)
            solver.parameters['newton_solver']['error_on_nonconvergence'] = False
            solver.parameters['newton_solver']['relaxation_parameter'] = 0.85
            solver.parameters['newton_solver']['report'] = False
            solver.solve()

        # Update previous solutions
        self.assigner_inv.assign([self.un_temp, self.u2n_temp, self.H0_c_temp, self.H0_temp, self.L0_temp], self.U)
             

        ### Accept the step by updating time
        ####################################################################

        if accept:
            # Update time
            self.t += float(self.dt)
            self.assigner_inv.assign([self.un, self.u2n, self.H0_c, self.H0, self.L0], self.U)
            return float(self.L0)
        else :
            return float(self.L0_temp)


    # Update model inputs
    def update(self, params = {}):

        # Update sea level
        if 'sea_level' in params:
            self.sea_level.assign(params['sea_level'])

        # Update calving parameter
        if 'q' in params:
            self.q.assign(params['q'])

        # Update model bed elevation
        self.B.assign(self.model_wrapper.input_functions['B'])
        # Update model basal traction
        self.beta2.assign(Constant(self.ice_constants['beta2']))
        # Update model surface
        self.S0_c.assign(project(self.Bhat + self.H0_c, self.V_cg))
        # Update model width
        self.width.assign(self.model_wrapper.input_functions['width'])
   
