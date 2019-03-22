from dolfin import *
from support.physical_constants import *
from support.momentum_form import *
from support.mass_form import *
from support.length_form_marine import *

parameters['form_compiler']['cpp_optimize'] = True
parameters["form_compiler"]["representation"] = "uflacs"
parameters['form_compiler']['quadrature_degree'] = 4
parameters['allow_extrapolation'] = True

class ForwardIceModel(object):

    def __init__(self, model_inputs, out_dir, checkpoint_file):

        # Model inputs object
        self.model_inputs = model_inputs
        # Mesh
        self.mesh = model_inputs.mesh
        # Physical constants / parameters
        self.constants = pcs
        # Max domain length
        self.domain_len = float(self.model_inputs.input_functions['domain_len'])
        # Model time
        self.t = float(self.model_inputs.input_functions['t0'])
        

        #### Function spaces
        ########################################################################

        # Define finite element function spaces.  Here we use CG1 for
        # velocity computations, DG0 (aka finite volume) for mass cons,
        # and "Real" (aka constant) elements for the length

        E_cg = self.model_inputs.E_cg
        E_dg = self.model_inputs.E_dg
        E_r =  self.model_inputs.E_r

        V_cg = self.model_inputs.V_cg
        V_dg = self.model_inputs.V_dg
        V_r =  self.model_inputs.V_r

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
        self.boundaries = model_inputs.boundaries


        ### Function initialization
        ########################################################################

        # Assign initial ice sheet length from data
        L0.vector()[:] = model_inputs.L_init
        # Initialize initial thickness
        H0.assign(model_inputs.input_functions['H0'])
        H0_c.assign(model_inputs.input_functions['H0_c'])
        # Initialize guesses for unknowns
        self.assigner.assign(U, [self.zero_guess, self.zero_guess, H0_c, H0, L0])


        ### Derived expressions
        ########################################################################

        # Ice surface
        S = B + H_c
        # Ice surface as DG function
        S_dg = B + H
        # Time derivatives
        dLdt = (L - L0) / dt
        dHdt = (H - H0) / dt
        # Overburden pressure
        P_0 = Constant(self.constants['rho']*self.constants['g'])*H_c
        # Overburden fraction
        P_frac = Constant(self.constants['P_frac'])
        # Water pressure
        P_w = P_frac*P_0
        # Effective pressure
        N = P_0 - P_w
        # Sea level
        self.sea_level = Constant(self.constants['sea_level'])
        # Minimum ice thickness
        self.min_thickness = Constant(15.0)
        # CG ice thickness at last time step
        self.S0_c = Function(self.V_cg)
        # SMB expression
        self.adot_prime = model_inputs.get_adot_exp(self.S0_c)
        # SMB as a function
        self.adot_prime_func = Function(self.V_cg)
        # Precip. (snowfall) as a function
        self.precip_func = Function(self.V_cg)

        self.S = S
        self.dLdt = dLdt
        self.dHdt = dHdt
        self.dt = dt
        self.P_frac = P_frac
        self.P_0 = P_0
        self.P_w = P_w
        self.N = N
        

        ### Temporary variables that store variable values before a step is accepted
        ########################################################################

        self.un_temp = Function(V_cg)
        self.u2n_temp = Function(V_cg)
        self.H0_c_temp = Function(V_cg)
        self.H0_temp = Function(V_dg)
        self.L0_temp = Function(V_r)


        ### Initialize inputs
        ########################################################################

        self.update_inputs(model_inputs.L_init, 0.)
        self.S0_c.assign(self.B + self.H0_c)
        self.update_inputs(model_inputs.L_init, 0.)


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

        # Variable length problem
        self.problem = NonlinearVariationalProblem(R, U, bcs=[], J=J, form_compiler_parameters = ffc_options)
        

        ### Setup the iterator for replaying a run
        ########################################################################

        # Get the time step from input file
        self.dt.assign(self.model_inputs.dt)
        # Iteration count
        self.i = 0


    # Assign input functions from model_inputs
    def update_inputs(self, L, precip_param = 0.0):
        self.S0_c.assign(self.B + self.H0_c)
        self.model_inputs.update_inputs(L, self.t, precip_param)
        self.B.assign(self.model_inputs.input_functions['B'])
        self.beta2.assign(self.model_inputs.input_functions['beta2'])
        self.adot_prime_func.assign(project(self.adot_prime, self.V_cg))
        self.width.assign(self.model_inputs.input_functions['width'])
        self.precip_func.assign(self.model_inputs.precip_func)


    # Step the model forward
    def step(self, precip_param = 0.0, accept = False, output = False):
        
        # Update length variable inputs
        self.update_inputs(float(self.L0), precip_param)

        
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
            self.i += 1
            self.assigner_inv.assign([self.un, self.u2n, self.H0_c, self.H0, self.L0], self.U)

            if output:
                print("real step: ", self.t, self.H0_c.vector().max(), self.H0_c.vector().get_local()[0], float(self.L0))
                
            return float(self.L0)
        else :
            return float(self.L0_temp)

        return float(self.L0)


    # Write the ice sheet model state to a file
    def write_state(self, output_file_name):
        
        output_file = HDF5File(mpi_comm_world(), output_file_name + '.h5', 'w')
        output_file.write(self.H0, "H0")
        output_file.write(self.H0_c, "H0_c")
        output_file.write(self.L0, "L0")
        output_file.write(interpolate(Constant(self.t), self.V_r), "t0")
        output_file.flush()
        output_file.close()
