"""
Implementation of the theory from Faltus Horák Doškář Rokoš 2024:
    Third Medium Finite Element Contact Formulation 
    for Pneumatically Actuated Systems
"""

import fenics as fe
import os
import time
from datetime import datetime

# create the mesh
ntol = 1e-3
C_tickness = 0.1
length = 1 # [m]
height = 0.5
mesh = fe.RectangleMesh.create(
    [fe.Point(0, 0),fe.Point(length + C_tickness/2, height)],
    [21,10],
    fe.CellType.Type.quadrilateral
)

# define the whole domain
materials = fe.MeshFunction("size_t", mesh, mesh.topology().dim())

# initialize marker for the whole domain
materials.set_all(1)

# Define third medium material subdomain
third_medium = fe.CompiledSubDomain(
    "x[0]>tt-tol && x[1]>tt-tol && x[1]<alt-tt+tol", 
    tt = C_tickness, 
    lar = length,
    alt = height,
    tol = ntol
)

# mark third medium material subdomain
third_medium.mark(materials, 2)

# Define third medium material subdomain
additional_third_medium = fe.CompiledSubDomain(
    "x[0]>lar-additional",
    additional = C_tickness/2,
    lar = length
)

# mark third medium material subdomain
additional_third_medium.mark(materials, 2)

# define integration measure based on materials domains
dx = fe.Measure('dx', domain=mesh, subdomain_data=materials)

# define zone for fix boundary condition
fix = fe.CompiledSubDomain("near(x[0], 0.0) && on_boundary")

# Define the whole boundary
boundaries = fe.MeshFunction("size_t", mesh, mesh.topology().dim()-1)

# initialize marker for the whole boundary
boundaries.set_all(0)

# mark fix zone boundary
fix.mark(boundaries, 1)

# define integration measure for fix zone boundary
ds_fix = fe.Measure("ds", subdomain_id=1, subdomain_data=boundaries)

# define zone for displacement input
load_point = fe.CompiledSubDomain(
    "near(x[0], lar-tt) && near(x[1], alt)",
    tt = C_tickness, 
    lar=length,
    alt=height
)

# define function space
SV = fe.VectorFunctionSpace(mesh, "Lagrange", 2)

# create expression for load increment
load = fe.Expression("incr", incr = 0.0, degree = 1)

# Define Dirichlet boundary conditions
bc1 = fe.DirichletBC(SV, fe.Constant((0.0, 0.0)), fix)
bc2 = fe.DirichletBC(SV.sub(1), load, load_point, method="pointwise")
bcs = [bc1, bc2]

# Define functions
utrial = fe.TrialFunction(SV)            # Incremental displacement
utest = fe.TestFunction(SV)             # Test function
usol = fe.Function(SV, name="Displacement")
uold = fe.Function(SV)

# vectorial space geometric dimension
dim = usol.geometric_dimension()

# Identity tensor
II = fe.Identity(dim)

# Deformation gradient
deformation = II + fe.grad(usol)

# Right Cauchy-Green tensor
CauchyGreen = deformation.T*deformation

# Invariants of deformation tensors
inv1_CauchyGreen = fe.tr(CauchyGreen)
inv3_deformation = fe.det(deformation)

# Stored strain energy density (compressible neo-Hookean model)
bulk_K = 2000e6
bulk_G = 10e6
psi_bulk = bulk_K/2*fe.ln(inv3_deformation)**2 + bulk_G/2*( inv3_deformation**(-2/3)*inv1_CauchyGreen - 3 )

# Third medium  material properties
cr = 1
gamma = 1
kr = 2e3

# Third medium contact fictions strain energy
psi_contact = fe.ln(inv3_deformation)**2 + ( inv3_deformation**(-2/3)*inv1_CauchyGreen - 3 )

# measure of rotation associated with material spin
ST = fe.TensorFunctionSpace(mesh, "Lagrange", 1, shape=(dim, dim, dim))
rotation_measure_old = fe.Function(ST)

def rotation_measure(unew, uold, rotation_measure_old):
    Deltau = unew - uold
    ff = fe.grad(Deltau)
    return 1/2*( fe.grad(ff) - fe.grad(ff.T) ) + rotation_measure_old

# Third medium regularization fictions strain energy
psi_regularization = 1/2*cr*( \
    fe.inner( rotation_measure(usol, uold, rotation_measure_old), rotation_measure(usol, uold, rotation_measure_old)) + \
    fe.inner( fe.grad(inv3_deformation), fe.grad(inv3_deformation))  
)

# assemble the third medium fictious strain energy
psi_third_medium = gamma*psi_contact + gamma*kr*psi_regularization

# Total potential energy
Pi = psi_bulk*dx(1) + psi_third_medium*dx(2)

# Compute first variation of Pi (directional derivative about u in the direction of utest)
Pi_variation = fe.derivative(Pi, usol, utest)

# Compute Jacobian of Pi_variation
Pi_Jacobian = fe.derivative(Pi_variation, usol, utrial)

# iterative procedure parameters
simulation_time = 0
dt = 0.001
tmax = 1
utot = -height
solution_number = 0
non_convergence_counter = 0

# set output folder
now = datetime.now()
date_string = now.strftime("%Y%m%dT%H%M")[2:]
outfolder = date_string + "/"

# setup field output file
ffile = fe.XDMFFile(outfolder + "third_medium_contact.xdmf")
ffile.parameters["flush_output"]=True
ffile.parameters["functions_share_mesh"]=True

# output material domains as field output
fe.XDMFFile( outfolder + "materials.xdmf" ).write(materials)

# copy the simulation file to the results folder for forensics
os.system("cp " + __file__ + " ./" + outfolder +
        __file__.split("/")[-1].replace(".py", "") +
        "_" + date_string + ".py" )

# redirect history output to terminal and file
def fprint(string):
    print(string)
    with open(outfolder + "README.md", "a") as file:
        file.writelines( "\n" + string)

# print header for history output
info_string =  f"{'cpu mins':^10}\t"
info_string += f"{'% complete':^10}\t"
info_string += f"{'dt':^10}\t"
info_string += f"{'sim sec':^10}\t"
info_string += f"{'incr disp':^10}\t"
info_string += f"{'E_tot':^10}\t"
info_string += f"{'E_bulk':^10}\t"
info_string += f"{'E_contact':^10}\t"
info_string += f"{'E_regul':^10}\t"
info_string += f"{'Reaction':^10}\t"
info_string += f"{'num conv':^10}\t"
info_string += f"{'solution':^10}\t"
fprint(info_string)

# set form compiler options
fe.parameters["form_compiler"]["cpp_optimize"] = True
fe.parameters['form_compiler']['quadrature_degree'] = 2
fe.parameters['form_compiler']['representation'] = 'uflacs'

# set solver compiler parameters
my_compiler_params = {
    "optimize" : True,
    "eliminate_zeros" : True,
    "precompute_basis_const" : True,
    "precompute_ip_const" : True
}

# set solver parameters
my_sol_params = {
    "newton_solver" : {
        "absolute_tolerance" : 1e-6,
        "relative_tolerance" : 1e-6,
        "krylov_solver" : {
            "nonzero_initial_guess" : True,
        },
        "maximum_iterations" : 30
    }
}

# small strain stress function for computation of reaction force
def sigma(u):
    return (bulk_K - 2*bulk_G/3)*fe.div(u)*fe.Identity(dim) + 2*bulk_G*fe.sym(fe.grad(u))

# keep track of cpu time 
tic = time.time()

# iterative solution procedure
while simulation_time+ntol < tmax:
    
    # update displacement input
    increment_disp = (simulation_time+dt)/tmax*utot
    load.user_parameters["incr"] = increment_disp
    
    # set previous solution as initial guess
    usol.vector().set_local( uold.vector().get_local() )
    usol.vector().apply("insert")

    # Solve variational problem
    convergence = True
    try:
        fe.solve(Pi_variation == 0, usol, bcs, J=Pi_Jacobian,
            form_compiler_parameters = my_compiler_params,
            solver_parameters = my_sol_params
        )
    except:
        convergence = False
        non_convergence_counter += 1

    # compute history variables
    E_tot = fe.assemble( Pi )
    E_bulk = fe.assemble( psi_bulk*dx(1) )
    E_contact = fe.assemble( psi_contact*dx(2) )
    E_regularization = fe.assemble( psi_regularization*dx(2) )
    reaction = fe.assemble( sigma(usol)[1,1]*ds_fix )
    
    # printout history output
    info_string =  f"{int((time.time() - tic)/60):^10.3g}\t"
    info_string += f"{simulation_time/tmax*100:^10.3}\t"
    info_string += f"{dt:^10.3g}\t"
    info_string += f"{simulation_time+dt:^10.3g}\t"
    info_string += f"{increment_disp:^10.3g}\t"
    info_string += f"{E_tot:^10.3g}\t"
    info_string += f"{E_bulk:^10.3g}\t"
    info_string += f"{E_contact:^10.3g}\t"
    info_string += f"{E_regularization:^10.3g}\t"
    info_string += f"{reaction:^10.3g}\t"
        
    if convergence:
        # reset non-converged iterations counter
        non_convergence_counter = 0
        
        # bring simulation time forward
        simulation_time += dt
        
        # increment time step
        dt = dt*1.5
        
        # printout solution identifier
        solution_number += 1
        info_string += f"{'-':^10}\t"
        info_string += f"{solution_number:^10g}\t"
        
        # Save field output
        ffile.write(usol, solution_number)
        
        # update incremental rotation measure 
        rotation_measure_old.assign( fe.project( rotation_measure(usol, uold, rotation_measure_old), ST))
        uold.assign( usol )
    else:
        # print number of non-converged iterations
        info_string += f"{non_convergence_counter:^10}\t"
        info_string += f"{'-':^10}\t"
        
        # cut time step
        dt = dt/2
        
    fprint(info_string)
    
    if non_convergence_counter >= 15:
        fprint(f"Exiting simulation, convergence judged unlikely")
        break
