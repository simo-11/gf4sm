# https://newfrac.github.io/fenicsx-fracture/
# %% common init
import importlib
import newfrac
import numpy as np
import dolfinx
from mpi4py import MPI
import dolfinx.fem as fem
import dolfinx.mesh as mesh
import ufl
import petsc4py 
# %% reload newfrac
importlib.reload(newfrac)
# %% set dimensions and parameters
Lx = 1.
Ly = 0.5
Lcrack = 0.3
lc =.05
dist_min = .1
dist_max = .3
E = 1. 
nu = 0.3 
mu = E / (2.0 * (1.0 + nu))
# %% Mesh generation and visualization
msh, mt = newfrac.generate_mesh_with_crack(
        Lcrack=Lcrack,
        Ly=Ly,
        lc=lc,  # caracteristic length of the mesh
        refinement_ratio=10,  # how much it is refined at the tip zone
        dist_min=dist_min,  # radius of tip zone
        dist_max=dist_max,  # radius of the transition zone
        verbosity=1
    )
# %% Linear Elasticity
element = ufl.VectorElement('Lagrange',msh.ufl_cell(),degree=1,dim=2)
V = fem.FunctionSpace(msh, element)
def bottom_no_crack(x):
    return np.logical_and(np.isclose(x[1], 0.0), 
                          x[0] > Lcrack)
def top(x):
    return np.isclose(x[1], Ly)
bottom_no_crack_facets = mesh.locate_entities_boundary(
    msh, msh.topology.dim-1, bottom_no_crack)
bottom_no_crack_dofs_y = fem.locate_dofs_topological(
    V.sub(1), msh.topology.dim-1, bottom_no_crack_facets)
top_facets = mesh.locate_entities_boundary(msh, msh.topology.dim-1, top)
top_dofs = fem.locate_dofs_topological(V, msh.topology.dim-1, top_facets)
bc_bottom = fem.dirichletbc(
    petsc4py.PETSc.ScalarType(0), bottom_no_crack_dofs_y, V.sub(1))
bc_top = fem.dirichletbc(
    np.array([0,1/2],dtype=petsc4py.PETSc.ScalarType), top_dofs, V)
bcs = [bc_bottom, bc_top]
dx = ufl.Measure("dx",domain=msh)
top_facets = mesh.locate_entities_boundary(msh, 1, lambda x : np.isclose(x[1], Ly))
ds = ufl.Measure("ds", subdomain_data=mt)
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
# this is for plane-stress
lmbda = 2 * mu * lmbda / ( lmbda + 2 * mu )
def eps(u):
    """Strain"""
    return ufl.sym(ufl.grad(u))
def sigma(eps):
    """Stress"""
    return 2.0 * mu * eps + lmbda * ufl.tr(eps) * ufl.Identity(2)
def a(u,v):
    """The bilinear form of the weak formulation"""
    return ufl.inner(sigma(eps(u)), eps(v)) * dx 
def L(v): 
    """The linear form of the weak formulation"""
    # Volume force
    b = fem.Constant(msh,petsc4py.PETSc.ScalarType((0, 0)))
    # Surface force on the top
    f = fem.Constant(msh,petsc4py.PETSc.ScalarType((0, 0.)))
    return ufl.dot(b, v) * dx + ufl.dot(f, v) * ds(1)
problem = fem.petsc.LinearProblem(a(u,v), L(v), bcs=bcs, 
       petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()
uh.name = "displacement"
energy = fem.assemble_scalar(fem.form(0.5 * a(uh, uh) - L(uh)))
print(f"The potential energy is {energy:2.3e}")
fn="../paraview/newfrac-elasticity.xdmf"
with dolfinx.io.XDMFFile(MPI.COMM_WORLD, fn, "w") as file:
    file.write_mesh(uh.function_space.mesh)
    file.write_function(uh)
    print("Wrote displacements to {0}".format(fn)) 
sigma_iso = 1./3*ufl.tr(sigma(eps(uh)))*ufl.Identity(len(uh))
sigma_dev =  sigma(eps(uh)) - sigma_iso
von_Mises = ufl.sqrt(3./2*ufl.inner(sigma_dev, sigma_dev))
V_von_mises = fem.FunctionSpace(msh, ("DG", 0))
stress_expr = fem.Expression(von_Mises, V_von_mises.element.interpolation_points())
vm_stress = fem.Function(V_von_mises)
vm_stress.interpolate(stress_expr)
todo="""
plotter = dolfinx.plot.warp_plot_2d(uh,cell_field=vm_stress,
                       field_name="Von Mises stress",
                       factor=.5,
                       show_edges=True,
                       clim=[0, 0.3])
import pyvista
if not pyvista.OFF_SCREEN:
    plotter.show()
"""    

# %%
