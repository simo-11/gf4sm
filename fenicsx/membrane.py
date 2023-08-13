# https://jsdokken.com/dolfinx-tutorial/chapter1/membrane_code.html
import gmsh
from dolfinx.io import gmshio
from mpi4py import MPI
from dolfinx import fem
import ufl
from petsc4py.PETSc import ScalarType
import numpy as np
import dolfinx.io
gdim = 2

def mesh():
    gmsh.initialize()
    membrane = gmsh.model.occ.addDisk(0, 0, 0, 1, 1)
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(gdim, [membrane], 1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin",0.05)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax",0.05)
    gmsh.model.mesh.generate(gdim)

def main():
    gmsh_model_rank = 0
    mesh_comm = MPI.COMM_WORLD
    domain, cell_markers, facet_markers = gmshio.model_to_mesh(gmsh.model, mesh_comm, gmsh_model_rank, gdim=gdim)   
    V = fem.FunctionSpace(domain, ("CG", 1))
    x = ufl.SpatialCoordinate(domain)
    beta = fem.Constant(domain, ScalarType(12))
    R0 = fem.Constant(domain, ScalarType(0.3))
    p = 4 * ufl.exp(-beta**2 * (x[0]**2 + (x[1] - R0)**2))
    def on_boundary(x):
        return np.isclose(np.sqrt(x[0]**2 + x[1]**2), 1)
    boundary_dofs = fem.locate_dofs_geometrical(V, on_boundary)
    bc = fem.dirichletbc(ScalarType(0), boundary_dofs, V)
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
    L = p * v * ufl.dx
    problem = fem.petsc.LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    uh = problem.solve()
    Q = fem.FunctionSpace(domain, ("CG", 5))
    expr = fem.Expression(p, Q.element.interpolation_points())
    pressure = fem.Function(Q)
    pressure.interpolate(expr)
    pressure.name = "Load"
    uh.name = "Deflection"
    fn="../paraview/membrane.xdmf"
    with dolfinx.io.XDMFFile(MPI.COMM_WORLD, fn, "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(pressure)
        xdmf.write_function(uh)
    print("Wrote {0}".format(fn))    