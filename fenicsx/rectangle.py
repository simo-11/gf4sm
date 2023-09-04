# Initially based on https://jsdokken.com/dolfinx-tutorial/chapter2/linearelasticity_code.html
# 
import numpy as np
import ufl
from mpi4py import MPI
from petsc4py.PETSc import ScalarType
from dolfinx import mesh, fem, plot, io
from enum import Enum, unique
@unique
class Load(Enum):
    Bending='Bending due to own weight',
    Torsion='Torsion at free end'

class TorsionLoad(fem.Expression):
    def __init__(self, mesh, **kwargs):
        super().__init__(**kwargs)
        self.mesh = mesh
    
    def eval(self, value, x):
        value[0]=x[0]
        value[1]=x[1]
        value[2]=x[2]


def main(l=3,order=1,h=0.1,w=0.1,nx=10,ny=1,nz=1,E=210e9,nu=0.3,
         rho=7800,g=9.8,loadCase=Load.Bending):
    domain = mesh.create_box(MPI.COMM_WORLD, [np.array([0,0,0]), np.array([l, h, w])],
                [nx,ny,nz], cell_type=mesh.CellType.hexahedron)
    V = fem.VectorFunctionSpace(domain, ("CG", order))
    def clamped_boundary(x):
        return np.isclose(x[0], 0)
    fdim = domain.topology.dim - 1
    boundary_facets = mesh.locate_entities_boundary(domain, fdim, clamped_boundary)
    u_D = np.array([0,0,0], dtype=ScalarType)
    bc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, boundary_facets), V)
    T = fem.Constant(domain, ScalarType((0, 0, 0)))
    ds = ufl.Measure("ds", domain=domain)
    lambda_=(E*nu)/((1+nu)*(1-2*nu))
    mu= E/(2*(1+nu))
    def epsilon(u):
        return ufl.sym(ufl.grad(u))
    def sigma(u):
        return lambda_ * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2*mu*epsilon(u)
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    if loadCase==Load.Bending:
        f = fem.Constant(domain, ScalarType((0, -rho*g, 0)))
    elif loadCase==Load.Torsion:
        def y_0(x):
            return np.isclose(x[1],0)
        def y_h(x):
            return np.isclose(x[1],h)
        x1 = mesh.locate_entities(domain,fdim,y_0) #TODO
        f1=fem.Expression('1',x1)
        x2 = mesh.locate_entities(domain,fdim,y_h) 
        f2=fem.Expression('-1',x2)
        f=TorsionLoad()
    else:
        raise Exception('loadCase is not defined')
    a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
    L = ufl.dot(f, v) * ufl.dx + ufl.dot(T, v) * ds
    problem = fem.petsc.LinearProblem(
    a, 
    L, 
    bcs=[bc], 
    petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    uh = problem.solve()
    print("Maximum displacement {0:.3g} for {1}".format(-uh.x.array.min(),loadCase))
    fn="../paraview/rectangle_{0}.xdmf".format(loadCase)
    with io.XDMFFile(domain.comm, fn, "w") as xdmf:
        xdmf.write_mesh(domain)
        uh.name = "Deformation"
        xdmf.write_function(uh)
    print("Wrote displacement results to {0}".format(fn))