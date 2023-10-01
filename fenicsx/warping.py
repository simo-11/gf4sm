# Initially based on https://jsdokken.com/dolfinx-tutorial/chapter1/fundamentals_code.html
# target is to solve warping function
# 
import numpy as np
import ufl
from mpi4py import MPI
from petsc4py.PETSc import ScalarType
from dolfinx import mesh, fem, plot, io
def main(order=1,h=0.1,w=0.1,nx=10,ny=10):
    domain = mesh.create_rectangle(MPI.COMM_WORLD, 
        [np.array([0.0, 0.0]), np.array([w, h])],
        [nx,ny], 
        cell_type=mesh.CellType.quadrilateral)
    V = fem.VectorFunctionSpace(domain, ("CG", order))
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    F = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx 
    x = ufl.SpatialCoordinate(domain)
    n = ufl.FacetNormal(domain)
    atol_for_corner=max(h/ny,w/nx)
    boundaries = [(1, lambda x: np.isclose(x[0], 0)),
              (2, lambda x: np.isclose(x[0], w)),
              (3, lambda x: np.isclose(x[1], 0)),
              (4, lambda x: np.isclose(x[1], h)),
              (5, lambda x: np.isclose(x[0]+x[1], 0, atol=atol_for_corner))]
    facet_indices, facet_markers = [], []
    fdim = domain.topology.dim - 1
    for (marker, locator) in boundaries:
        facets = mesh.locate_entities(domain, fdim, locator)
        if len(facets)<1:
            raise ValueError(f"No facets found for marker {marker}"
                             ", redefine boundaries")
        facet_indices.append(facets)
        facet_markers.append(np.full_like(facets, marker))
    facet_indices = np.hstack(facet_indices).astype(np.int32)
    facet_markers = np.hstack(facet_markers).astype(np.int32)
    sorted_facets = np.argsort(facet_indices)
    facet_tag = mesh.meshtags(domain, fdim, 
                              facet_indices[sorted_facets], 
                              facet_markers[sorted_facets])
    domain.topology.create_connectivity(domain.topology.dim-1, domain.topology.dim)
    fn="../paraview/warping_facet_tags.xdmf"
    with io.XDMFFile(domain.comm, fn, "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_meshtags(facet_tag)
        print(f"Wrote tags to {fn}.")
    ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tag)
    class BoundaryCondition():
        def __init__(self, type, marker, values):
            self._type = type
            if type == "Dirichlet":
                u_D = fem.Constant(domain,ScalarType(values))
                facets = facet_tag.find(marker)
                dofs = fem.locate_dofs_topological(V, fdim, facets)
                self._bc = fem.dirichletbc(u_D, dofs)
            elif type == "Neumann":
                self._bc = ufl.inner(values, v) * ds(marker)
            elif type == "Robin":
                self._bc = values[0] * ufl.inner(u-values[1], v)* ds(marker)
            else:
                raise TypeError("Unknown boundary condition: {0:s}".format(type))
        @property
        def bc(self):
            return self._bc

        @property
        def type(self):
            return self._type    
    boundary_conditions = [BoundaryCondition("Neumann", 1, n),
        BoundaryCondition("Neumann", 2, n),
        BoundaryCondition("Neumann", 3, n),
        BoundaryCondition("Neumann", 4, n),
        BoundaryCondition("Dirichlet", 5, 0.)] 
    bcs = []
    for condition in boundary_conditions:
        if condition.type == "Dirichlet":
            bcs.append(condition.bc)
        else:
            F += condition.bc
    a = ufl.lhs(F)
    L = ufl.rhs(F)
    problem = fem.petsc.LinearProblem(a, L, bcs=bcs, 
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    uh = problem.solve()
    fn="../paraview/warping.xdmf"
    with io.XDMFFile(domain.comm, fn, "w") as xdmf:
        xdmf.write_mesh(domain)
        uh.name = "Warping"
        xdmf.write_function(uh)
    print("Wrote warping to {0}".format(fn))