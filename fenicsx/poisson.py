# https://jsdokken.com/dolfinx-tutorial/chapter1/fundamentals_code.html
from mpi4py import MPI
from dolfinx import mesh
from dolfinx.fem import FunctionSpace
from dolfinx import fem
import numpy
import ufl
import dolfinx
from petsc4py.PETSc import ScalarType
import pyvista
from dolfinx import plot
from dolfinx import io
def run(order=1):
    domain = mesh.create_unit_square(
        MPI.COMM_WORLD, 8, 8, mesh.CellType.quadrilateral)
    V = FunctionSpace(domain, ("CG", order))
    uD = fem.Function(V)
    uD.interpolate(lambda x: 1 + x[0]**2 + 2 * x[1]**2)
    # Create facet to cell connectivity required to determine boundary facets
    tdim = domain.topology.dim
    fdim = tdim - 1
    domain.topology.create_connectivity(fdim, tdim)
    boundary_facets = mesh.exterior_facet_indices(domain.topology)
    boundary_dofs = fem.locate_dofs_topological(V, fdim, boundary_facets)
    bc = fem.dirichletbc(uD, boundary_dofs)
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    f = fem.Constant(domain, ScalarType(-6))
    a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
    L = f * v * ufl.dx
    problem = fem.petsc.LinearProblem(
        a, L, bcs=[bc], 
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    uh = problem.solve()
    V2 = fem.FunctionSpace(domain, ("CG", 2))
    uex = fem.Function(V2)
    uex.interpolate(lambda x: 1 + x[0]**2 + 2 * x[1]**2)
    L2_error = fem.form(ufl.inner(uh - uex, uh - uex) * ufl.dx)
    error_local = fem.assemble_scalar(L2_error)
    error_L2 = numpy.sqrt(domain.comm.allreduce(error_local, op=MPI.SUM))
    error_max = numpy.max(numpy.abs(uD.x.array-uh.x.array))
    # Only print the error on one process
    if domain.comm.rank == 0:
        print(f"Error_L2 : {error_L2:.2e}")
        print(f"Error_max : {error_max:.2e}")
    print(pyvista.global_theme.jupyter_backend)
    pyvista.start_xvfb()
    topology, cell_types, geometry = plot.create_vtk_mesh(domain, tdim)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
    plotter = pyvista.Plotter()
    plotter.add_mesh(grid, show_edges=True)
    plotter.view_xy()
    if not pyvista.OFF_SCREEN:
        plotter.show()
    else:
        figure = plotter.screenshot("fundamentals_mesh.png")
    u_topology, u_cell_types, u_geometry = plot.create_vtk_mesh(V)
    u_grid = pyvista.UnstructuredGrid(u_topology, u_cell_types, u_geometry)
    u_grid.point_data["u"] = uh.x.array.real
    u_grid.set_active_scalars("u")
    u_plotter = pyvista.Plotter()
    u_plotter.add_mesh(u_grid, show_edges=True)
    u_plotter.view_xy()
    if not pyvista.OFF_SCREEN:
        u_plotter.show()
    warped = u_grid.warp_by_scalar()
    plotter2 = pyvista.Plotter()
    plotter2.add_mesh(warped, show_edges=False, show_scalar_bar=True)
    if not pyvista.OFF_SCREEN:
        plotter2.show()
    if dolfinx.cpp.common.has_adios2:
        fn="../paraview/poisson.bp"   
        with io.VTXWriter(domain.comm, fn, [uh]) as vtx:
            vtx.write(0.0)
            print(f"Wrote {fn}")
    fn="../paraview/poisson.xdmf"  
    with io.XDMFFile(domain.comm, fn, "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(uh)
        print(f"Wrote {fn}")
       