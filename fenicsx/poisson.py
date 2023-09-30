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

def plot_results(domain,V,uh,title,verbose=2):
    tdim = domain.topology.dim
    topology, cell_types, geometry = plot.create_vtk_mesh(domain, tdim)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
    plotter = pyvista.Plotter()
    plotter.add_mesh(grid, show_edges=True)
    plotter.view_xy()
    if verbose>2:
        if not pyvista.OFF_SCREEN:
            plotter.show()
        else:
            figure = plotter.screenshot(f"{title}_mesh.png")
    u_topology, u_cell_types, u_geometry = plot.create_vtk_mesh(V)
    u_grid = pyvista.UnstructuredGrid(u_topology, u_cell_types, u_geometry)
    u_grid.point_data["u"] = uh.x.array.real
    u_grid.set_active_scalars("u")
    u_plotter = pyvista.Plotter()
    u_plotter.add_mesh(u_grid, show_edges=True)
    u_plotter.view_xy()
    if not pyvista.OFF_SCREEN and verbose>2:
        u_plotter.show()
    warped = u_grid.warp_by_scalar()
    plotter2 = pyvista.Plotter()
    plotter2.add_mesh(warped, show_edges=False, show_scalar_bar=True)
    if not pyvista.OFF_SCREEN and verbose>0:
        plotter2.show()
    if dolfinx.cpp.common.has_adios2 and verbose>0:
        fn=f"../paraview/{title}.bp"   
        with io.VTXWriter(domain.comm, fn, [uh]) as vtx:
            vtx.write(0.0)
            print(f"Wrote {fn}")
    if verbose>0:        
        fn=f"../paraview/{title}.xdmf"  
        with io.XDMFFile(domain.comm, fn, "w") as xdmf:
            xdmf.write_mesh(domain)
            xdmf.write_function(uh)
            print(f"Wrote {fn}")

def dirichlet(order=1,verbose=1):
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
    max=uh.x.array.max()
    min=uh.x.array.min()
    print(f"max : {max:.2g}, min : {min:.2g}")   
    if verbose<1:
        return    
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
    if verbose>1:
        plot_results(domain,tdim,V,uh,"dirichlet",
                     verbose=verbose)

import numpy as np
import pyvista

from dolfinx.fem import (Constant, Function, FunctionSpace, 
                         assemble_scalar, 
                         dirichletbc, form, locate_dofs_geometrical)
from dolfinx.fem.petsc import LinearProblem
from dolfinx.mesh import create_unit_square
from mpi4py import MPI
from petsc4py.PETSc import ScalarType
from ufl import (SpatialCoordinate, 
    TestFunction, TrialFunction, dot, ds, dx, grad)

def dirichlet_and_neumann(order=1,verbose=1):
    mesh = create_unit_square(MPI.COMM_WORLD, 10, 10)
    V = FunctionSpace(mesh, ("CG", 1))
    u = TrialFunction(V)
    v = TestFunction(V)
    a = dot(grad(u), grad(v)) * dx
    def u_exact(x):
        return 1 + x[0]**2 + 2*x[1]**2
    def boundary_D(x):
        return np.logical_or(np.isclose(x[0], 0), np.isclose(x[0],1))
    dofs_D = locate_dofs_geometrical(V, boundary_D)
    u_bc = Function(V)
    u_bc.interpolate(u_exact)
    bc = dirichletbc(u_bc, dofs_D)
    x = SpatialCoordinate(mesh)
    g = -4 * x[1]
    f = Constant(mesh, ScalarType(-6))
    L = f * v * dx - g * v * ds
    problem = LinearProblem(a, L, 
        bcs=[bc], petsc_options={
            "ksp_type": "preonly", "pc_type": "lu"})
    uh = problem.solve()
    max=uh.x.array.max()
    min=uh.x.array.min()
    print(f"max : {max:.2g}, min : {min:.2g}")   
    if verbose<1:
        return
    V2 = FunctionSpace(mesh, ("CG", 2))
    uex = Function(V2)
    uex.interpolate(u_exact)
    error_L2 = assemble_scalar(form((uh - uex)**2 * dx))
    error_L2 = np.sqrt(MPI.COMM_WORLD.allreduce(error_L2, op=MPI.SUM))
    u_vertex_values = uh.x.array
    uex_1 = Function(V)
    uex_1.interpolate(uex)
    u_ex_vertex_values = uex_1.x.array
    error_max = np.max(np.abs(u_vertex_values - u_ex_vertex_values))
    error_max = MPI.COMM_WORLD.allreduce(error_max, op=MPI.MAX)
    print(f"Error_L2 : {error_L2:.2e}")
    print(f"Error_max : {error_max:.2e}")
    if verbose>1:
       plot_results(mesh,V,uh,"dirichlet_and_neuman",
            verbose=verbose)
