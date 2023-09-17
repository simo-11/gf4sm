import gmsh
import numpy as np
from mpi4py import MPI
from dolfinx.io import gmshio
import dolfinx.plot
import ufl
import dolfinx.fem as fem
import dolfinx.mesh as mesh
import pyvista
import petsc4py 

my_globals = {'gmsh_initialized': False}

def generate_mesh_with_crack(Lx = 1.,
    Ly = .5,
    Lcrack = 0.3,
    lc = 0.1,
    dist_min = .1,
    dist_max = .3,
    refinement_ratio = 10,
    verbosity =3):
    """
    Lcrack, length of crack if <=0, no crack will be made
    Ly, height
    lc, caracteristic length of the mesh
    refinement_ratio, how much it is refined at the tip zone
    dist_min, radius of tip zone, not used if there is no crack
    dist_max, radius of the transition zone, not used if there is no crack
    verbosity, >2 will create figures
    """
    gdim = 2
    mesh_comm = MPI.COMM_WORLD
    model_rank = 0
    if not my_globals['gmsh_initialized']:
        gmsh.initialize()
        my_globals['gmsh_initialized']=True
    facet_tags = {"left": 1, "right": 2, 
                  "top": 3, "bottom_no_crack": 5}
    cell_tags = {"all": 20}
    if mesh_comm.rank == model_rank:
        model = gmsh.model()
        model.add("Rectangle")
        model.setCurrent("Rectangle")
        # Create the points
        p1 = model.geo.addPoint(0.0, 0.0, 0, lc)
        if Lcrack>0:
            p2 = model.geo.addPoint(Lcrack, 0.0, 0, lc)
            facet_tags["crack"]=4
        p3 = model.geo.addPoint(Lx, 0, 0, lc)
        p4 = model.geo.addPoint(Lx, Ly, 0, lc)
        p5 = model.geo.addPoint(0, Ly, 0, lc)
        # Create the lines
        if Lcrack>0:
            l1 = model.geo.addLine(p1, p2, tag=facet_tags["crack"])
            l2 = model.geo.addLine(p2, p3, tag=facet_tags["bottom_no_crack"])
        else:
            l2 = model.geo.addLine(p1, p3, tag=facet_tags["bottom_no_crack"])                
        l3 = model.geo.addLine(p3, p4, tag=facet_tags["right"])
        l4 = model.geo.addLine(p4, p5, tag=facet_tags["top"])
        l5 = model.geo.addLine(p5, p1, tag=facet_tags["left"])
        # Create the surface
        if Lcrack>0:
            cloop1 = model.geo.addCurveLoop([l1, l2, l3, l4, l5])
        else:
            cloop1 = model.geo.addCurveLoop([l2, l3, l4, l5])
        surface_1 = model.geo.addPlaneSurface([cloop1])      
        # Define the mesh size and fields for the mesh refinement
        model.mesh.field.add("Distance", 1)
        # SizeMax -                   / ------------------
        #                            /
        # SizeMin -o----------------/
        #          |                |  |
        #        Point        DistMin   DistMax
        model.mesh.field.add("Threshold", 2)
        if Lcrack>0:
            model.mesh.field.setNumbers(1, "NodesList", [p2])
        model.mesh.field.setNumber(2, "IField", 1)
        model.mesh.field.setNumber(2, "LcMin", lc / refinement_ratio)
        model.mesh.field.setNumber(2, "LcMax", lc)
        model.mesh.field.setNumber(2, "DistMin", dist_min)
        model.mesh.field.setNumber(2, "DistMax", dist_max)
        model.mesh.field.setAsBackgroundMesh(2)
        gmsh.option.setNumber("General.Verbosity", 3)
        model.geo.synchronize()
        # Assign mesh and facet tags
        surface_entities = [entity[1] for entity in model.getEntities(2)]
        model.addPhysicalGroup(2, surface_entities, tag=cell_tags["all"])
        model.setPhysicalName(2, 2, "Rectangle surface")
        model.mesh.generate(gdim)
        for (key,value) in facet_tags.items():
            model.addPhysicalGroup(1, [value], tag=value)
            model.setPhysicalName(1, value, key)
    msh, cell_tags, facet_tags = gmshio.model_to_mesh(
        model, mesh_comm, model_rank, gdim=gdim)
    msh.name = "rectangle"
    dx = ufl.Measure("dx",domain=msh)
    top_facets = mesh.locate_entities_boundary(msh, 1, 
            lambda x : np.isclose(x[1], Ly))
    mt = mesh.meshtags(msh, 1, top_facets, 1)
    if verbosity>2:
        cell_tags.name = f"{msh.name}_cells"
        facet_tags.name = f"{msh.name}_facets" 
        fn="../paraview/newfrac-mesh.xdmf"
        with dolfinx.io.XDMFFile(MPI.COMM_WORLD, fn, "w") as file:
            file.write_mesh(msh)
            msh.topology.create_connectivity(1, 2)
            file.write_meshtags(cell_tags, 
                geometry_xpath=f"/Xdmf/Domain/Grid[@Name='{msh.name}']/Geometry")
            file.write_meshtags(facet_tags, 
                geometry_xpath=f"/Xdmf/Domain/Grid[@Name='{msh.name}']/Geometry")
        print("Wrote {0}, Use Wireframe in paraview to view mesh".format(fn))  
    if verbosity>2:
        pyvista.start_xvfb()
        # Extract topology from mesh and create pyvista mesh
        topology, cell_types, x = dolfinx.plot.create_vtk_mesh(msh)
        grid = pyvista.UnstructuredGrid(topology, cell_types, x)
        plotter = pyvista.Plotter()
        plotter.add_mesh(grid, show_edges=True)
        plotter.camera_position = 'xy'
        if not pyvista.OFF_SCREEN:
            plotter.show()
        else:
            plotter.screenshot("mesh.png")
    return msh, mt    

def warp_plot_2d(u,cell_field=None,
        field_name="Field",factor=1.,backend="none",**kwargs):
    #"ipyvtklink", "panel", "ipygany", "static", "pythreejs", "none"
    msh = u.function_space.mesh
    # Create plotter and pyvista grid
    plotter = pyvista.Plotter()
    topology, cell_types, geometry = dolfinx.plot.create_vtk_mesh(msh)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
    # Attach vector values to grid and warp grid by vector
    values = np.zeros((geometry.shape[0], 3), dtype=np.float64)
    values[:, :len(u)] = u.x.array.real.reshape((geometry.shape[0], len(u)))
    grid["u"] = values
    warped_grid = grid.warp_by_vector("u", factor=factor)
    if cell_field is not None:
        warped_grid.cell_data[field_name] = cell_field.vector.array
        warped_grid.set_active_scalars(field_name)
    plotter.add_mesh(warped_grid,**kwargs)
    plotter.camera_position = 'xy'
    return plotter

def solve_elasticity(Lx = 1.,
    Ly = .5,
    Lcrack = 0.3,
    lc = 0.1,
    dist_min = .1,
    dist_max = .3,
    refinement_ratio = 10,
    E=1,
    nu=0.3,
    load=None,
    loadX=None,
    verbosity =1):
    msh, mt = generate_mesh_with_crack(
        Lcrack=Lcrack,
        Ly=Ly,
        lc=lc,  
        refinement_ratio=refinement_ratio,  
        dist_min=dist_min,
        dist_max=dist_max,
        verbosity=verbosity
    )
    element = ufl.VectorElement('Lagrange',msh.ufl_cell(),degree=1,dim=2)
    V = fem.FunctionSpace(msh, element)
    bottom_no_crack_facets = mesh.locate_entities_boundary(
        msh, 1, lambda x : np.logical_and(np.isclose(x[1], 0.0), 
                            x[0] > Lcrack))
    bottom_no_crack_dofs_y = fem.locate_dofs_topological(
        V.sub(1), msh.topology.dim-1, bottom_no_crack_facets)
    right_facets = mesh.locate_entities_boundary(
        msh, 1, lambda x :np.isclose(x[0], Lx))
    right_dofs_x = fem.locate_dofs_topological(
        V.sub(0), msh.topology.dim-1, right_facets)
    bc_bottom = fem.dirichletbc(
        petsc4py.PETSc.ScalarType(0), bottom_no_crack_dofs_y, V.sub(1))
    bc_right = fem.dirichletbc(
        petsc4py.PETSc.ScalarType(0), right_dofs_x, V.sub(0))
    bcs = [bc_bottom, bc_right]
    dx = ufl.Measure("dx",domain=msh)
    top_facets = mesh.locate_entities_boundary(msh, 1, 
        lambda x : np.isclose(x[1], Ly))
    if loadX!=None:
        load_facets = mesh.locate_entities_boundary(msh, 1, 
            lambda x : np.logical_and(np.isclose(x[1], Ly),
                        np.isclose(x[0],loadX,atol=lc)))
        if len(load_facets)<1:
            raise ValueError('No facets found for loadX={0:.3g}'
                             .format(loadX))
        elif len(load_facets)>1:
            print(f"{len(load_facets)} cells will be loaded")
        mt = mesh.meshtags(msh, 1, load_facets, 1)
    ds = ufl.Measure("ds", subdomain_data=mt)
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    # this is for plane-stress
    mu = E / (2.0 * (1.0 + nu))
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
        #b = fem.Constant(msh,petsc4py.PETSc.ScalarType((0, 0)))
        # Surface force on the top
        f = fem.Constant(msh,petsc4py.PETSc.ScalarType(load))
        #return ufl.dot(b, v) * dx + ufl.dot(f, v) * ds(1)
        return ufl.dot(f, v) * ds(1)
    problem = fem.petsc.LinearProblem(a(u,v), L(v), bcs=bcs, 
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    uh = problem.solve()
    uh.name = "displacement"
    energy = fem.assemble_scalar(fem.form(0.5 * a(uh, uh) - L(uh)))
    print(f"The potential energy is {energy:2.3e}")
    if verbosity>0:
        fn="../paraview/newfrac-elasticity.xdmf"
        with dolfinx.io.XDMFFile(MPI.COMM_WORLD, fn, "w") as file:
            file.write_mesh(uh.function_space.mesh)
            file.write_function(uh)
            print("Wrote displacements to {0}".format(fn)) 
    if verbosity>0:
        sigma_iso = 1./3*ufl.tr(sigma(eps(uh)))*ufl.Identity(len(uh))
        sigma_dev =  sigma(eps(uh)) - sigma_iso
        von_Mises = ufl.sqrt(3./2*ufl.inner(sigma_dev, sigma_dev))
        V_von_mises = fem.FunctionSpace(msh, ("DG", 0))
        stress_expr = fem.Expression(von_Mises, 
            V_von_mises.element.interpolation_points())
        vm_stress = fem.Function(V_von_mises)
        vm_stress.interpolate(stress_expr)
        plotter = warp_plot_2d(uh,cell_field=vm_stress,
                            field_name="Von Mises stress",
                            factor=.05,
                            show_edges=True,
                            clim=[0, 0.3])
        if not pyvista.OFF_SCREEN:
            plotter.show()