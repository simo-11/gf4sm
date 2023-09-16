import gmsh
import numpy as np
from mpi4py import MPI
from dolfinx.io import gmshio
import dolfinx.plot
import ufl
import dolfinx.fem as fem
import dolfinx.mesh as mesh
import pyvista

my_globals = {'gmsh_initialized': False}

def generate_mesh_with_crack(Lx = 1.,
    Ly = .5,
    Lcrack = 0.3,
    lc = 0.1,
    dist_min = .1,
    dist_max = .3,
    refinement_ratio = 10,
    verbosity =1):
    gdim = 2
    mesh_comm = MPI.COMM_WORLD
    model_rank = 0
    if not my_globals['gmsh_initialized']:
        gmsh.initialize()
        my_globals['gmsh_initialized']=True
    facet_tags = {"left": 1, "right": 2, "top": 3, "crack": 4, "bottom_no_crack": 5}
    cell_tags = {"all": 20}
    if mesh_comm.rank == model_rank:
        model = gmsh.model()
        model.add("Rectangle")
        model.setCurrent("Rectangle")
        # Create the points
        p1 = model.geo.addPoint(0.0, 0.0, 0, lc)
        p2 = model.geo.addPoint(Lcrack, 0.0, 0, lc)
        p3 = model.geo.addPoint(Lx, 0, 0, lc)
        p4 = model.geo.addPoint(Lx, Ly, 0, lc)
        p5 = model.geo.addPoint(0, Ly, 0, lc)
        # Create the lines
        l1 = model.geo.addLine(p1, p2, tag=facet_tags["crack"])
        l2 = model.geo.addLine(p2, p3, tag=facet_tags["bottom_no_crack"])
        l3 = model.geo.addLine(p3, p4, tag=facet_tags["right"])
        l4 = model.geo.addLine(p4, p5, tag=facet_tags["top"])
        l5 = model.geo.addLine(p5, p1, tag=facet_tags["left"])
        # Create the surface
        cloop1 = model.geo.addCurveLoop([l1, l2, l3, l4, l5])
        surface_1 = model.geo.addPlaneSurface([cloop1])
        
        # Define the mesh size and fields for the mesh refinement
        model.mesh.field.add("Distance", 1)
        model.mesh.field.setNumbers(1, "NodesList", [p2])
        # SizeMax -                   / ------------------
        #                            /
        # SizeMin -o----------------/
        #          |                |  |
        #        Point        DistMin   DistMax
        model.mesh.field.add("Threshold", 2)
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
    top_facets = mesh.locate_entities_boundary(msh, 1, lambda x : np.isclose(x[1], Ly))
    mt = mesh.meshtags(msh, 1, top_facets, 1)
    ds = ufl.Measure("ds", subdomain_data=mt)
    if verbosity>1:
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
    if verbosity>0:
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