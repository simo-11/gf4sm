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
import pyvista
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
load=1
refinement_ratio=3
# %% solve_elasticity
newfrac.solve_elasticity(
        Lcrack=Lcrack,
        Ly=Ly,
        lc=lc,  # caracteristic length of the mesh
        refinement_ratio=refinement_ratio,  # how much it is refined at the tip zone
        dist_min=dist_min,  # radius of tip zone
        dist_max=dist_max,  # radius of the transition zone
        E=E,
        nu=nu,
        load=load,
        verbosity=1
    )
# %% generate mesh
msh, mt = newfrac.generate_mesh_with_crack(
        Lcrack=Lcrack,
        Ly=Ly,
        lc=lc,  # caracteristic length of the mesh
        refinement_ratio=refinement_ratio,  # how much it is refined at the tip zone
        dist_min=dist_min,  # radius of tip zone
        dist_max=dist_max,  # radius of the transition zone
        verbosity=3
    )