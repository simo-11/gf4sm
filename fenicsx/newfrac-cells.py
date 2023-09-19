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
pyvista.set_jupyter_backend('client')
# %% reload newfrac
importlib.reload(newfrac)
# %% set dimensions and parameters
Lx = 1.
Ly = 0.5
Lcrack = 0.3
lc =.01*Lx
dist_min = .1
dist_max = .3
E = 21.E10 
nu = 0.3 
load=(0,-600e3)
#loadX=None
loadX=(1-0.5)*Lx-0.5*lc
refinement_ratio=2
# %% solve_elasticity with crack
crack_results=newfrac.solve_elasticity(Lx=Lx,
        Lcrack=Lcrack,
        Ly=Ly,
        lc=lc,  # characteristic length of the mesh
        refinement_ratio=refinement_ratio,  # how much it is refined at the tip zone
        dist_min=dist_min,  # radius of tip zone
        dist_max=dist_max,  # radius of the transition zone
        E=E,
        nu=nu,
        load=load,
        loadX=loadX,
        verbosity=1
    )
# %% solve_elasticity without crack
no_crack_results=newfrac.solve_elasticity(Lx=Lx,
        Lcrack=-1,
        Ly=Ly,
        lc=lc,  
        refinement_ratio=refinement_ratio, 
        dist_min=dist_min,
        dist_max=dist_max,
        E=E,
        nu=nu,
        load=load,
        loadX=loadX,
        verbosity=1
    )
# %% generate mesh
msh, mt = newfrac.generate_mesh_with_crack(Lx=Lx,
        Lcrack=Lcrack,
        Ly=Ly,
        lc=lc,  # caracteristic length of the mesh
        refinement_ratio=refinement_ratio,  # how much it is refined at the tip zone
        dist_min=dist_min,  # radius of tip zone
        dist_max=dist_max,  # radius of the transition zone
        verbosity=3
    )