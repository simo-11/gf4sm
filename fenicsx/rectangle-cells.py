# %%
import importlib
import rectangle
import membrane
import timeit
import numpy as np
import dolfinx
from mpi4py import MPI
# %%
importlib.reload(rectangle)
# %%
l=10
h=1
w=1
nx=10
ny=2
nz=2
mesh = dolfinx.mesh.create_box(MPI.COMM_WORLD, 
                               [np.array([0,0,0]), np.array([l, h, w])],
                [nx,ny,nz], cell_type=dolfinx.mesh.CellType.hexahedron)
V = dolfinx.fem.VectorFunctionSpace(mesh, ("CG", 1))
f = rectangle.TorsionLoad(mesh=mesh)
x=np.array([1, 2, 3])
print(f(x))
# %%
timeit.timeit(lambda:rectangle.main(order=1,nx=100,ny=2,loadCase=rectangle.Load.Torsion),number=1)
# %%
timeit.timeit(lambda:rectangle.main(order=1,nx=100,ny=2),number=1)
# %%
rectangle.main(order=3)
# %%
import importlib
importlib.reload(membrane)
# %%
# %%
membrane.mesh()
# %%
membrane.main()

# %%
