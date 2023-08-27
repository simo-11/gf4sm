# %%
import importlib
import rectangle
import membrane
import timeit
# %%
importlib.reload(rectangle)
# %%
timeit.timeit(lambda:rectangle.main(order=1,nx=100,ny=2,loadCase=rectangle.Load.Torsion),number=1)
# %%
timeit.timeit(lambda:rectangle.main(order=1,nx=100,ny=2),number=1)
# %%
rectangle.main(order=3)
# %%
import importlib
import warping
import timeit
# %%
import importlib
importlib.reload(warping)
# %%
timeit.timeit(lambda:warping.main(order=1,nx=100,ny=2),number=1)
# %%
warping.main(order=1)
# %%
import importlib
importlib.reload(membrane)
# %%
# %%
membrane.mesh()
# %%
membrane.main()

# %%
