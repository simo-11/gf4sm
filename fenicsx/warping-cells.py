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
