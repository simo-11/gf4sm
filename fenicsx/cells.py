# %%
import importlib
import rectangle
import timeit
# %%
importlib.reload(rectangle)
# %%
timeit.timeit(lambda:rectangle.main(order=1,nx=100,ny=2),number=1)
# %%
rectangle.main(order=3)