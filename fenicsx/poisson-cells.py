# %% common init
import pyvista
import poisson
import importlib
pyvista.set_jupyter_backend('client')
# %% reload newfrac
importlib.reload(poisson)
# %%
poisson.dirichlet(order=2,verbose=0)
# %%
poisson.dirichlet_and_neumann(order=1,verbose=3)
# %%
