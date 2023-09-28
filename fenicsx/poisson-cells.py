# %% common init
import pyvista
import poisson
import importlib
pyvista.set_jupyter_backend('client')
# %% reload newfrac
importlib.reload(poisson)
# %%
poisson.run(2)
# %%