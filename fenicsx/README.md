# Fenicsx

# Running cases

https://jsdokken.com/dolfinx-tutorial/chapter2/linearelasticity_code.html
```
simo@MSI:~/github/gf4sm/fenicsx$ python3 tutorial_lineareleasticity.py
Wrote results to ../paraview/tutorial_linearelasticity_deformation.xdmf
```

https://docs.fenicsproject.org/dolfinx/main/python/demos/demo_elasticity.html
```
simo@MSI:~/github/gf4sm/fenicsx$ python3 demo_elasticity.py
```
## Issues

### AdjacencyList_int32 object has no attribute flatten
```
  File "/home/simo/github/gf4sm/fenicsx/demo_elasticity.py", line 66, in <listcomp>
    dofs = [V.sub(i).dofmap.list.flatten() for i in range(3)]
AttributeError: 'dolfinx.cpp.graph.AdjacencyList_int32' object has no attribute 'flatten'
```


# Developing and debugging
```
simo@MSI:~/github/gf4sm/fenicsx$ code .
```

# Links
 * Tutorial - https://jsdokken.com/dolfinx-tutorial/
 * Documentation - https://docs.fenicsproject.org/dolfinx/main/python/
 * Unified Form Language - creates c-code based on python
   * Partial element list in form language - https://fenics.readthedocs.io/projects/ufl/en/latest/manual/form_language.html
 * Book on fenics (old version) - https://launchpadlibrarian.net/83776282/fenics-book-2011-10-27-final.pdf
   * Applications in solid mechanics starting at page 505(517 in pdf)
 * Fracture mechanics examples https://newfrac.github.io/fenicsx-fracture/
 * dolphiny high-level wrappers for DOLFINx - https://github.com/michalhabera/dolfiny/
 * Discource
   * https://fenicsproject.discourse.group/t/plasticity-in-fenicsx/11855 
   

# Installation on WSL Ubuntu
Based on https://github.com/FEniCS/dolfinx#ubuntu-packages, 453 packages, 2169 MB
```
root@MSI:/home/simo# apt install fenicsx
root@MSI:/home/simo# apt install gmsh python3-gmsh
```
