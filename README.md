# gf4sm
Studies for using general fem packages for structural mechanics

# Julia update
https://stackoverflow.com/questions/30555225/how-to-upgrade-julia-to-a-new-release
```
]add UpdateJulia
using UpdateJulia
update_julia()
```
update PATH
```
using Pkg
Pkg.update()
```
# 1.8.5 -> 1.9.0
Update fails for these but torsion works
 * P4est_wrapper
 * GridapP4est
 * GridapPETSc


# References
 * Gridap
   * https://gridap.github.io/Tutorials/dev/
