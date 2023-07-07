# gf4sm
Studies for using general fem packages for structural mechanics

# Julia notes
## Julia update
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
### 1.8.5 -> 1.9.0
Update fails for these but torsion works
 * P4est_wrapper
 * GridapP4est
 * GridapPETSc

### 1.9.2
code integrarion broken, possibly too early update. 

Workflow updated so that only editing is done in vs-code and code is executed on distinct julia REPL.

PATH updated

## Using Debugger
 * https://julialang.org/blog/2019/03/debuggers/

```
julia> using Debugger
julia> @enter Mesh.main(domain=domain,partition=partition)
1|debug> ?
  Debugger commands...
```


# References
 * [Gridap](https://github.com/gridap/Gridap.jl)
   * https://gridap.github.io/Tutorials/dev/
   * [The software design of Gridap](https://arxiv.org/pdf/2109.12818v1.pdf)
