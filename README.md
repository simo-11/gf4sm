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

 P4est_wrapper is due to issue of build.jl not in detecting libp4est.dll.a as DLL library
 ```
 ERROR: Error building `P4est_wrapper`:
[ Info: p4est directory found at: D:\julia\artifacts\89a337ea6f60a4fd58999ab73dea099e41032138
[ Info: p4est lib directory found at: D:\julia\artifacts\89a337ea6f60a4fd58999ab73dea099e41032138\lib
┌ Warning: P4EST lib (.dll) not found at: D:\julia\artifacts\89a337ea6f60a4fd58999ab73dea099e41032138\lib
└ @ Main D:\julia\packages\P4est_wrapper\v17J6\deps\build.jl:44
[ Info: P4EST include directory found at: D:\julia\artifacts\89a337ea6f60a4fd58999ab73dea099e41032138\include
┌ Info: P4EST configuration:
│ ==============================================
│   - P4EST_FOUND           = false
│   - P4EST_DIR             = D:\julia\artifacts\89a337ea6f60a4fd58999ab73dea099e41032138
│   - P4EST_LIB_DIR         = D:\julia\artifacts\89a337ea6f60a4fd58999ab73dea099e41032138\lib
│   - P4EST_INCLUDE_DIR     = D:\julia\artifacts\89a337ea6f60a4fd58999ab73dea099e41032138\include
│   - P4EST_LIB             = .dll
└   - P4EST_ENABLE_MPI      = false
ERROR: LoadError: P4EST library not found
```

Attemp to fix
```
pkg> add P4est_jll # did not help

julia> ENV["P4EST_LIB"]="D:/julia/artifacts/89a337ea6f60a4fd58999ab73dea099e41032138/lib/libp4est.dll.a"
julia> Pkg.build("Tutorials")

julia> import GridapMakie
[ Info: Precompiling GridapMakie [41f30b06-6382-4b60-a5f7-79d86b35bf5d]
ERROR: LoadError: UndefVarError: `point_iterator` not defined # commented it out
Stacktrace:
 [1] getproperty(x::Module, f::Symbol)
   @ Base .\Base.jl:31
 [2] top-level scope
   @ D:\julia\packages\GridapMakie\oS6vj\src\recipes.jl:158
```

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
