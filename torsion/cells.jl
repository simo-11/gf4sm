## In visual studio code: Run selected cell by ALT-Enter
# https://www.julia-vscode.org/docs/dev/userguide/runningcode/
# The-Julia-REPL
# @enter println("Test")
## (re)load rectangle
include("rectangle.jl")
## run rectangle simulation
using Gridap
force = VectorValue(0,-7800*9.8,0)#gravity
# l>7 will fails to converge using default parameters to NLSolver
Rectangle.main(l=1,n=1,nsteps=1,order=1,load=force)
#
## (re)load mesh
include("mesh.jl")
## 1d mesh
domain=(0,1)
partition=(2)
model,reffe,V,U,omega,d_omega,solver,uh=
    Mesh.main(domain=domain,partition=partition,order=1,load=1)
## 2d mesh
domain=(0,1,0,1)
partition=(2,2)
model,reffe=Mesh.main(domain=domain,partition=partition,order=1)
## 3d mesh
domain=(0,1,0,1,0,1)
partition=(2,2,2)
model,reffe=Mesh.main(domain=domain,partition=partition,order=1)
