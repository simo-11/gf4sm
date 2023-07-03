## In visual studio code: Run selected cell by ALT-Enter
# https://www.julia-vscode.org/docs/dev/userguide/runningcode/
# The-Julia-REPL
# @enter println("Test")
## (re)load rectangle
include("rectangle.jl")
## run rectangle simulation
using Gridap
force = VectorValue(0,-7800*9.8,0)#gravity
# l>7 will fails to converge, using default parameters to NLSolver
Rectangle.main(l=100,n=2,nsteps=10,order=1,load=force)
#