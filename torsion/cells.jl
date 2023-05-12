## In visual studio code: Run selected cell by ALT-Enter
# https://www.julia-vscode.org/docs/dev/userguide/runningcode/
# The-Julia-REPL
# @enter println("Test")
## (re)load rectangle
include("rectangle.jl")
## run rectangle simulation
Rectangle.main(r=3,n=6,nsteps=1)
#