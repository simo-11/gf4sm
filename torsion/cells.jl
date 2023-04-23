## In visual studio code: Run selected cell by ALT-Enter
# https://www.julia-vscode.org/docs/dev/userguide/runningcode/#The-Julia-REPL@enter println("Test")
## reload rectangle
include("rectangle.jl")
## run rectangle simulation
Rectangle.main(n=6,nsteps=20)
##