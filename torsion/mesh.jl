module Mesh
using Gridap
using LinearAlgebra
using Printf
function main(;domain,partition,
  nsteps=1,order=1,load,xtol=0.001,ftol=0.001,iterations=10)
  dim=length(partition)
  model = CartesianDiscreteModel(domain,partition)
  labeling = get_face_labeling(model)
  fn=@sprintf("../paraview/%dd_mesh_model",length(partition))
  vtk_file=writevtk(model,fn)
  labeling = get_face_labeling(model)
  reffe = ReferenceFE(lagrangian,VectorValue{dim,Float64},order)
  add_tag_from_tags!(labeling,"fixed",1)
  V = TestFESpace(model,reffe,labels=labeling
  ,dirichlet_tags=["fixed"]
  )
  fx(x)=VectorValue(0)
  U = TrialFESpace(V,[fx])
  degree = 2*order
  omega = Triangulation(model)
  d_omega = Measure(omega,degree)
  nls = NLSolver(show_trace=true
  , extended_trace=false
  , xtol=xtol
  , ftol=ftol
  , iterations=iterations
  , method=:newton)
  solver = FESolver(nls)
  function step(uh_in,factor,cache,b_max)
    b = factor*b_max
    res(u,v) = ∫(  ε(v) ⊙ (sigma∘(ε(u),r,d))  - v⋅b )*d_omega
    jac(u,du,v) = ∫(  ε(v) ⊙ (d_sigma∘(ε(du),ε(u),new_state∘(r,d,ε(u))))  )*d_omega
    op = FEOperator(res,jac,U,V)
    uh_out, cache = solve!(uh_in,solver,op,cache)
  end
  factors = collect(1:nsteps)*(1/nsteps)
  uh = zero(V)
  cache = nothing
  println("Domain=$domain, maximum load=$load, order=$order, xtol=$xtol, ftol=$ftol")
  for (istep,factor) in enumerate(factors)
    @printf("\nSolving for load factor %.2f in step %d of %d\n",
    factor,istep,nsteps)
    uh,cache = step(uh,factor,cache,load)
    writevtk(
      omega,"../paraview/torsion_$(lpad(istep,3,'0'))",
      cellfields=["uh"=>uh,"epsi"=>ε(uh),"sigma_elast"=>sigma_e∘ε(uh)])
  end
  model,reffe,V,U,omega,d_omega,solver,uh
end
end
