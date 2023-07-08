module Mesh
using Gridap
using LinearAlgebra
using Printf
using Match
const E = 210e9
const ν = 0.3
const λ = (E*ν)/((1+ν)*(1-2*ν))
const μ = E/(2*(1+ν))
function main(;domain,partition,
  nsteps=1,order=1,load,xtol=0.001,ftol=0.001,iterations=10,A=1,t=1)
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
  if dim===1
    σ(ε) = E*ε
    a(u,v)=A*∫( ε(v) ⊙ (σ∘ε(u)) )*d_omega
  else  
    σ(ε) = λ*tr(ε)*one(ε) + 2*μ*ε
    a(u,v)=∫( ε(v) ⊙ (σ∘ε(u)) )*d_omega
  end
  l(v)=load
  op = AffineFEOperator(a,l,U,V)
  if dim===1
    v=op.op.vector
    v[length(v)]=load
  end
  uh = solve(op)
  fn=@sprintf("../paraview/%dd_mesh_results",length(partition))
  writevtk(omega,fn,cellfields=["uh"=>uh,"epsi"=>ε(uh),"sigma"=>σ∘ε(uh)])
  model,reffe,V,U,omega,d_omega,op,uh
end
end
