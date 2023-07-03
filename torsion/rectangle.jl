module Rectangle
using Gridap
using LinearAlgebra
const E = 210e9 # Pa - Steel
const nu = 0.3 # dim-less
const lambda = (E*nu)/((1+nu)*(1-2*nu))
const my = E/(2*(1+nu))
sigma_e(e) = lambda*tr(e)*one(e) + 2*my*e # Pa
tau(e) = sqrt(e ⊙ sigma_e(e)) # Pa^(1/2), ⊙ = inner, Gridap: src/TensorValues/Operations.jl
const sigma_u = 800e6 # Pa - High Strength
const r_0 = sigma_u / sqrt(E) # Pa^(1/2)
const H = 0.5 # dim-less, find reference

function d(r)
  1 - q(r)/r
end

function q(r)
  r_0 + H*(r-r_0)
end

function new_state(r_in,d_in,e_in)
    tau_in = tau(e_in)
    if tau_in <= r_in
      r_out = r_in
      d_out = d_in
      damaged = false
    else
      r_out = tau_in
      d_out = d(r_out)
      damaged = true
    end
    damaged, r_out, d_out
  end

  function sigma(e_in,r_in,d_in)
    _, _, d_out = new_state(r_in,d_in,e_in)
    (1-d_out)*sigma_e(e_in)
  end
  
  function d_sigma(de_in,e_in,state)
    damaged, r_out, d_out = state
    if ! damaged
      return (1-d_out)*sigma_e(de_in)
    else
      c_inc = ((q(r_out) - H*r_out)*(sigma_e(e_in) ⊙ de_in))/(r_out^3)
      return (1-d_out)*sigma_e(de_in) - c_inc*sigma_e(e_in)
    end
  end

  function project(q,model,d_omega,order)
    reffe = ReferenceFE(lagrangian,Float64,order)
    V = FESpace(model,reffe,conformity=:L2)
    a(u,v) = ∫( u*v )*d_omega
    l(v) = ∫( v*q )*d_omega
    op = AffineFEOperator(a,l,V,V)
    qh = solve(op)
    qh
  end

  function main(;r,n,nsteps,order,load)
    domain = (0,r,0,1,0,1)
    partition = (r*n,n,n)
    model = CartesianDiscreteModel(domain,partition)
    labeling = get_face_labeling(model)
    # entities at lower x values
    # 1,3,5,7 are corners
    # 13,15,17,19 edges (excluding corners)
    # 25 is center area (excluding corners and edges)
    x0c=range(start=1,step=2,length=4)
    x0e=range(start=13,step=2,length=4)
    x0f=25
    fixed=vcat(x0c,x0e,x0f) 
    add_tag_from_tags!(labeling,"fixed",fixed)
    reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
    V = TestFESpace(model,reffe,labels=labeling,dirichlet_tags=["fixed"])
    U = TrialFESpace(V)
    degree = 2*order
    omega = Triangulation(model)
    d_omega = Measure(omega,degree)
    r = CellState(r_0,d_omega)
    d = CellState(0.0,d_omega)
    nls = NLSolver(show_trace=true, method=:newton)
    solver = FESolver(nls)
    function step(uh_in,factor,cache,b_max)
      b = factor*b_max
      res(u,v) = ∫(  ε(v) ⊙ (sigma∘(ε(u),r,d))  - v⋅b )*d_omega
      jac(u,du,v) = ∫(  ε(v) ⊙ (d_sigma∘(ε(du),ε(u),new_state∘(r,d,ε(u))))  )*d_omega
      op = FEOperator(res,jac,U,V)
      uh_out, cache = solve!(uh_in,solver,op,cache)
      update_state!(new_state,r,d,ε(uh_out))
      uh_out, cache
    end
    factors = collect(1:nsteps)*(1/nsteps)
    uh = zero(V)
    cache = nothing
    println("Maximum load=$load, order=$order")
    for (istep,factor) in enumerate(factors)
      println("\n+++ Solving for load factor $factor in step $istep of $nsteps +++\n")
      uh,cache = step(uh,factor,cache,load)
      dh = project(d,model,d_omega,order)
      rh = project(r,model,d_omega,order)
      writevtk(
        omega,"../paraview/torsion_$(lpad(istep,3,'0'))",
        cellfields=["uh"=>uh,"epsi"=>ε(uh),"damage"=>dh,
                    "threshold"=>rh,"sigma_elast"=>sigma_e∘ε(uh)])
    end
  end
end
