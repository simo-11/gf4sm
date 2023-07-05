module Mesh
using Gridap
using LinearAlgebra
using Printf
  function main(;domain,partition)
    model = CartesianDiscreteModel(domain,partition)
    labeling = get_face_labeling(model)
    fn=@sprintf("../paraview/%dd_mesh_model",length(partition))
    vtk_file=writevtk(model,fn)
    domain,partition,model,vtk_file
  end
end
