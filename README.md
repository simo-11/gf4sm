# gf4sm
Studies for using general fem packages for structural mechanics

## Bending

[Beam theory based results](https://docs.google.com/spreadsheets/d/1350EPOZFU3kTkPZUV8PogOySmf1ne2bMS18DN-5WTs4/)
 * E=210 GPa
 * rho=7800 kg/m^3
 * nu=0.3 (Not effective in beam theory)

100x100 mm (0.1x0.1 m) solid, own weight, L= 3 m, All dofs for all nodes locked at fixed end.

| Simulation | deflection at tip [m] | Solve time [s]|
|------------|-----------------------|----------------|
|  | g/f | g/f |   
| Beam theory | 0.0044| |
| b_1_9_3_3 | 0.0008/0.0008 | 0.05/0.11 | 
| b_1_30_1_1 | 0.0028 | 0.1 | 
| b_1_400_1_1 | 0.0038 | 0.1 | 
| b_1_100_4_2 | 0.0042/0.0042 | 0.19/0.12 |
| b_1_400_4_2 (nu=0) | 0.0044/0.00424 | 1.0/0.13 |
| b_2_1_1_1 | 0.0029/0.0029 | 0.16/0.1 | 
| b_2_2_1_1 | 0.0039 | 0.18 | 
| b_2_3_1_1 | 0.0041 | 0.17 | 
| b_2_4_1_1 | 0.0042 | 0.18 | 
| b_2_5_1_1 | 0.0043 | 0.18 | 
| b_3_1_1_1 | 0.0041 | 1.5 |
| b_3_2_1_1 | 0.0043/0.00425 | 1.6/0.13 |

cell notation
 * g_ means Gridap(Julia), f_ means fenicsx(Python)
 * b_ means Bending
 * next number is order.
   * First order elements are not efficient for bending but large amounts give good results
 * next numbers nx, ny and nz are number of elements in each direction

Gridap solve time is about doubled if material constants are variables. 
This can be handled by refactoring code in same ways as is done in python.
Gridap (as of 2023/07) is quite sensitive to increasing of number or order of elements. 
Runtime and memory allocations increase heavily.

### Gridap sample runs
See also (Sample command cells)[julia/cells.jl]
```
julia> 
julia> Rectangle.setNu(0)
julia> @time _,=Rectangle.main(l=3,h=0.1,w=0.1,nx=300,ny=4,nz=2,order=1,load=force)
Maximum displacement 0.0044
  1.026924 seconds (16.17 M allocations: 938.163 MiB, 12.75% gc time)
```

### Fenicsx sample runs
See also (Sample command cells)[fenicsx/cells.py]
```
>>> timeit.timeit(lambda:rectangle.main(order=3,nx=2,ny=1,nz=1),number=1)
WARNING:py.warnings:/usr/lib/python3/dist-packages/ffcx/element_interface.py:82: UserWarning: Number of integration points per cell is: 125. Consider using 'quadrature_degree' to reduce number.
```
