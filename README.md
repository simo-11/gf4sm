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
| Beam theory | 0.0044| |
| g_b_1_9_3_3 | 0.0008 | 0.05 | 
| g_b_1_30_1_1 | 0.0028 | 0.1 | 
| g_b_1_400_1_1 | 0.0038 | 0.1 | 
| g_b_1_100_4_2 | 0.0042 | 0.19 |
| g_b_1_400_4_2 (nu=0) | 0.0044 | 1.0 |
| g_b_2_1_1_1 | 0.0029 | 0.16 | 
| g_b_2_2_1_1 | 0.0039 | 0.18 | 
| g_b_2_3_1_1 | 0.0041 | 0.17 | 
| g_b_2_4_1_1 | 0.0042 | 0.18 | 
| g_b_2_5_1_1 | 0.0043 | 0.18 | 
| g_b_3_1_1_1 | 0.0041 | 1.5 |
| g_b_3_2_1_1 | 0.0043 | 1.6 |

cell notation
 * g_b_ means Gridap(Julia) Bending
 * next number is order.
   * First order elements are not efficient for bending but large amounts give good results
 * next numbers nx, ny and nz are number of elements in each direction

Gridap solve time is about doubled if material constants are variables.

```
julia> Rectangle.setNu(0)
julia> @time _,=Rectangle.main(l=3,h=0.1,w=0.1,nx=300,ny=4,nz=2,order=1,load=force)
Maximum displacement 0.0044
  1.026924 seconds (16.17 M allocations: 938.163 MiB, 12.75% gc time)
```

