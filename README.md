# gf4sm
Studies for using general fem packages for structural mechanics

## Bending

(Beam theory based results)[https://docs.google.com/spreadsheets/d/1350EPOZFU3kTkPZUV8PogOySmf1ne2bMS18DN-5WTs4/]
 * E=210 GPa
 * rho=7800 kg/m^3
 * nu=0.3 (Not effective in beam theory)

100x100 mm (0.1x0.1 m) solid, own weight, L= 3 m
|| Simulation || deflection at tip [m] || Solve time [s]||
| Beam theory | 0.0042| |
| julia(Gridap) cell g_b_1_9_3_3 | 0.0008 | 0.05 | 
| julia(Gridap) cell g_b_1_30_1_1 | 0.0028 | 0.1 | 
| julia(Gridap) cell g_b_1_400_1_1 | 0.0038 | 0.1 | 
| julia(Gridap) cell g_b_1_100_4_2 | 0.0042 | 0.19 |
| julia(Gridap) cell g_b_2_1_1_1 | 0.0029 | 0.16 | 
| julia(Gridap) cell g_b_2_2_1_1 | 0.0039 | 0.18 | 
| julia(Gridap) cell g_b_2_3_1_1 | 0.0041 | 0.17 | 
| julia(Gridap) cell g_b_2_4_1_1 | 0.0042 | 0.18 | 
| julia(Gridap) cell g_b_2_5_1_1 | 0.0043 | 0.18 | 
| julia(Gridap) cell g_b_3_1_1_1 | 0.0041 | 1.5 |
| julia(Gridap) cell g_b_3_2_1_1 | 0.0043 | 1.6 |

