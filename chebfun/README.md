# chebfun
https://www.chebfun.org/ is required for these examples.
Clone https://github.com/chebfun/chebfun.git beside gf4sm to get started easily.

Modify parameters function to test other values and get plots, timings and verifications. 

## Bending
Pure [Bending](beam.m) problems work quite nicely.

## Torsion
[Torsion](torsion.m) does not seem usable. If there is better way, please let me know.
Rotation gets quite correct value

![image](https://github.com/simo-11/gf4sm/assets/1210784/3aad5281-f943-42e1-a5c5-b3c2e8c1f11e)

but warping (derivate of rotation) oscillates badly.
Quick convergence to 0.32 is expected

![image](https://github.com/simo-11/gf4sm/assets/1210784/cac909fe-6bbe-43f5-a46d-bb9c8744515e) 

but constant oscillation between 0 and 0.66 is obtained

![image](https://github.com/simo-11/gf4sm/assets/1210784/fb8f5c47-32a8-4865-99c7-aed0019c9346)



