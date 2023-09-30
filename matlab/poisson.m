%% poisson example
% https://jsdokken.com/dolfinx-tutorial/chapter1/fundamentals.html
syms f(x,y)
f(x,y)=1+x^2+2*y^2 %#ok<NOPTS> %
lf=laplacian(f) %#ok<NOPTS>
%mlf=diff(f,x,2)+diff(f,y,2) %#ok<NOPTS>
%diff(f,x,1) 
%diff(f,x,2) 
% neuman bc, derivative of u
% in the outward normal direction 
% on the boundary.
% https://jsdokken.com/dolfinx-tutorial/chapter3/neumann_dirichlet_code.html
df(x,y)=-diff(f,y,1) %#ok<NOPTS>
df(x,1)
df(x,0)
%diff(f,y,2)
%% 

