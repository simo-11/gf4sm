%% poisson example
syms f(x,y)
f(x,y)=1+x^2+2*y^2 %#ok<NOPTS> %
lf=laplacian(f) %#ok<NOPTS>
%mlf=diff(f,x,2)+diff(f,y,2) %#ok<NOPTS>
%diff(f,x,1) 
%diff(f,x,2) 
%diff(f,y,1)
%diff(f,y,2)
%% 

