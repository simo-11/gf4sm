%% path and init Formatting
p=path;
if ~contains(p,'chebfun')
    addpath("..\..\chebfun");
end
ODEformats
%% Set parameters
E=210E9;
rho=7800;
g=9.8;
b=0.1;
h=0.1;
I=b*h^3/12;
l=3;
q=-b*h*rho*g;
ymax=q*l^4/(8*E*I);
mmax=q*l^2/2;
fprintf("Analytical:\nymax=%.2g\nmmax=%.2g\n",ymax,mmax);
%% Solve using chebop
L = chebop(0,l); % domain 0..L
L.op = @(x,y) E*I*diff(y,4)-q; % Euler-Bernoulli
L.lbc = [0;0]; % zero displacement and slope at fixed end
L.rbc=@(y)[diff(y,2);diff(y,3)]; % zero moment and force at free end
y = L\-1; % solve
plot(y,CO,bvp,LW,3); 
%axis([0 l -5e-3 0])
ymax=y(l);
M=E*I*diff(y,2);
fprintf("Chebfun:\nymax=%.2g\n",ymax);
fprintf("moment(0)=%.2g\n",M(0));
title_text=sprintf("Length=%.2g ymax=%.2g",l,ymax);
title(title_text)
