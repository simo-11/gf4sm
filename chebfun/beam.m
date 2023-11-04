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
sx=0.0001*l;
q=-b*h*rho*g;
%% constant load, fixed and free ends
tic
L = chebop(0,l); % domain 0..L
L.op = @(x,y) E*I*diff(y,4); % Euler-Bernoulli
L.lbc = [0;0]; % zero displacement and slope at fixed end
L.rbc=@(y)[diff(y,2);diff(y,3)]; % zero moment and force at free end
y = L\q; % solve
toc
figure(11)
%axis([0 l -5e-3 0])
ymax=y(l);
M=E*I*diff(y,2);
fprintf("Constant load q=%.3g, fixed and free ends\n",q);
fprintf("%10s: ymax=%.3g mmax=%.3g\n","Chebfun",ymax,M(0));
title_text=sprintf("y(l)=%.3g",ymax);
title(title_text)
ymax=q*l^4/(8*E*I);
mmax=q*l^2/2;
fprintf("%10s: ymax=%.3g mmax=%.3g\n","Analytical",ymax,mmax);
%% above + simple support at sx - free rotation at 0
tic
dom = [0 l];
N = chebop(@(x,u) E*I*diff(u,4), dom);
N.bc = @(x,u) [u(0)
u(sx)    
%feval(diff(u),0)
feval(diff(u,2),l)
feval(diff(u,3),l)
];
rhs = q;
y = solvebvp(N, rhs);
figure(12)
plot(y)
toc
