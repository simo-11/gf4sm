%% init path
% https://en.wikipedia.org/wiki/Euler%E2%80%93Bernoulli_beam_theory
p=path;
if ~contains(p,'chebfun')
    addpath("..\..\chebfun");
end
%% cantilevered beam with constant load
p=parameters;
p.case=1;
fprintf("Case %d: Constant load q=%.3g, clamped and free ends\n",...
        p.case,p.q);
tic
L = chebop(0,p.l); % domain 0..L
L.op = @(x,y) p.E*p.I*diff(y,4); % Euler-Bernoulli
L.lbc = [0;0]; % zero displacement and slope at fixed end
L.rbc=@(y)[diff(y,2);diff(y,3)]; % zero moment and force at free end
y = L\p.q; % solve
if p.toc
    toc
end
ymax=p.q*p.l^4/(8*p.E*p.I);
mmax=p.q*p.l^2/2;
report_results(y,p,p.l,0,ymax,mmax)
%% above but two simple supports at left end
p=parameters;
p.case=2;
fprintf("Case %d: Constant load q=%.3g, %s\n",...
    p.case,p.q,"two simple support and free end");
tic
dom = [0 p.l];
N = chebop(@(x,u) p.E*p.I*diff(u,4), dom);
N.bc = @(x,u) [u(0)
u(p.l/2000)    
feval(diff(u,2),p.l)
feval(diff(u,3),p.l)
];
rhs = p.q;
y = solvebvp(N, rhs);
if p.toc
    toc
end
ymax=p.q*p.l^4/(8*p.E*p.I);
mmax=p.q*p.l^2/2;
report_results(y,p,p.l,0,ymax,mmax)
%% symmetry at right end
p=parameters;
p.case=3;
fprintf("Case %d: Constant load q=%.3g, symmetry\n",...
    p.case,p.q);
tic
dom = [0 p.l];
N = chebop(@(x,u) p.E*p.I*diff(u,4), dom);
N.bc = @(x,u) [u(0)
feval(diff(u),p.l)    
feval(diff(u,2),0)
feval(diff(u,3),p.l)
];
rhs = p.q;
y = solvebvp(N, rhs);
if p.toc
    toc
end
ymax=5*p.q*(2*p.l)^4/(384*p.E*p.I);
mmax=-p.q*(2*p.l)^2/8;
report_results(y,p,p.l,p.l,ymax,mmax)
%% cantilever with free end supported on a roller
p=parameters;
p.case=4;
fprintf("Case %d: %s q=%.3g %s\n",...
        p.case,"Constant load",p.q,...
        "cantilever with free end supported on a roller");
tic
L = chebop(0,p.l);
L.op = @(x,y) p.E*p.I*diff(y,4);
L.lbc = [0;0]; % zero displacement and slope at fixed end
L.rbc=@(y)[y;diff(y,2)]; % zero displacement and moment at right end
y = L\p.q; % solve
if p.toc
    toc
end
ymax=(39+55*sqrt(33))*p.q*p.l^4/(65536*p.E*p.I);
mmax=p.q*p.l^2/8;
report_results(y,p,((15-sqrt(33))/16)*p.l,0,ymax,mmax)
%% simply supported beam with asymmetric load
% 
p=parameters;
rtol=0.011; % higher rtol is needed here unless better method is found
p.case=5;
P=p.q*p.l;
a=0.8*p.l;
b=p.l-a;
dx=0.02;
fprintf("Case %d: %s P=%.3g at %.3g(%.3g-%.3g) %s\n",...
        p.case,"Point load",P,a,a-dx,a+dx,...
        "with simple supports");
tic
dom = [0 p.l];
N = chebop(@(x,u) p.E*p.I*diff(u,4), dom);
N.bc = @(x,u) [u(0)
u(p.l)    
feval(diff(u,2),0)    
feval(diff(u,2),p.l)    
];
q=P/(2*dx);
rhs = chebfun({@(x)0,q,0},[0 a-dx a+dx p.l]);
y = solvebvp(N, rhs);
if p.toc
    toc
end
if p.plot || 1
    figure(100+p.case)
    plot(rhs)
    title('loading')
end
x_for_max_d=sqrt((p.l^2-b^2)/3);
ymax=sqrt(3)*P*b*(p.l^2-b^2)^(3/2)/(27*p.l*p.E*p.I);
mmax=-P*a*b/p.l;
report_results(y,p,x_for_max_d,a,ymax,mmax,rtol)
%% local functions
function p=parameters
    %% Set parameters
    p.E=210E9;
    p.rho=7800;
    p.g=9.8;
    p.b=0.1;
    p.h=0.1;
    p.I=p.b*p.h^3/12;
    p.l=3;
    p.q=-p.b*p.h*p.rho*p.g;
    p.case=0;
    p.plot=0;
    p.verify=0;
    p.toc=0;
end
function report_results(y,p,x_for_max_d,x_for_max_m,...
    expected_ymax,expected_mmax,rtol)
    arguments
       y
       p
       x_for_max_d
       x_for_max_m
       expected_ymax
       expected_mmax
       rtol=0.001
    end
    import matlab.unittest.TestCase
    import matlab.unittest.constraints.IsEqualTo
    import matlab.unittest.constraints.RelativeTolerance
    ymax=y(x_for_max_d);
    M=p.E*p.I*diff(y,2);
    mmax=M(x_for_max_m);
    fprintf("Case %d: ymax=%.3g(%.3g) mmax=%.3g(%.3g) rtol=%.2g%%\n",...
        p.case,ymax,expected_ymax,mmax,expected_mmax,100*rtol);
    if p.plot
        figure(10*p.case+1)
        plot(y)
        title_text=sprintf("case %d y(%.3g)=%.3g",p.case,x_for_max_d,ymax);
        title(title_text)
        figure(10*p.case+2)
        plot(M)
        title_text=sprintf("case %d M(%.3g)=%.3g",p.case,x_for_max_m,mmax);
        title(title_text)
        Q=diff(M);
        figure(10*p.case+3)
        plot(Q)
        title_text=sprintf("case %d Q",p.case);
        title(title_text)
    end
    if p.verify
        testCase = TestCase.forInteractiveUse;
        testCase.verifyThat(ymax,IsEqualTo(expected_ymax, ...
        "Within",RelativeTolerance(rtol)))
        testCase.verifyThat(mmax,IsEqualTo(expected_mmax, ...
        "Within",RelativeTolerance(rtol)))
    end
end