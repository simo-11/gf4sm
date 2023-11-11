%% init path
% https://en.wikipedia.org/wiki/Euler%E2%80%93Bernoulli_beam_theory
p=path;
if ~contains(p,'chebfun')
    addpath("..\..\chebfun");
end
%% cantilever beam with torsional moment at free end based on 
% differential equation
p=parameters;
p.case=1;
fprintf("Case %d: T=%.3g, clamped and free ends\n",...
        p.case,p.T);
tic
N = chebop(0,p.l);
N.op = @(x,rx) p.G*p.It*diff(rx)-p.E*p.Iw*diff(rx,3); % 
N.bc = @(x,rx) [rx(0)
feval(diff(rx),0)
feval(diff(rx,2),p.l)
];
rhs = p.T;
rx = solvebvp(N, rhs);
if p.toc
    toc
end
rxmax=0.162;
wmax=p.T/(p.G*p.It);
report_results(rx,p,p.l,p.l,rxmax,wmax);
%% check how well chebfun can describe solution for warping torsion
p=parameters;
p.case=2;
fprintf("Case %d: T=%.3g, chebfun for known warping solution\n",...
            p.case,p.T);
w1 = chebfun(@(x) p.T/(p.G*p.It)*(tanh(p.k*p.l)*sinh(p.k*x)...
    -cosh(p.k*x)+1),[0 p.l]);
figure(10*p.case+1)
plot(w1)
w1 %#ok<NOPTS>
w2 = chebfun(@(x) p.T/(p.G*p.It)*(tanh(p.k*p.l)*sinh(p.k*x)...
    -cosh(p.k*x)+1),[0 p.l],'splitting','on');
figure(10*p.case+2)
plot(w2)
w2 %#ok<NOPTS>
%% local functions
%% https://t.co/proHJo1CZA
function p=parameters
    %% Set parameters
    p.E=71.7E9;
    p.nu=0.31;
    p.G=p.E/(2*(1+p.nu));
    p.It=114E-12; % table 2
    p.Iw=328E-18;
    p.k=sqrt((p.G*p.It)/(p.E*p.Iw));
    p.l=0.508;
    p.T=1;
    p.case=0;
    p.plot=1;
    p.verify=1;
    p.toc=1;
end
function report_results(rx,p,x_for_max_rx,x_for_max_w,...
    expected_rxmax,expected_wmax,rtol)
    arguments
       rx
       p
       x_for_max_rx
       x_for_max_w
       expected_rxmax
       expected_wmax
       rtol=0.001
    end
    import matlab.unittest.TestCase
    import matlab.unittest.constraints.IsEqualTo
    import matlab.unittest.constraints.RelativeTolerance
    rxmax=rx(x_for_max_rx);
    w=diff(rx);
    wmax=w(x_for_max_w);
    fprintf("Case %d: rxmax=%.3g(%.3g) wmax=%.3g(%.3g) rtol=%.2g%%\n",...
        p.case,rxmax,expected_rxmax,wmax,expected_wmax,100*rtol);
    if p.plot
        figure(10*p.case+1)
        plot(rx)
        title_text=sprintf("case %d rx(%.3g)=%.3g",...
            p.case,x_for_max_rx,rxmax);
        title(title_text)
        figure(10*p.case+2)
        plot(w)
        title_text=sprintf("case %d w(%.3g)=%.3g",...
            p.case,x_for_max_w,wmax);
        title(title_text)
        BM=p.E*p.Iw*diff(w);
        figure(10*p.case+3)
        plot(BM)
        title_text=sprintf("case %d BM",p.case);
        title(title_text)
    end
    if p.verify
        testCase = TestCase.forInteractiveUse;
        testCase.verifyThat(rxmax,IsEqualTo(expected_rxmax, ...
        "Within",RelativeTolerance(rtol)))
        testCase.verifyThat(wmax,IsEqualTo(expected_wmax, ...
        "Within",RelativeTolerance(rtol)))
    end
end