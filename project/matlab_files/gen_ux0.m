function u_x0 = gen_ux0(n,cycles)
tlimit=cycles*2*pi;tlen=n; tstepsize=tlimit/tlen;
xlimit=cycles*2*pi;xlen=n; xstepsize=xlimit/xlen;
tspan=linspace(0,tlimit,tlen);
xmesh=linspace(0,xlimit,xlen);
% For recovering g(x) we meed f(0)!=0
% Define f(t)=cos(t)
f=cos(tspan);
% f=exp(-tspan);
% Define g(x)=sin(2*pi*x)
g=sin(2*pi*xmesh);
% g=exp(-xmesh).*sin(2*pi*xmesh);

sig = conv(f,g);
sig = sig(1:n);

% normalize u_x0
% u_x0 = u_x0 ./ max(u_x0);
u_x0 = sig/max(sig);

end