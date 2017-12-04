function u_x0 = gen_ux0(n)
tlimit=2*pi;tlen=n; tstepsize=tlimit/tlen;
xlimit=2*pi;xlen=n; xstepsize=xlimit/xlen;
tspan=linspace(0,tlimit,tlen);
xmesh=linspace(0,xlimit,xlen);
% For recovering g(x) we meed f(0)!=0
% Define f(t)=cos(t)
f=cos(tspan);
% Define g(x)=sin(2*pi*x)
g=sin(2*pi*xmesh);
% g=exp(-xmesh).*sin(2*pi*xmesh);

% Generate the w(x,t)
w=zeros(xlen,tlen);
for t=1:tlen
    for x=1:xlen
        if (x<=t)
            w(x,t) =  f(t-x+1);
        else
            w(x,t) = 0;
        end
    end
end

% u_x(0,t)=integral(g(x)w(x,t)dt) from 0->t
% Generate the u_x0(t)
u_x0=zeros(1,tlen);
for t=2:tlen
    u_x0(t)=u_x0(t-1); % start with the previous iteration because we build over time
    for i=1:xlen
        u_x0(t)=u_x0(t) + g(i)*w(i,t);
    end
end


% normalize u_x0
u_x0 = u_x0 ./ max(u_x0);

end