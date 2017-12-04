%% THIS PART IS THE FORWARD PROBLEM
% We have that w(x,t) = f(t-x) + integral(k(x,s)*f(t-s)ds) from s->t for
% 0<x<=t
% Assuming that q(x) = 0 at first for simplicity which causes k(x,s) = 0;
tlimit=2*pi;tlen=100; tstepsize=tlimit/tlen;
xlimit=2*pi;xlen=100; xstepsize=xlimit/xlen;
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

% We want to test f(x) and g(x) that are from the sinusoidal and
% decaying exponential families.

% For this noise, calculate the power of the signal then choose
% a range of SNR.  The goal is to have an exponential curve of error % vs
% SNR
% Loop over SNR and collect data, then plot MSE vs SNR
idx=1; snr_lim=30; snr_step=.2; clear mean_squared_err;
for snr=.1:snr_step:snr_lim;
    u_x0_noisy = awgn(u_x0,snr,'measured');
    
    %% THIS PART STARTS THE INVERSE PROBLEM
    % this solves the following volterra equation for G(t):
    %   u_x(0,t)= f(0)G(t) + integral(f'(t-x)G(x)dx)(0->t)
    % then we differentiate G(t) to get g(t)=g(x) for t<x
    
    % calculate f'(t) = (f(x+h)-f(x))/h
    df=gradient(f)/xstepsize;
    
    % rearrange the above volterra eq for G(t)
    %   G(t)= u_x(0,t)/f(0) + integral(-f'(t-x)G(x)dx)(0->t)
    [G,tt]=voltra(tlen, 1, 0, tstepsize, u_x0_noisy./f(1), -df./f(1));
    
    % get g(t) from G(t)
    recovered_g=gradient(G);
    drg=gradient(recovered_g); %drg is our actual recovered g
   
    mean_squared_err(idx)=immse(drg,g);
    idx=idx+1;
end
snr=.1:snr_step:snr_lim;

plot(snr,mean_squared_err,'.'); 
xlabel('SNR');
ylabel('MSE');
str=sprintf('MSE vs SNR for f(t)=cos(t), g(x)=sin(2*pi*x), q(x)=0; 0.1<=SNR<=%d',snr_lim);
title(str)