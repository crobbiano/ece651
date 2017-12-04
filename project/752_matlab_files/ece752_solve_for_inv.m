clear all
clc
addpath('./fortan_conversions/')

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
snr=15;
u_x0_noisy = awgn(u_x0,snr,'measured');

% plotting the calculated u_x(0,t) and noisy version and some other stuff
str=sprintf('u_x(0,t) and u_x(0,t) with noise added, SNR=%d',snr);
figure(1); subplot(3,1,1); plot(tspan,u_x0_noisy);title(str);xlabel('t');hold on; plot(tspan, u_x0);grid minor; legend('noise','noiseless');hold off; xlim([0 2*pi]);
% figure(1); subplot(3,1,1);  plot(tspan, u_x0);title('Generated u_x(0,t)');xlabel('time');grid minor; xlim([0 2*pi]);
subplot(3,1,2); plot(xmesh,g);title('g(x)=sin(2*pi*x) versus x'); xlabel('x'); grid minor; xlim([0 2*pi]);
subplot(3,1,3); plot(tspan,f);title('w(0,t)=f(t)=cos(t)'); xlabel('time');grid minor; xlim([0 2*pi]);


%% THIS PART STARTS THE INVERSE PROBLEM
% this solves the following volterra equation for G(t):
%   u_x(0,t)= f(0)G(t) + integral(f'(t-x)G(x)dx)(0->t)
% then we differentiate G(t) to get g(t)=g(x) for t<x

% calculate f'(t) = (f(x+h)-f(x))/h
df=gradient(f)/xstepsize;

% plotting pretty pictures of f(x) and f'(x)
figure(2);subplot(2,1,1);plot(xmesh,f);ylabel('f(t)');grid minor;xlim([0 2*pi]);
title('visualizing the derivitives');
subplot(2,1,2);plot(xmesh,df);ylabel('df(t)');xlabel('time');grid minor;xlim([0 2*pi]);

% rearrange the above volterra eq for G(t)
%   G(t)= u_x(0,t)/f(0) + integral(-f'(t-x)G(x)dx)(0->t)
% [G,tt]=voltra(tlen, 1, 0, tstepsize, u_x0./f(1), -df./f(1));
[G,tt]=voltra(tlen, 1, 0, tstepsize, u_x0_noisy./f(1), -df./f(1));

% get g(t) from G(t)
recovered_g=gradient(G);
drg=gradient(recovered_g); %drg is our actual recovered g

% Plotting pretty pictures
figure(3); subplot(3,1,1);plot(tspan,G);
title('Recovered G(t)'); grid minor;xlim([0 tlimit]);
subplot(3,1,2); plot(tspan,recovered_g);
title('Recovered g(t)');xlabel('t'); grid minor;xlim([0 tlimit]);
subplot(3,1,3); plot(tspan,drg);
title('Recovered dg(t)');xlabel('t'); grid minor;xlim([0 tlimit]);

mean_squared_err=immse(drg,g);
str=sprintf('Plot of recovered g versus known g for SNR=%d; MSE=%f',snr,mean_squared_err);
figure(4);  plot(tspan,g); 
hold on;plot(tspan,drg);hold off; 
xlabel('time');legend('known g','recovered g');title(str);grid minor;

% str=sprintf('u_x(0,t) and u_x(0,t) with noise added, SNR=%d',snr);
% figure(5);subplot(2,1,1);plot(tspan,u_x0_noisy);title(str);xlabel('time');hold on; plot(tspan, u_x0);grid minor; legend('noise','noiseless');hold off; xlim([0 2*pi]);
% str=sprintf('Plot of recovered g versus known g for SNR=%d; MSE=%f',snr,mean_squared_err);
% subplot(2,1,2);plot(tspan,g); hold on;plot(tspan,drg);hold off; xlabel('time');legend('known g','recovered g');title(str);grid minor;

% str=sprintf('u_x(0,t)');
% figure(6);subplot(2,1,1);plot(tspan,u_x0_noisy);title(str);xlabel('time');hold on; plot(tspan, u_x0);grid minor; legend('noise','noiseless');hold off; xlim([0 2*pi]);
% str=sprintf('Plot of recovered g versus known g; MSE=%f',mean_squared_err);
% subplot(2,1,2);plot(tspan,g); hold on;plot(tspan,drg);hold off; xlabel('time');legend('known g','recovered g');title(str);grid minor;
