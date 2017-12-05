function [p,W]=SRpoi(data,a,s)
% Poisson SR parameters
n=1;
Wstop=700;
del=3;
rho=(a+del*sqrt(abs(a)))./a;

% SR statistic calculation
W=zeros(length(data),1);
for m=2:length(data)
    s1=0;
    wm=0;
    for i=1:m-1
        s1=s1+data(m-i+1);
        wm=wm+exp(-i*del*sqrt(abs(a(m)))+log(rho(m))*s1);
    end
    W(m)=wm;
end

% Alarm location
for p=1:size(W,1)
    if W(p)>Wstop
        p;
        return;
    end
end