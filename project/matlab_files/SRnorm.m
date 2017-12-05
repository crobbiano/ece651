function [p,W]=SRnorm(data,a,s)
% Gaussian SR parameters
n=1;
Wstop=700;
del=1;
d2=del^2;
s2=s^2;

% SR statistic calculation
W=zeros(length(data),1);
for m=2:length(data)
    wm=0;
    sum=0;
    for i=1:m-1
        sum=sum+data(m-i+1)-a(m-i+1);
        wm=wm+exp(-i*n*d2/(2*s2)+n*del*sum/s2);
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
end