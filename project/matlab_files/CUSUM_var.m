function [k,c]=CUSUM_var(data,a,s,h,m)
% CUSUM parameters
k=.7;
t=size(data,1);
c=zeros(t,1);
for i=m+1:t
    sig_est = (1/m)*sum((data(i-m:i)-a(i-m:i)).^2);
    c(i) = min([0,.5*log(s) - .5 * log(sig_est) + .5*(1/sig_est - 1/s)*sum((data(i-m:i)-a(i-m:i)).^2)]);
end

% Alarm location
for k=1:size(data,1)
    if c(k)<h
        step=k;
        return;
    end
end
end


