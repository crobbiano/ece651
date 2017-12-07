function [k,c]=CUSUM(data,a,s,h)
% CUSUM parameters
k=.7;
t=size(data,1);
c=zeros(t,1);
ndata=(data-a)./s; %normalize data
for j=1:t
    if j==1
        c(j)=max(0,(ndata(j)-k)); % c(0)=0
%         c(j)=((ndata(j))); % c(0)=0
    else
        c(j)=max(0,(ndata(j)-k+c(j-1)));
%         c(j)=((ndata(j))+c(j-1));
    end
end

% Alarm location
for k=1:size(data,1)
    if c(k)>h
        step=k;
        return;
    end
end
end
