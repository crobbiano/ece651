function [i]=Shewhart(data,a,s)
% Turn counts into count rate
data=data./100;
a=a./100;
s=s/100;

% Shewhart parameters
d=3;

% Alarm location
for i=1:size(data,1)
    if data(i)>a(i)+d*s
        step=i;
        return;
    end
end
end