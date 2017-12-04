function [b]=lubksb(a,n,np,indx,b)
%INTEGER n,np,indx(n)
%REAL a(np,np),b(n)
%INTEGER i,ii,j,ll
%REAL sum
ii=0;
for i=1:n
    ll=indx(i);
    sum=b(ll);
    b(ll)=b(i);
    if (ii>0)
        for j=ii:i-1
            sum=sum-a(i,j)*b(j);
        end
    elseif (sum>0)
        ii=i;
    end
    b(i)=sum;
end
for i=n:-1:1
    sum=b(i);
    for j=i+1:n
        sum=sum-a(i,j)*b(j);
    end
    b(i)=sum/a(i,i);
end
return