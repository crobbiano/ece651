function [a,indx,d]= ludcmp(a, n, np)
%PARAMETER (NMAX=100,TINY=1.0E-20)
%DIMENSION A(NP,NP),INDX(N),VV(NMAX)
d=1;
tiny=1.0e-20;
for i=1:n
    aamax=0;
    for j=1:n
        if (abs(a(i,j))>aamax)
            aamax=abs(a(i,j));
        end
    end
    %if (aamax == 0) PAUSE 'Singular matrix.'
    vv(i)=1.0/aamax;
end
for j=1:n
   % if (j>1)
        for i=1:j-1
            sum=a(i,j);
            %if (i>1)
                for k=1:i-1
                    sum=sum-a(i,k)*a(k,j);
                end
                a(i,j)=sum;
            %end
        end
   % end
    aamax=0;
    for i=j:n
        sum=a(i,j);
        %if (j>1)
            for k=1:j-1
                sum=sum-a(i,k)*a(k,j);
            end
            a(i,j)=sum;
        %end
        dum=vv(i)*abs(sum);
        if (dum>=aamax)
            imax=i;
            aamax=dum;
        end
    end
    if (j~=imax)
        for k=1:n
            dum=a(imax,k);
            a(imax,k)=a(j,k);
            a(j,k)=dum;
        end
        d=-d;
        vv(imax)=vv(j);
    end
    indx(j)=imax;
    if(a(j,j)==0)
        a(j,j)=tiny;
    end
    if(j~=n)
        dum=1.0/a(j,j);
        for i=j+1:n
            a(i,j)=a(i,j)*dum;
        end
    end
end
if(a(n,n)==0)
    a(n,n)=tiny;
end
return
