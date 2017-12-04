function [f,t]= voltra(n,m,t0,h,g,ak)
NMAX=5;
% Solves a set of m linear Volterra equations of the second kind using the extended trapezoidal
% rule. On input, t0 is the starting point of the integration and n-1 is the number of steps
% of size h to be taken. g(k,t) is a user-supplied external function that returns gk(t), while
% ak(k,l,t,s) is another user-supplied external function that returns the (k; l) element
% of the matrix K(t; s). The solution is returned in f(1:m,1:n), with the corresponding
% abscissas in t(1:n).

t(1)=t0;
for  k=1:m %Initialize.
%     f(k,1)=g(k,t(1));
    f(k,1)=g(1); % +1 to offset for 1 based indexes since t0==0
end
for i=2:n %Take a step h.
    t(i)=t(i-1)+1;
    for k=1:m
%         sum=g(k,t(i)); 
        sum=g(t(i)); %Accumulate right-hand side of linear equations in sum
        for l=1:m 
%             sum=sum+0.5*h*ak(k,l,t(i),t(1))*f(l,1);
            sum=sum+0.5*h*ak(1)*f(l,1);
            for j=2:i-1
%                 sum=sum+h*ak(k,l,t(i),t(j))*f(l,j);
                sum=sum+h*ak(i-j)*f(l,j);
            end
            if(k==l)%Left-hand side goes in matrix a.
                a(k,l)=1;
            else
                a(k,l)=0;
            end
            a(k,l)=a(k,l)-0.5*h*ak(i);
        end
        b(k)=sum;
    end
    [a,indx,d]=ludcmp(a,m,NMAX); %Decompose the matrix.
    
    b=lubksb(a,m,NMAX,indx,b); %Backsubstitute.
    for k=1:m
        f(k,i)=b(k);
    end
end
return