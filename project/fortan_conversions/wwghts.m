function [wghts]=wwghts(n,h,x)
% n
% wghts(n),h
% Constructs in wghts(1:n) weights for the n-point equal-interval quadrature from 0 to
% (n?1)h of a function f(x) times an arbitrary (possibly singular) weight function w(x) whose
% indefinite-integral moments Fn(y) are provided by the user-supplied subroutine kermom.
% INTEGER j,k
% DOUBLE PRECISION wold(4),wnew(4),w(4),hh,hi,c,fac,a,b
hh=h; %Double precision on internal calculations even though
hi=1/hh; %the interface is in single precision.
for j=1:n %Zero all the weights so we can sum into them.
    wghts(j)=0;
end
wold=kermom(0,x,4); %Evaluate indefinite integrals at lower end.
if (n>=4) %then Use highest available order.
    b=0; %For another problem, you might change this lower
    for j=1:n-3 %limit.
        c=j-1; %This is called k in equation (18.3.5).
        a=b; %Set upper and lower limits for this step.
        b=a+hh;
        if (j==n-3)
            b=(n-1)*hh; %Last interval: go all the way to end.
        end
        wnew=kermom(b,x,4);
        fac=1;
        for k=1:4 %Equation (18.3.4).
            w(k)=(wnew(k)-wold(k))*fac;
            fac=fac*hi;
        end
        %Equation (18.3.5).
        wghts(j)=wghts(j)+  ((c+1)*(c+2)*(c+3)*w(1) -(11+c*(12+c*3))*w(2) +3*(c+2)*w(3)-w(4))/6;
        wghts(j+1)=wghts(j+1)+ (-c*(c+2)*(c+3)*w(1) +(6+c*(10+c*3))*w(2) -(3*c+5)*w(3)+w(4))*.5;
        wghts(j+2)=wghts(j+2)+ (c*(c+1)*(c+3)*w(1) -(3+c*(8+c*3))*w(2) +(3*c+4)*w(3)-w(4))*.5;
        wghts(j+3)=wghts(j+3)+ (-c*(c+1)*(c+2)*w(1) +(2+c*(6+c*3))*w(2) -3*(c+1)*w(3)+w(4))/6;
        for k=1:4 %Reset lower limits for moments.
            wold(k)=wnew(k);
        end
    end
elseif (n==3) %Lower-order cases; not recommended.
    wnew=kermom(hh+hh,x,3);
    w(1)=wnew(1)-wold(1);
    w(2)=hi*(wnew(2)-wold(2));
    w(3)=hi^2*(wnew(3)-wold(3));
    wghts(1)=w(1)-1.5*w(2)+0.5*w(3);
    wghts(2)=2*w(2)-w(3);
    wghts(3)=0.5*(w(3)-w(2));
elseif (n==2)
    wnew=kermom(hh,x,2);
    wghts(2)=hi*(wnew(2)-wold(2));
    wghts(1)=wnew(1)-wold(1)-wghts(2);
end
return