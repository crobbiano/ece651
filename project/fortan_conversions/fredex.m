function [g,x]=fredex(n)
% INTEGER NMAX
% REAL PI
% PARAMETER (NMAX=100,PI=3.14159265)
% INTEGER indx(NMAX),j,n
% REAL a(NMAX,NMAX),g(NMAX),x,d
% C USES quadmx,ludcmp,lubksb
NMAX=100;
%This sample program shows how to solve a Fredholm equation of the second kind using
%the product Nystrom method and a quadrature rule especially constructed for a particular,
%singular, kernel.
%n=40; %Here the size of the grid is specified.
a=quadmx(n,NMAX); %Make the quadrature matrix; all the action is here.

[a,indx,d]=ludcmp(a,n,NMAX); %Decompose the matrix.
for j=1:n %Construct the right hand side, here sin x.
    x=(j-1)*2*pi/(n-1);
    g(j)=sin(x);
end
g=lubksb(a,n,NMAX,indx,g); %Backsubstitute.
for j=1:n %Write out the solution.
    x(j)=(j-1)*pi/(n-1);
    g(j);
    %write (*,*) j,x,g(j)
end
%write (*,*) ’normal completion’
return