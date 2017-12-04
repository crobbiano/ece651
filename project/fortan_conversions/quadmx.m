function [a]=quadmx(n,np)
%INTEGER n,np,NMAX
%REAL a(np,np),PI
%DOUBLE PRECISION xx
%PARAMETER (PI=3.14159265,NMAX=257)
NMAX=257
%COMMON /momcom/ xx
%EXTERNAL kermom
%C USES wwghts,kermom
%Constructs in a(1:n,1:n) the quadrature matrix for an example Fredholm equation of the
%second kind. The nonsingular part of the kernel is computed within this routine, while the
%quadrature weights which integrate the singular part of the kernel are obtained via calls
%to wwghts. An external routine kermom, which supplies indefinite-integral moments of the
%singular part of the kernel, is passed to wwghts.
%INTEGER j,k
%REAL h,wt(NMAX),x,cx,y
h=2*pi/(n-1);
for j=1:n
    x=(j-1)*h;
    %xx=x; %Put x in COMMON for use by kermom.
    wt=wwghts(n,h,x);
    cx=cos(x); %Part of nonsingular kernel.
    for k=1:n
        y=(k-1)*h;
        a(j,k)=wt(k)*cx*cos(y); %Put together all the pieces of the kernel.
    end % 11
    a(j,j)=a(j,j)+1; %Since equation of the second kind, there is diagonal
end %12 piece independent of h.
return
