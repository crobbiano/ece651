function [w] = kermom(y,x,m)
% Returns in w(1:m) the first m indefinite-integral moments of one row of the singular part
% of the kernel. (For this example, m is hard-wired to be 4.) The input variable y labels the
% column, while x (in COMMON) is the row.

% w(m),y,x,d,df,clog,x2,x3,x4

% We can take x as the lower limit of integration. Thus, we return the moment integrals either
% purely to the left or purely to the right of the diagonal.
if (y>=x)
    d=y-x;
    df=2*sqrt(d)*d;
    w(1)=df/3;
    w(2)=df*(x/3+d/5);
    w(3)=df*((x/3 + 0.4*d)*x + d^2/7);
    w(4)=df*(((x/3 + 0.6*d)*x + 3*d^2/7)*x + d^3/9);
else
    x2=x^2;
    x3=x2*x;
    x4=x2*x2;
    d=x-y;
    clog=log(d);
    w(1)=d*(clog-1);
    w(2)=-0.25*(3*x+y-2*clog*(x+y))*d;
    w(3)=(-11*x3+y*(6*x2+y*(3*x+2*y)) +6*clog*(x3-y^3))/18;
    w(4)=(-25*x4+y*(12*x3+y*(6*x2+y*(4*x+3*y)))+12*clog*(x4-y^4))/48;
end
return
