function [a,b]=Givens(x)

% a,b are computed such that G=[a',b';-b,a] is unitary satisfying (G*x)(2)=0

%m=norm(x(1:2)); a=x(1)/m; b=x(2)/m;

% new Givens' computation (rescaled of above)
if abs(x(2))>abs(x(1)), x=x/abs(x(2)); else
x=x/abs(x(1)); end
if x(2)==0, a=1; b=0; else
m=norm(x(1:2)); a=x(1)/m; b=x(2)/m; end