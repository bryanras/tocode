function Dx = deefee(x,lam,aa)

tre = (x(1)^2+x(2)^2)^(3/2);
cinq = (x(1)^2+x(2)^2)^(5/2);

Dx(1,1) = lam - 3*x(1)^2 - x(2)^2 + 3*aa*x(1)^2/tre - 3*aa*x(1)^4/cinq;
Dx(1,2) = -2*x(1)*x(2) - 3*aa*x(1)^3*x(2)/cinq -1;
Dx(2,2) = -2*x(1)*x(2) - +2*aa*x(1)*x(2)/tre -3*aa*x(1)^3*x(2)/cinq +1;
Dx(2,2) = lam - x(1)^2 - 3*x(2)^2 + aa*x(1)^2/tre - 3*aa*(x(1)^2)*(x(2)^2)/cinq;
