function [T,U,PT]=AlgL1(A)

% The result is triangular matrices Ti and unitary U such that
% A1*A2*...*Ap=U*T1*T2*...*Tp*U'

ns=size(A); n=ns(1); p=ns(3);
[T,U]=Prod2HessTri(A);

%form product into upper Hessenberg form
H=eye(n);
for j=p:-1:1, H=T(:,:,j)*H; end;
[QH,PT]=schur(H); U=U*QH;