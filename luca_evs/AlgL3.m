function [T,U,PT]=AlgL3(A)

% The result is triangular matrices Ti and unitary U such that
% A1*A2*...*Ap=U*T1*T2*...*Tp*U'

ns=size(A); n=ns(1); p=ns(3);

[q(:,:,p),T(:,:,p)]=qr(A(:,:,p));
for j=p-1:-1:1, 
   [q(:,:,j),T(:,:,j)]=qr(A(:,:,j)*q(:,:,j+1));
end

%form product of triangulars and orthogonal
H=eye(n);
for j=p:-1:1, H=T(:,:,j)*H; end; H=q(:,:,1)*H;
[QH,PT]=schur(H); U=QH; 
