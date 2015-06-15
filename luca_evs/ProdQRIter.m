function [T,U]=ProdQRIter(A,tol,Imax)

% This is a QR iteration for the product A1*A2*...*Ap
% The result is triangular matrices Ti and unitary U such that
% A1*A2*...*Ap=U*T1*T2*...*Tp*U'

ns=size(A); n=ns(1); p=ns(3);
[T,U]=Prod2HessTri(A);

for k=n:-1:2 , iters=0;
%   while abs(T(k,k-1,1))>tol & iters< Imax ,
%new relative error check below 
%D=eye(n); for j=2:p, D=D*diag(diag(T(:,:,j))); end
%   while abs(T(k,k-1,1)*D(k-1,k-1)) > tol*(abs(T(k,k,1)*D(k,k))+... %continue  
%   abs(T(k-1,k-1,1)*D(k-1,k-1))) & iters < Imax ,
    k2=max(1,k-2); E=T(k2:k,k-1:k,p);
    for j=p-1:-1:2, E=T(k2:k,k2:k,j)*E; end, E=T(k-1:k,k2:k,1)*E; 
    while abs(T(k,k-1,1)*E(1,1)) > tol*(abs(T(k,k,1)*E(2,2)+T(k,k-1,1)*E(1,2))+... 
             abs(T(k-1,k-1,1)*E(1,1))) & iters < Imax ,
   [T(1:k,:,:),W]=ProdQR(T(1:k,:,:)); iters=iters+1; 
   disp(sprintf('%g  %g  abserr=%g',k,iters,abs(T(k,k-1,1)))) ,
   U(1:n,1:k)=U(1:n,1:k)*W;  end, end
