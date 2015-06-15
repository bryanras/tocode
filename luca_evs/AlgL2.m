function [T,U]=AlgL2(A,tol,Imax)

%this is alg=3 (2nd Luca's algorithm)

%The result is triangular PT1 and PT2 and unitary U 
%such that A1*A2*...*Ap=U*PT1*PT2*U'

ns=size(A); n=ns(1); p=ns(3);
[T,U]=Prod2HessTri(A);

%form product of T2, ..., Tp into unique triangular 
TT=eye(n);
for j=p:-1:2, TT=T(:,:,j)*TT; end;

%now can proceed as in the algorithm 1 (Timo's original) 
%pretending that p=2

p=2; T(:,:,2)=TT;
for k=n:-1:2 , iters=0;
    k2=max(1,k-2); E=T(k2:k,k-1:k,p); E=T(k-1:k,k2:k,1)*E; 
    while abs(T(k,k-1,1)*E(1,1)) > tol*(abs(T(k,k,1)*E(2,2)+T(k,k-1,1)*E(1,2))+... 
             abs(T(k-1,k-1,1)*E(1,1))) & iters < Imax ,
   [T(1:k,:,1:2),W]=ProdQR(T(1:k,:,1:2)); iters=iters+1; 
   disp(sprintf('%g  %g  abserr=%g',k,iters,abs(T(k,k-1,1)))) ,
   U(1:n,1:k)=U(1:n,1:k)*W;  end, end
