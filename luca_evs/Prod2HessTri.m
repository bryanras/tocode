function [T,U]=Prod2HessTri(A)

% Unitary transformations are applied such that 
% A1*A2*...*Ap=U*T1*T2*...*Tp*U'
% where matrices T2,...,Tp are upper triangular and T1 is Hessenberg

T=A; ns=size(A); n=ns(1); p=ns(3); U=eye(n);
for k=1:n-1 ,
for j=p:-1:2
   w=HouseVect(T(k:n,k,j)); T(k:n,k:n,j)=T(k:n,k:n,j)-w*(w'*T(k:n,k:n,j));
   T(:,k:n,j-1)=T(:,k:n,j-1)-(T(:,k:n,j-1)*w)*w'; end
   if k<n-1 ,
     w=HouseVect(T(k+1:n,k,1)); 
     T(k+1:n,k:n,1)=T(k+1:n,k:n,1)-w*(w'*T(k+1:n,k:n,1));
     U(:,k+1:n)=U(:,k+1:n)-(U(:,k+1:n)*w)*w';
     T(:,k+1:n,p)=T(:,k+1:n,p)-(T(:,k+1:n,p)*w)*w'; end, end
% force the correct structure, just in case
for k=2:p, T(:,:,k)=triu(T(:,:,k),0); end; T(:,:,1)=triu(T(:,:,1),-1);
