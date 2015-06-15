function [T,U]=Prod2TriHess(A)

% Unitary transformations are applied such that 
% A1*A2*...*Ap=U*T1*T2*...*Tp*U'
% where matrices T1,...,T_(p-1) are upper triangular and Tp is Hessenberg
% 
% this choice makes it easier to choose shifts

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
for k=1:n-1, 
[a,b]=Givens(T(k:k+1,k,1)); W=[a',b';-b,a];
T(k:k+1,k:n,1)=W*T(k:k+1,k:n,1);
U(1:n,k:k+1)=U(1:n,k:k+1)*W'; T(1:k+1,k:k+1,p)=T(1:k+1,k:k+1,p)*W';
end