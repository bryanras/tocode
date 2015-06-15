function [T,U]=ProdQR(S)

% One QR-iteration step is performed to the product S1*S2*...*Sp
% matrices T2,...,Tp and S2,...,Sp are upper triangular 
% and T1, S1 are Hessenberg. 

T=S; ns=size(S); m=ns(1); n=ns(2); p=ns(3); U=eye(m); Sp=S(:,:,p);
m2=max(1,m-2); E=T(m2:m,m-1:m,p);
for j=p-1:-1:2, E=T(m2:m,m2:m,j)*E; end
E=T(m-1:m,m2:m,1)*E; 
%ev=eig(E); mu=ev(1); if abs(mu)>abs(ev(2)) , mu=ev(2); end, mc=mu;
%disp(sprintf('old shift=%g',mc)) ,
%% Wilkinson shift to make it closest eigenvalue of E to E(2,2) %%
ev=eig(E); mu=ev(1); 
if abs(mu-E(2,2))>abs(ev(2)-E(2,2)), mu=ev(2); end, mc=mu;
%disp(sprintf('new shift=%g',mc)) ,
% last element shift (Rayleigh quotient shift)
%mu=E(2,2); mc=mu;
% no shift
%mu=0; mc=mu;

for k=1:m-1 , k1=max(1,k-1); k2=min(m,k+2);
    w=Sp(k1:k,k);
    for j=p-1:-1:2, w=T(k1:k,k1:k,j)*w; end, w=T(k:k+1,k1:k,1)*w;
    [a,b]=Givens([w(1)-mc,w(2)]); 
    W=[a',b';-b,a]; mc=a*mu;
    T(k:k+1,1:n,1)=W*T(k:k+1,1:n,1); 
    U(1:k+1,k:k+1)=U(1:k+1,k:k+1)*W';
    T(1:k+1,k:k+1,p)=T(1:k+1,k:k+1,p)*W'; 
    for j=p:-1:2
       [a,b]=Givens(T(k:k+1,k,j)); W=[a',b';-b,a]; 
       if j==p , Sp(k:k+1,k:n)=W*Sp(k:k+1,k:n); end
       T(k:k+1,k:n,j)=W*T(k:k+1,k:n,j); 
       T(1:k2,k:k+1,j-1)=T(1:k2,k:k+1,j-1)*W';end, 
end