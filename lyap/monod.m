% this function gets the monodromy matrix
% The syntax is the following
% function M=monod(YY)
%
% Input: 
% YY: values of transition matrices at the ts(i+1), i=1:N,
%     so, YY is dimensioned YY(n,n,N)
% Output:
% M: monodromy matrix
function M=monod(YY)

% MS matrix and RHS.
n=size(YY,1); N=size(YY,3); A=sparse(n*(N+1),n*(N+1)); rhs=sparse(n*(N+1),n);

id=eye(n); rhs(1:n,1:n)=id; A(1:n,1:n)=id; 

% Big loop, load up A.  This is more stable numerically than a long product.
for i=1:N, 
	im1n=(i-1)*n; in=i*n; A(in+1:in+n,im1n+1:im1n+n)=YY(:,:,i);  
	A(in+1:in+n,in+1:in+n)=-id;
end;

rhs=A\rhs; M=rhs(n*N+1:n*N+n,:);

