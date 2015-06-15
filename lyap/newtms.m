% This function performs the Newton loop for the
% multiple shooting solution of a periodic DE
% The syntax is the following
% function [s,tau,YY,flag]=newtms(s,tau,N,h,v0,f,Df,tol,maxit)
%
% Input: 
% s: initial guesses of solution at the points i/N 
% tau: initial guess of period 
% N: number of equispaced multiple shooting points
% h: integration stepsize to be used to go from point i to point (i+1)
%    Careful: h must be compatible with N.  That is, an integer multiple
%             of h must be equal to 1/N  
% v0: reference vector for phase condition
% f: string specifying the DE we want to solve
% Df: string specifying the Jacobian of the DE to solve
% tol: tolerance value for convergence of Newton's method
% maxit: max number of Newton iterations before giving up
% Output:
% s: values of solution at the points i/N
% tau: value of the period
% YY: values of transition matrices at the points i/N, i=1:N,
%     YY will be dimensioned YY(n,n,N)
% flag: success/insuccess indicator
%       flag='fail', failed to converge within Maxit iterations
%       otherwise flag is a vector [#Its,why] specifying first the 
%       number of iterations needed for convergence, then the 
%       reason the method converged: if norm(rhs)<tol --> why=1
%                                  if norm(update)<tol --> why=2
%       
function [s,tau,YY,flag]=newtms(s,tau,N,h,v0,f,Df,tol,maxit)

n=size(s(:,1),1);  k=1; % k is Newton's loop counter

while k<=maxit,

	[upd,YY,rhs]=musstep(tau,s,N,h,v0,f,Df);
	tau=tau+upd(n*N+1);

	for i=1:N, im1n=(i-1)*n; s(:,i)=s(:,i)+upd(im1n+1:im1n+n); end; 

	if (norm(rhs)<=tol), flag=[k;1]; return; end;

	if(norm(upd)<=tol), flag=[k;2]; return; end; 

	k=k+1;

end

flag='fail';
