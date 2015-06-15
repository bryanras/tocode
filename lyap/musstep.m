% This function sets up the multiple shooting matrix
% arising during Newton's method for periodic DEs
% It further solves the system for the Newton's correction
% and returns the transition matrices at the shooting points
% The syntax is the following
% function [upd,YY,rhs]=musstep(tau,s,N,h,v0,f,Df)
%
% Input: 
% tau: value of period (current guess)
% s: current initial conditions at the points i/N (shooting guesses)
% N: number of equispaced multiple shooting points
% h: integration stepsize to be used to go from point i to point (i+1)
%    Careful: h must be compatible with N.  That is, an integer multiple
%             of h must be equal to 1/N  
% v0: reference vector for phase condition
% f: string specifying the DE we want to solve
% Df: string specifying the Jacobian of the DE to solve
% Output:
% upd: updates for the shooting guesses and update for tau
%      [first the shooting guesses, then tau)
% YY: values of transition matrices at the ts(i+1), i=1:N,
%     so, YY is dimensioned YY(n,n,N)
% rhs: right-hand-side of system to monitor if it is small
function [upd,YY,rhs]=musstep(tau,s,N,h,v0,f,Df)

% MS matrix and RHS
n=size(s(:,1),1);  M=sparse(n*N+1,n*N+1); rhs=zeros(n*N+1,1); 

% ICs for transitions
id=eye(n); YY=zeros(n,n,N); for i=1:N, YY(:,:,i)=id; end;

% ICs for derivative w.r.t. tau
v=zeros(n,N);

% ICs for solution (shooting guesses).
x=s(:,1:N);

m=1/(h*N); % Number of h-steps to cover one subinterval.

% Big loop.
for i=1:N

	% Inner integration loop.
	for k=1:m 

		x0=x(:,i); Y0=YY(:,:,i); vz=v(:,i);

		[x1,Y1,v1]=rk38(h,tau,x0,Y0,vz,f,Df); 

		x(:,i)=x1; YY(:,:,i)=Y1; v(:,i)=v1;
		
	end % obtained solutions at point (i+1)

	im1n=(i-1)*n;
	M(im1n+1:im1n+n,im1n+1:im1n+n)=YY(:,:,i);
	M(im1n+1:im1n+n,n*N+1)=v(:,i); 

	mm=mod(i,N);
	rhs(im1n+1:im1n+n)=-(x(:,i)-s(:,mm+1));
	mm1=mm*n;
	M(im1n+1:im1n+n,mm1+1:mm1+n)=-id;

end

ff=feval(f,v0); M(n*N+1,1:n)=ff'; rhs(n*N+1)=-ff'*(s(:,1)-v0);

upd=M\rhs;

