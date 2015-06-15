function [ll,u1,u2,pp] = luit(JJ)

% This is a function for pivoting and decomposing the 
% jacobian equation for a 1-torus in 2-D.

% Pivot and solve.
[N,Np]=size(JJ);
pp=[1:N]'; ll=zeros(N-1,1); u1=zeros(N,1); u2=zeros(N-1,1);

% Reduce.
cc=JJ(N,1); dd=JJ(N,N);
for ii=1:N-1

	% Permute, if necessary.
	if abs(JJ(ii,ii)) < abs(cc)
		pp([ii, N]) = pp([N, ii]);
		ll(ii) = JJ(ii,ii)/cc;
		u1(ii)=cc; u2(ii)=dd;
		cc=JJ(ii,ii+1); dd=JJ(ii,N)-ll(ii)*dd;

	else
		ll(ii)=cc/JJ(ii,ii);
		u1(ii)=JJ(ii,ii); u2(ii) = JJ(ii,ii+1);
		cc=-JJ(ii,ii+1)*ll(ii);  dd=dd-JJ(ii,N)*ll(ii);
	end

end

% Take care of the last element.
u1(N) = dd;

