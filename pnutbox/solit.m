function yy = solit(JJ,ff)

% This is a function for pivoting and back-substituting the 
% jacobian equation for a 1-torus in 2-D.

% Pivot and solve.
[N,Np]=size(ff);
yy=zeros(N,1);

% Reduce.
ii=1;
while abs(JJ(ii,ii)) < abs(JJ(N,ii))

	% Permute.
	JJ([ii,N],:) = JJ([N,ii],:);
	ff([ii,N]) = ff([N,ii]);

	% Reduce.
	cc = JJ(N,ii)/JJ(ii,ii); 
	JJ(N,:) = JJ(N,:)-JJ(ii,:)*cc;
	ff(N) = ff(N)-ff(ii)*cc;
	
	ii=ii+1;
end

% Take care of the last element.
ii=ii-1;
cc = JJ(N,ii)/JJ(ii,ii); 
JJ(N,:) = JJ(N,:)-JJ(ii,:)*cc;
ff(N) = ff(N)-ff(ii)*cc;

% Now back-substitute.
yy(N) = ff(N)/JJ(N,N);
yy(N-1) = (ff(N-1)-JJ(N-1,N)*yy(N))/JJ(N-1,N-1);
for ii=N-2:-1:1
	yy(ii) = (ff(ii) - JJ(ii,ii+1)*yy(ii+1) - JJ(ii,N)*yy(N))/JJ(ii,ii);
end

