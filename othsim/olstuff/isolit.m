function xx = isolit(JJ,ff,pts)

% This is an iterative function for solving a periodic, 
% block bi-diagonal system, JJyy = ff.

% Set the tolerance and maximum # of iterations.
cut=5; tolly=1e-4; 

% Get the numbers of points and all.
[tp tmp] = size(JJ);
bksz=round(tp/pts(1));

% Get these only if we are not doing block stuff.
L = tril(JJ,-1);
U = tril(JJ,1);
D = diag(diag(JJ));

% Declare some memory. 
xx=zeros(size(ff)); xxt=xx; xxt(1)=xxt(1)+2*tolly;

% As long as we don't get carried away, keep iterating.
ci=0;
while (ci<cut) & (norm(xxt-xx,inf)>tolly)


	xxt = xx;				% Save the old to check for convergence.

	% Block Gauss-Seidel.
%	xx = umult(xx,JJ,pts)+ff;
%	xx = dmlsub(xx,JJ,pts);

	% Block Jacobi.
	xx = uplmult(xx,JJ,pts)+ff;
	xx = Asub(xx,JJ,pts,0);

	% Gauss-Seidel.
% 	xx = U*xx+ff;
%	xx= (D-L)\xx;

	% Jacobi.
%	xx=(L+U)*xx+ff;
%	xx=D\xx;

	ci=ci+1;

	norm(xxt-xx,inf)

end


