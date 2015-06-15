function xx=isolit(JJ,FF)


% This is an iterative solution for JJxx=FF.
% It uses simple iteration with a diagonal pre-conditioner.

% First, check if the |a|'s are bigger than the |b|'s.
[ptsx ptsy] = size(JJ);

% Permute the columns of JJ if necessary.
flip=0;
if abs(JJ(1,1))<abs(JJ(1,2))
	flip=1;
	JJ=JJ(:,[2:ptsx,1]);
end

A=diag(diag(JJ));			% Just what it looks like.

cut=100;					% Max iters.
tol=1e-8;					% Tolerance.

xx=FF; xxt=xx; ko=2*tol; jj=0;
while (ko>tol) & (jj<cut)

	jj=jj+1;		% Update out cut-off;
	
	% Iterate.
	tmp=A\(FF-JJ*xx);
	ko=norm(tmp,inf);
	xx=xx+tmp;

end

% Permute the elements of y if necessary.
if flip==1
	xx=xx([ptsx,1:ptsx-1]);
end
