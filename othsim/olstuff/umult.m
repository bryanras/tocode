function xx=umult(xx,JJ,pts)

% This is just a function for multiplying by U, where
% U is the strictly upper-block-diagonal part of the Jacobian.

% Get the necessary stuff.
[tp tmp] = size(JJ);
bksz=round(tp/pts(1));

% Just substitute along the diagonal blocks.
scnA=1:bksz;
for ii=1:pts(1)-1
	
	scnB = scnA+bksz;

	xx(scnA,:) = JJ(scnA,scnB)*xx(scnB,:);

	scnA=scnB;	% Set the section.

end
	
% Take care of the last row block.
xx(scnA,:) = 0*xx(scnA,:);


