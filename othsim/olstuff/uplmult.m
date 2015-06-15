function xx=uplmult(xx,JJ,pts)

% This is just a function for multiplying by U+L, where
% U and L are the strictly upper- and lower-
% block-diagonal parts of the Jacobian respectively.

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
xx(scnA,:) = JJ(scnA,1:bksz)*xx(1:bksz,:);


