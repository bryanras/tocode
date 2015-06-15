function xx=dioffsub(xx,JJ,pts)

% This is just a function for dividing by the 
% block-upper-triangular part of the Jacobian.

% Get the necessary stuff.
[tp tmp] = size(JJ);
bksz=round(tp/pts(1));

% Just substitute along the diagonal blocks.
for ii=1:pts(1)-1

	scnA=(ii-1)*bksz+1:ii*bksz;		% Set the sections.
	scnB=scnA+bksz;
	
	xx(scnA) = JJ(scnA,scnA)\(xx(scnA)-JJ(scnA,scnB)*xx(scnB));

end
	
xx(scnB) = JJ(scnB,scnB)\xx(scnB);
