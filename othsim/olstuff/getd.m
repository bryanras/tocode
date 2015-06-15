function D=getd(JJ,pts)

% This is just a function for getting the block-diagonal 
% part of the Jacobian.

% Allocate the memory.
[tp tmp]=size(JJ);
bksz = round(tp/pts(1));

% Just a bunch of zeros.
D=spalloc(tp,tp,25*prod(pts)*4);

% Simple as can be.
for ii=1:pts(1)

	scn=(ii-1)*bksz+1:ii*bksz;

	D(scn,scn) = JJ(scn,scn);

end

	
