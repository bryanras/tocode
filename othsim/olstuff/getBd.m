function BB=getBd(JJ,pts)

% This is just a function for getting the off-block-diagonal 
% part of the Jacobian AND PUTTING IT ON THE DIAGONAL.

% Allocate the memory.
[tp tmp]=size(JJ);
bksz = round(tp/pts(1));
cdm=bksz/(pts(2)*pts(3));

% Just a bunch of zeros.
BB=spalloc(tp,tp,8*prod(pts)*cdm^2);

%  Do the last row (and hence the first column) first.
scnA = tp-bksz+1:tp; scnB=1:bksz;
BB(scnA,scnA)=JJ(scnA,scnB);

% Now do the rest. This looks bad, but it saves memory swapping.
for ii=1:pts(1)-1

	scnA=scnB; scnB=scnB+bksz;

	BB(scnA,scnA) = JJ(scnA,scnB);

end
	
