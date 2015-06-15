function UL=getupml(JJ,pts,sgn)

% This is just a function for getting U+-L, where L is the strictly 
% lower-block diagonal part of the Jacobian, and U is the strictly
% upper-block diagonal.

% sgn is +- 1, depending on U+L or U-L.

% Allocate the memory.
[tp tmp]=size(JJ);
bksz = round(tp/pts(1));
cdm=bksz/(pts(2)*pts(3));

% Just a bunch of zeros.
UL=spalloc(tp,tp,8*prod(pts)*cdm^2);

%  Do the last row (and hence the first column) first.
scnA = 1:bksz; scnB=scnA;
UL(tp-bksz+1:tp,scnA)=sgn*JJ(tp-bksz+1:tp,scnA);

% Now do the rest. This looks bad, but it saves memory swapping.
for ii=1:pts(1)-1

	scnA=scnB; scnB=scnB+bksz;

	UL(scnA,scnB) = JJ(scnA,scnB);

end
	
