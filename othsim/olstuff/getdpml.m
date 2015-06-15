function DL=getdpml(JJ,pts,sgn)

% This is just a function for getting the block-diagonal 
% D-+L.

% sgn indicates D+L or D-L. It should be 1 or -1.

% Allocate the memory.
[tp tmp]=size(JJ);
bksz = round(tp/pts(1));

% Just a bunch of zeros.
DL=spalloc(tp,tp,25*pts(2)*pts(3)*(pts(1)+1)*4);

% Simple as can be.
DL(tp-bksz+1:tp,1:bksz) = sgn*DL(tp-bksz+1:tp,1:bksz);
for ii=1:pts(1)

	scn=(ii-1)*bksz+1:ii*bksz;

	DL(scn,scn) = JJ(scn,scn);

end

	
