function xx=dsub(xx,JJ,pts)

% This is just a function for dividing by the block 
% diagonal part of the Jacobian.

% Get the necessary stuff.
[tp tmp] = size(JJ);
bksz=round(tp/pts(1));

% Just substitute along the diagonal blocks.
for ii=1:pts(1)

	scn=(ii-1)*bksz+1:ii*bksz;		% Set the section.
	xx(scn) = JJ(scn,scn)\xx(scn);

end
	
