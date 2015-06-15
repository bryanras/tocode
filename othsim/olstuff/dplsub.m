function xx=dplsub(xx,JJ,pts)

% This is just a function for dividing by D+L, where
% D is the block-diagonal part of the Jacobian, and
% L is just the lower left corner.

% Get the necessary stuff.
[tp tmp] = size(JJ);
bksz=round(tp/pts(1));

% Just substitute along the diagonal blocks.
scn=1:bksz;
for ii=1:pts(1)-1

	xx(scn,:) = JJ(scn,scn)\xx(scn,:);
	scn=scn+bksz;	% Set the section.

end
	
% Take care of the last row block.
xx(scn,:) = JJ(scn,scn)\(xx(scn,:)-JJ(scn,1:bksz)*xx(1:bksz,:));


