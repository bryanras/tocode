function xx=Asub(xx,JJ,pts,typ)

% This is just a function for multiplying or dividing by the block 
% diagonal part of the Jacobian.

% typ=0  => divide
% typ=1  => multiply

% Get the necessary stuff.
[tp tmp] = size(JJ);
bksz=round(tp/pts(1));

% Just substitute along the diagonal blocks.
scn=-bksz+1:0;
for ii=1:pts(1)

	scn=scn+bksz;	% Set the section.

	if typ == 0
		xx(scn,:) = JJ(scn,scn)\xx(scn,:);
	else
		xx(scn,:) = JJ(scn,scn)*xx(scn,:);
	end

end
	
