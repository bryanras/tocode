function [cc,dd]=csnds(JJ,pts)

% This function returns the diagonal and off-diagonal blocks
% of a block periodic bi-diagonal J.
[tp,tmp]=size(JJ);

% Get the stats.
cdm = tp/prod(pts);
ibk = pts(2)*pts(3)*cdm;

% cc and dd are as in the notes from 3/19/03.
for ii=1:pts(1)

	% Update the sections.
	secAi=(ii-1)*ibk+1:ii*ibk;

	tmp = mod(ii-2,pts(1));
	secBim=tmp*ibk+1:(tmp+1)*ibk;

	tmp = mod(ii,pts(1));
	secBi=tmp*ibk+1:(tmp+1)*ibk;

	% Add on to the D and C vectors.
	cc(secAi,:) = JJ(secAi,secAi)'*JJ(secAi,secAi)+ ...
						JJ(secBim,secAi)'*JJ(secBim,secAi);
	dd(secAi,:) = JJ(secAi,secAi)'*JJ(secAi,secBi);

end
