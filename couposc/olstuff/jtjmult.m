function yy=jtjmult(xx,cc,dd)

% This function multiplies by J^TJ, given the diagonals and off-diagonals
% of J^TJ in cc and dd respectively.
[tp,ibk]=size(cc);
bks=tp/ibk;					% Should be N_i-1.

% Just loop down. The ends are a little different.
yy=zeros(size(xx));

secm=tp-ibk+1:tp; sec=1:ibk; secp=sec+ibk;
yy(sec,:) = cc(sec,:)*xx(sec,:) + dd(sec,:)*xx(secp,:);

for ii=1:bks

	% Multiply.
	yy(sec,:) = dd(secm,:)'*xx(secm,:)+ cc(sec,:)*xx(sec,:) + ...
					dd(sec,:)*xx(secp,:);

	% Update the sections.
	secm=sec; sec=secp; secp=secp+ibk;

	% Rather than a mod function, ...
	if secp(1)>tp secp=1:ibk; end

end

