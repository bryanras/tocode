function xx=cdiv(xx,cc,dd)

% This is a pre-conditioning function for the J^TJ-solution.
% This preconditioner is simply the block diagonal.

% Block sizes, as before.
[tp, ibk]=size(cc); [tp, tmp]=size(xx); bks=tp/ibk;

for ii=1:bks

	scn=(ii-1)*ibk+1:ii*ibk;
	xx(scn,:)=cc(scn,:)\xx(scn,:);

end
	
	
