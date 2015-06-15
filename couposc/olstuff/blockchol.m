function [cc,dd]=blockchol(cc,dd)

% This does a block cholesky on a block-tri-diagonal system.
% cc is the diagonal, while dd is the upper-off-diagonal.

% There may be a bug or two in this.  See G & V.L., p. 146.

[tp,ibk]=size(cc);
bks=tp/ibk;

% First column of blocks.
sec=1:ibk;
cc(sec,:)=chol(cc(sec,:))';				% Diagonal.
dd(sec,:)=dd(sec,:)'/cc(sec,:)' ;		% Off-diagonal. 

% Note that we are filling in the transpose of the diagonal. 


% Now get the others.
for jj=2:bks

	sec=sec+ibk;

	% Diagonal.
	cc(sec,:)=cc(sec,:)-dd(sec-ibk,:)*dd(sec-ibk,:)';	
	cc(sec,:)=chol(cc(sec,:))';

	% Off-diagonal.
	dd(sec,:)=dd(sec,:)'/cc(sec,:)';

end
	

