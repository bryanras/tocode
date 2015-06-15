function [yy,npi,npj] = chopit( xx, ptsi, ptsj, ifac, jfac )

% Usage: [yy,npi,npj] = chopit( xx, ptsi, ptsj, ifac, jfac )
% 
% This is a function for reducing huge tori.
% ifac and jfac are integers that represent the sampling rate.
%
% For example, to get roughly every fifth point in the i direction and
% every 3rd point in the j direction, say
%
% [yy,npi,npj] = chopit( xx, ptsi, ptsj, 5, 3 )
%
% The function goes by rows, so the number of columns in xx is irrelevant.
%
% npi is of course "new ptsi". Ditto for npj.


npi=length(1:ifac:ptsi);

yy=[];
for ii=1:ifac:ptsi

	scn=(ii-1)*ptsj+1:jfac:ii*ptsj;

	yy=[yy; xx(scn,:)];

end

npj=length(scn);

