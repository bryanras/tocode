function [xx,nb]=reord(xx,nb,lam)

%
% Usage: [xx,nb]=reord(xx,nb,lam)
%
% This function reorders xx so that the point with the smallest vector
% field in norm is xx(1,:);
%
% nb receives the same ordering.
%

N = size(xx,1);
ev=1; nmy=norm(funcy(xx(1,:),lam));

% Just check the norms.  Not much to this.
for ii=2:N

	nmyt = norm(funcy(xx(ii,:),lam));
	if nmyt < nmy
		ev=ii;
		nmy=nmyt;
	end

end

% Well, that's it.
if ev > 1
	scn = [ev:N,1:ev-1];
	xx=xx(scn,:);
	nb=nb(scn,:);
end
