function xx=spkill(xx,tolly)

% This function blanks out entries that have absolute
% value less than tolly.

ugly = size(xx);

for ii=1:ugly(1)

	for jj=1:ugly(2)

		if abs(xx(ii,jj))<tolly
			xx(ii,jj)=0;

		end
	end
end
