function A=cknorm(B,dim)

% This is a file for checking the normality of a space O'vectors.
%
% Usage: A=cknorm(B,dim);
%
%
% Function returns progressive projections.  See 8/4/03 p. 1.

[tp,tot]=size(B);

vcs = tot/dim;			% Total number of vectors.

A = zeros(tp,vcs-1);		% Declare memory.

scn = reshape(1:tot,dim,vcs)';		% Sections.

% Go through the points.
for ii = 1:tp

	% Go through the vectors.
	for jj=2:vcs

		% Go down the daisy chain.
		for kk=jj-1:-1:1
			
			A(ii,jj-1) = A(ii,jj-1) + abs( B(ii,scn(jj,:))*B(ii,scn(kk,:))' );

			if jj==8 & ii<10 & (abs(B(ii,scn(jj,:))*B(ii,scn(kk,:))')>0.001)
				disp(sprintf('%d',kk))
			end

		end
	end
end
