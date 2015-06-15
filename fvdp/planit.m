function xy = planit(huge)

%
% Usage: xy = planit(huge)
%
% This is a function for computing the planar sections, given the output
% of a run.
%
% See the notes on 7/14/03, p.4 if you are as confused as I am.
%

[pts, scs] = size(huge);

% Declare xy in order to prevent excessive memory allocation.
xy = zeros(pts,scs*2/3);

% Each lambda.
for ii=1:round(scs/3)

	% Do it column-by column.
	xy(:,2*ii) = huge(:,3*ii); 
	xy(:,ii*2-1) = sqrt( huge(:,ii*3-2).^2 + huge(:,ii*3-1).^2 );

end
		

		

