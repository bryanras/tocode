function [xx, enns] = whup(pts)

% 
% Usage: [xx, enns] = whup(pts)
%
% This is a file for generating the original x_bar 
% for the uncoupled system.

% Declare some memory.
xx=zeros(prod(pts),4);
enns=zeros(prod(pts),8);

pos=0;
for ii=1:pts(1)
	phi1 = (ii-1)*2*pi/pts(1);
	cs1 = cos(phi1); si1 = sin(phi1);

	for jj=1:pts(2)
		phi2 = (jj-1)*2*pi/pts(2);
		cs2 = cos(phi2); si2 = sin(phi2);

		% There is the initial guess with normal vector.
		pos=pos+1;
		xx(pos,:) = [cs1, si1, cs2, si2];
		enns(pos,:) = [cs1, si1, 0, 0, 0, 0, cs2, si2];
							
	end
end
	
