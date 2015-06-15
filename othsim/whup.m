function xx=whup(pts,alp)

% This is a file for generating the original x_bar 
% for the uncoupled system.
%
% Usage : xx=whup(pts,alp)
%

aa=sqrt(alp);

% Declare some memory.
xx=zeros(prod(pts),8);

% Just go through the loop.
pos=0;
for ii=1:pts(1)
	phi1 = (ii-1)*2*pi/pts(1);
	cs1 = aa(1)*cos(phi1); si1 = aa(1)*sin(phi1);

	for jj=1:pts(2)
		phi2 = (jj-1)*2*pi/pts(2);
		cs2 = aa(2)*cos(phi2); si2 = aa(2)*sin(phi2);

		for kk=1:pts(3)
			phi3 = (kk-1)*2*pi/pts(3);
			cs3 = aa(3)*cos(phi3); si3 = aa(3)*sin(phi3);

			pos=pos+1;

			% Here it is.
			xx(pos,:) = [cs1, si1, cs2, si2, cs3, si3, ...
							(si1+cs1+si2+cs2+si3+cs3)/6, ...
							(si1-cs1+si2-cs2+si3-cs3)/6];
							
		end
	end
end
	
