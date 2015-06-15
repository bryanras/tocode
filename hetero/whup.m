function xpy=whup(numpts)

% This is just an initialization function.
xpy=zeros(numpts,2);
for ii=0:numpts-1
	theta=ii*2*pi/numpts;
	xpy(ii+1,1:2) = [ cos(theta) , sin(theta) ];
end

