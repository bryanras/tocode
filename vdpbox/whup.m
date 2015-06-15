function xpy=whup(numpts)

% This is just an initialization function.
xpy=zeros(numpts,2);
for ii=0:numpts-1
	theta = -2*pi*ii/numpts;
	xpy(ii+1,1:2) = [ 2*cos(theta) , 2*sin(theta) ];
end

