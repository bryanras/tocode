function xpy=whup(numpts)

% This is just an initialization function.
xpy=zeros(numpts,2);
for ii=0:numpts-1
	xpy(ii+1,1:2) = 2*[ cos(ii*2*pi/numpts) , sin(ii*2*pi/numpts) ];
end

