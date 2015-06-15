function outmat = ordrtst(N)

%
% Usage: outmat = ordrtst(N)
%
% This is just a function for testing an order argument that I have.
% See the notes on 9/25/03.

% All these are as in the notes.
b = 1-1/N;
a= -b+1/N;
c= a;
d = b;
h=1/N;

outmat=eye(2);
for ii=1:N

	outmat=outmat*[d, -h; h, b]\[c, -h; h, a];

end
	
