function M = rspcon(JJ,pts)

% This is a function for generating a row-section-wise preconditioner.

% Get the numbers.
[tp, jnk] = size(JJ);
cdm = round(tp/prod(pts));

% Declare the memory.
M=spalloc(tp,tp,tp*cdm);

scn=-cdm+1:0;					% Initialize the section.

% Go through all the boxes column-wise.
% Note that there is no reason to think that M will be invertible.
for ii=1:prod(pts)

	scn=scn+cdm;						% Get the section.
	
	% Trust in the first element.
	[y, ind] = max(abs(JJ(:,scn(1))));
	M( ind:ind+cdm-1 , scn ) = JJ( ind:ind+cdm-1 , scn );

end

% Invert that bad boy.
M=inv(M);

