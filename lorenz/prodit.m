function prod = prodit(JJ,sz)

% This is a routine for checking the product of the 'JJ' Jacobian against
% the product of the a's and b's.
%
% To get the product with n X n blocks, give it 
%
% prod = prodit(JJ,2).


N=size(JJ,1)/sz;

prod=eye(sz);

% Go from bottom to top.
if 1 == 0
	for ii= N:-1:1

		% Get the i+1 section.
		sec=((ii-1)*sz+1):(ii*sz);
		ipsec = (mod(ii,N)*sz+1):(mod(ii,N)*sz+sz);

		% Calculate.
		prod = -JJ(sec,ipsec)\JJ(sec,sec)*prod;

	end

end

	
% This code does it in a slightly more stable way.
if 1 == 1
	PP=JJ;
	PP( ((N-1)*sz+1):(N*sz), (N*sz+1):((N+1)*sz) ) = ...
			JJ( ((N-1)*sz+1):(N*sz), 1:sz );
	PP( ((N-1)*sz+1):(N*sz), 1:sz ) = zeros(sz,sz);
	PP( (N*sz+1):((N+1)*sz), 1:sz ) = eye(sz);

	rhs = zeros((N+1)*sz,sz);
	rhs((N*sz+1):((N+1)*sz),1:sz) = eye(sz);

	prod = PP\rhs;

	prod = prod((N*sz+1):((N+1)*sz),:);

end

