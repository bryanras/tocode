function [h, tau, A, B, M] = colloc(xx,lam)

%
% Usage: [h, tau, A, B, M] = colloc(xx,lam)
%
%
% This is a script for generating all the matrices and values associated
% with collocation, given an orthogonality solution.
%
% See the notes on 17 Sept. for details.
%

pts=size(xx,1);
dm = size(xx,2);
II = eye(dm);
tau = 0;

% Set aside some memory.
fham=spalloc(dm*(pts+1), dm*(pts+1), (dm^2)*2*(pts+1));


% Run through everything.
for ii=1:pts

	ip = mod(ii,pts)+1;

	% Get the half-point and tangent.
	xtan=xx(ip,:)-xx(ii,:);
	xh = (xx(ip,:)+xx(ii,:))/2;

	% Get the stuff associated with the vector field.
	Phih = funcy(xh,lam);
	dPhih = dfunc(xh,lam);

	% Now go for it.
	h(ii)=norm(xtan)/norm(Phih);  % Make this minus if reversing field.
	B(:,:,ii) = -(II + h(ii)*dPhih/2);
	A(:,:,ii) =  (II - h(ii)*dPhih/2);

	% This will be necessary when solving for the monodromy.
	fham( (ii*dm-1):(ii*dm), (ii*dm-1):(ii*dm) ) = A(:,:,ii);
	fham( (ii*dm-1):(ii*dm), (ii*dm+1):(ii*dm+2) ) = B(:,:,ii);

	tau=tau+h(ii);

end

% Now solve for the monodromy in a numerically stable way.
fham(dm*pts+1:dm*(pts+1),1:dm) = II;
rhs = zeros((pts+1)*dm,dm);
rhs(1:dm,1:dm) = II;

M = fham\rhs;
M = M(pts*dm+1:(pts+1)*dm,:);

