function G=lcfunc(xxl,oln,lam,alp,eps0)

% This is a local function using the box scheme.

% The rows of xxl, in order, are  ijk, i+jk, i+j+k, ij+k, ijk+, i+jk+,
% i+j+k+, i,j+,k+

% Get the tangents.
tti=sum(xxl([2,3,6,7],:)-xxl([1,4,5,8],:));
ttj=sum(xxl([3,4,7,8],:)-xxl([2,1,6,5],:));
ttk=sum(xxl([5,6,7,8],:)-xxl([1,2,3,4],:));

% Orthonormalize.
[Q,R]=qr([tti;ttj;ttk;oln]');
oln(1,:) = sign(R(4,4))*Q(:,4)';
oln(2,:) = sign(R(5,5))*Q(:,5)';
oln(3,:) = sign(R(6,6))*Q(:,6)';
oln(4,:) = sign(R(7,7))*Q(:,7)';
oln(5,:) = sign(R(8,8))*Q(:,8)';


% Here is the function.
G = oln*funcy(mean(xxl,1),lam,alp,eps0);

