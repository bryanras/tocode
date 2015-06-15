function G=lcfunc(xxl,oln,alp,lam)

%
% Usage: G = lcfunc(xxl,oln,alp,lam)
%
% This is a local function using the box scheme.
%
% The ordered rows of xxl are  ij, i+j, ij+, i+j+.
% oln contains [n1(i,j); n2(i,j)]
% 

% Get the tangents.
tt1=xxl(4,:)-xxl(1,:);
tt2=xxl(2,:)-xxl(3,:);

% Orthonormalize and adjust the sign.
[Q,R]=qr([tt1;tt2;oln]');
oln(1,:) = sign(R(3,3))*Q(:,3)';
oln(2,:) = sign(R(4,4))*Q(:,4)';

% Now we have the new normals. Here is the function.
G = oln*funcy(mean(xxl,1),alp,lam);

