function G=lcfunc(xxl,oln,lam,pars)

%
% Usage: G = lcfunc(xxl,oln,lam,pars)
%
% This is a local function using the box scheme.
%
% xxl contains xx([ii,ii+],:).
% oln contains nb([ii,ii+],:).
% 
% pars = [ww,alp,ofst]
%

% Get the tangent.
tn=xxl(2,:)-xxl(1,:);

% Average the normals.
oln=mean(oln,1);

% Orthonormalize and adjust the sign.
[Q,R]   = qr( [tn; oln(1:3); oln(4:6)]' );
nn(1,:) = sign(R(2,2))*Q(:,2)';
nn(2,:) = sign(R(3,3))*Q(:,3)';

% Now we have the new normals. Here is the function.
G = nn*funcy(mean(xxl,1),lam,pars);

