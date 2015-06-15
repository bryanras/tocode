
% Script file.  Don't ask.
[pts h1] = size(huge);
xx=huge(:,h1-1:h1);
[pts h2] = size(xx);

% Get the new jacobian.
JJ = njacob(xx,nb,lam,aa);

% Getting the xminuses.
xm = [xx(2:pts,:);xx(1,:)] - xx;
xh = ([xx(2:pts,:);xx(1,:)] + xx)/2;

% Useful matrices.
II=eye(2);
AA=[0,1;-1,0];

% Calculte all this mess.
for ii=1:pts
	
	% Get the norms.
	nm=1/norm(xm(ii,:));
	
	% Projection and vector field.
	pN = II-(nm^2)*(xm(ii,:)'*xm(ii,:));
	dPhi = deefee(xh(ii,:),lam,aa);
	Phi = funcy(xh(ii,:),lam,aa);
	
	% Here are the derivatives.
	dgri(ii,1) = ...
			nm*(-Phi'*AA*pN - 0.5*xm(ii,:)*AA*dPhi )*nb(ii,:)';

	dgrip(ii,1) = ...
			nm*( Phi'*AA*pN - 0.5*xm(ii,:)*AA*dPhi )*nb(mod(ii,pts)+1,:)';
	

	% These are some more things to keep track of.
	lefty(ii,1) = Phi'*AA*pN*(nb(mod(ii,pts)+1,:)'-nb(ii,:)');
	righty(ii,1) = 0.5*xm(ii,:)*AA*dPhi*(nb(mod(ii,pts)+1,:)'+nb(ii,:)');
end



	
