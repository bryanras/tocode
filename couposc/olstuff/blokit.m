function [ijk,ipjk,ipjpk,ipjkp,ipjpkp,ijpk,ijpkp,ijkp]=blokit(rbk,pts,JJ)


% This is a function for getting the 
% rbkth row block of the Jacobian.
% Usage:
% [ijk, ipjk, ipjpk, ipjkp, ipjpkp, ijpk, ijpkp, ijkp]=blokit(rbk,pts,JJ)

% Get the sections. This isn't well-tested.
scn=(rbk-1)*5+1:rbk*5;

% i+ section:
if rbk/(pts(2)*pts(3)) > pts(1)-1
	ipa = -5*pts(3)*pts(2)*(pts(1)-1);
else
	ipa = 5*pts(3)*pts(2);
end

% j+ section:
if mod(rbk,pts(3)*pts(2))/pts(3) > pts(2)
	jpa = -5*pts(3)*(pts(2)-1);
else
	jpa =5*pts(3);
end

% k+ section:
if mod(rbk,pts(3)) == 0
	kpa = -5*(pts(3)-1);
else
	kpa = 5;
end

ijk		=	full(JJ(scn,scn));
ipjk	=	full(JJ(scn,scn+ipa));
ipjpk	=	full(JJ(scn,scn+ipa+jpa));
ipjpkp	=	full(JJ(scn,scn+ipa+jpa+kpa));
ipjkp	=	full(JJ(scn,scn+ipa+kpa));
ijpkp	=	full(JJ(scn,scn+jpa+kpa));
ijpk	=	full(JJ(scn,scn+jpa));
ijkp	=	full(JJ(scn,scn+kpa));

