
% Another script file.  Again, don't ask.

% Calculating the a's and b's.

[rows, cols] = size(JJ);
aob=1;
boa=1;
ugh=1;
for ii=1:rows

	ip = mod(ii,rows)+1;

	% From the Jacobian.
	aob = aob*JJ(ii,ii)/JJ(ii,ip);
	boa = boa*JJ(ii,ip)/JJ(ii,ii);

	% From what I think the limit should be.
	xtan=xx(ip,:)-xx(ii,:);
	xh=mean(xx([ii,ip],:),1);
	nn=[xtan(2), -xtan(1)];nn=nn/norm(nn);

	DP=dPhi(xh,lam);
	aay = 0.5*nn*(DP*nb(ii,:)') - (nn*nb(ii,:)'/norm(xtan))* ...
			[nn(2), -nn(1)]*funcy(xh,lam);

	bee = 0.5*nn*(DP*nb(ip,:)') + (nn*nb(ip,:)'/norm(xtan))* ...
			[nn(2), -nn(1)]*funcy(xh,lam);

	if ii<5
		disp(sprintf('\t %g \t%g', JJ(ii,ii)-aay, JJ(ii,ip)-bee))
	end
	
	ugh=ugh*bee/aay;
end

aob 
boa
ugh
