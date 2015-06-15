
% Fill in a matrix based on C's and D's.
ibk=5*pts(2)*pts(3);

M=sparse(ibk*(pts(1)-1),ibk*(pts(1)-1),2*nnz(dd)+nnz(cc));
secp=1:ibk;

for ii=1:pts(1)-2

	ii
	sec=secp;
	secp=secp+ibk;

	M(sec,sec)=cc(sec,:);
	M(secp,sec)=dd(sec,:);
	
end


M(secp,secp)=cc(secp,:);
