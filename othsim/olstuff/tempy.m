
ibk=5*pts(2)*pts(3);
multyi = eye(ibk);
multyii = multyi;

for ii=1:pts(1)

	scn=(ii-1)*ibk+1:ii*ibk;
	scnp=scn+ibk;
	if scnp(1)>pts(1)*ibk scnp=scnp-pts(1)*ibk; end

	multyi=(JJ(scn,scnp)\JJ(scn,scn))*multyi;
	multyii=(JJ(scn,scn)\JJ(scn,scnp))*multyii;

end
