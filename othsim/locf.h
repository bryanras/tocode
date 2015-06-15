#include "mex.h"
#include<math.h>

/*
 * This header file contains the vector field generator and the 
 * local function (used in locfunc.c and njacob).
 *
 * The rows of xx, in order, are  ijk, i+jk, i+j+k, ij+k, ijk+, i+jk+,
 * i+j+k+, i,j+,k+
 */


/* This function returns the vector field. */
/* It is basically the same thing as funcy.m */
double *fcy(double *xx, double lam, double *alp, double eps0)
{
	double *xpy;
	
	xpy=mxMalloc(sizeof(double)*8);

	/* This is everything. */
	xpy[0] = alp[0]*xx[0]-xx[1]-xx[0]*(xx[0]*xx[0]+xx[1]*xx[1])+
					lam*(xx[6]-xx[0]);
	xpy[1] = xx[0]+alp[0]*xx[1]-xx[1]*(xx[0]*xx[0]+xx[1]*xx[1])+
					lam*(xx[7]-xx[1]);
	xpy[2] = alp[1]*xx[2]-xx[3]-xx[2]*(xx[2]*xx[2]+xx[3]*xx[3])+
					lam*(xx[6]-xx[2]);
	xpy[3] = xx[2]+alp[1]*xx[3]-xx[3]*(xx[2]*xx[2]+xx[3]*xx[3])+
					lam*(xx[7]-xx[3]);
	xpy[4] = alp[2]*xx[4]-xx[5]-xx[4]*(xx[4]*xx[4]+xx[5]*xx[5])+
					lam*(xx[6]-xx[4]);
	xpy[5] = xx[4]+alp[2]*xx[5]-xx[5]*(xx[4]*xx[4]+xx[5]*xx[5])+
					lam*(xx[7]-xx[5]);

	xpy[6] = eps0*((xx[0]+xx[2]+xx[4])/3.0 - xx[6]);
	xpy[7] = eps0*((xx[1]+xx[3]+xx[5])/3.0 - xx[7]);

	return xpy;
}


/* Main computational function. */
void lcfunc(double *xx, double *oln, double lam, double *alp, 
			double eps0, double *outy)
{
	double *space, *xmn, *ff;
	double nm;
	int ii, jj, kk;

	/* Allocate some memory. */
	space=mxMalloc(64*sizeof(double));		/* Tangents, then normals. */
	xmn=mxMalloc(8*sizeof(double));			/* Average x. */

	/* Set the tangent and average vectors. */
	for (ii=0; ii<8; ii++)	{

		jj=8*ii;
		space[ii]	=	xx[jj+1]+xx[jj+2]+xx[jj+5]+xx[jj+6]-
						xx[jj]-xx[jj+3]-xx[jj+4]-xx[jj+7];

		space[ii+8]	=	xx[jj+2]+xx[jj+3]+xx[jj+6]+xx[jj+7]-
						xx[jj+1]-xx[jj]-xx[jj+5]-xx[jj+4];

		space[ii+16]=	xx[jj+4]+xx[jj+5]+xx[jj+6]+xx[jj+7]-
						xx[jj]-xx[jj+1]-xx[jj+2]-xx[jj+3];

		xmn[ii]		=	( xx[jj]+xx[jj+1]+xx[jj+2]+xx[jj+3]+
							xx[jj+4]+xx[jj+5]+xx[jj+6]+xx[jj+7] )/8;
	}


	/* Fill in the rest of space with the normal vectors. */
	jj=23;
	for (ii=0; ii<5; ii++)	{

		for (kk=ii; kk<ii+40; kk+=5)	{
			jj++;
			space[jj]=oln[kk];
		}
	}

	/* Do my own QR factorization here. */
	for (ii=0; ii<64; ii+=8)	{

		/* Subtract off the projections. */
		for (jj=0; jj<ii; jj+=8)	{

			nm=0.0;
			for (kk=0; kk<8; kk++)						/* Project. */
				nm+=space[ii+kk]*space[jj+kk];

			for (kk=0; kk<8; kk++)						/* Divide. */
				space[kk+ii] -= nm*space[kk+jj];
		}

		/* Get the norm and divide. */
		nm=0.0;
		for (kk=0; kk<8; kk++)							/* Get norm. */
			nm+=space[ii+kk]*space[ii+kk];

		nm = 1/sqrt(nm);								/* Divide. */
		for (kk=0; kk<8; kk++)
			space[kk+ii] *= nm; 
	}

	/* Calculate the vector field. */
	ff=fcy(xmn, lam, alp, eps0);

	/* Multiply it all out. */
	for (ii=0; ii<5; ii++)	{
		outy[ii]=0.0;
		for (kk=0; kk<8; kk++)
			outy[ii]+=space[24+ii*8+kk]*ff[kk];
	}

	/* Clean up. */
	mxFree(space);
	mxFree(xmn);

}

