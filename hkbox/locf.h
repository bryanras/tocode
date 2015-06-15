#include "mex.h"

/* This is the vector field calculation. */

/* This function returns the vector field. */
double *fcy(double *xx, double b, double c, double d, double lam)
{
    double dit;
	double *outy;

	/* Declare the memory. */
	outy=mxMalloc(3*sizeof(double));

    /* Don't do this calculation more than once. */
    dit=xx[2]+d*(1-xx[2]*xx[2]);

    /* This is everything. */
    outy[0] = (lam-b)*xx[0]-c*xx[1]+xx[0]*dit;
    outy[1] = c*xx[0]+(lam-b)*xx[1]+xx[1]*dit;
    outy[2] = lam*xx[2]-xx[0]*xx[0]-xx[1]*xx[1]-xx[2]*xx[2];


	return outy;

}

