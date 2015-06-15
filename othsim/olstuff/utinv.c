#include "mex.h"
#include "matrix.h"

/******
	This is a function for inverting an U-T matrix.
*******/
void utinv(double *ut, const int dim, double *outy)
{
	int i, j, k;			/* Iterators. */
	double tmp;
	int col;

	/* Note that I'm not initializing the lower triangular part of outy.*/
	/* Fill in the diagonals. */
	for (i=0; i<dim; i++){
		col=i+i*dim;
		outy[col]=1.0/ut[col];	}

	/* Fill in the rest. */
	for (i=0; i<dim; i++)  {		/* Row.*/

		for (j=i+1; j<dim; j++)  {		/* Column.*/

			col=j*dim;

			/* Get the thing in parentheses. 3/26/03, p.1*/
			tmp=0.0;
			for (k=i; k<j; k++)
				tmp-=outy[i+k*dim]*ut[k+col];

			outy[i+col] = outy[j+col]*tmp;
		}
	}

			
	
}

/******
	This is the "gateway" routine for the upper-triangular-inverse function.
*******/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  	double *ut, *outy;
	int dim;
  
	/* Assign pointers to input and get the scalar. */
	ut = mxGetPr(prhs[0]);
	dim = mxGetM(prhs[0]);

	/* Set the output pointer to the output matrix. */
	plhs[0]=mxCreateDoubleMatrix(dim,dim,mxREAL);
	outy = mxGetPr(plhs[0]);
  
	/* Call the computational subroutine. */
	utinv(ut, dim, outy);
	
}

