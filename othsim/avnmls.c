#include "mex.h"

/* Here is the main function for averaging the normal vectors. */
/* Note that this routine requires enns, as well as pts */

void avnmls(double *pts, double *enns, double *avenns)
{

	int ii, jj, kk, nn, pp;
	int iisc, ipsc, jjsc, jpsc, kp;
	int pos;

	/* Stupid type forcing. */
	int pt0 = (int)pts[0];
	int pt1 = (int)pts[1];
	int pt2 = (int)pts[2];
	int ibk = pt1*pt2;
	int aps = pt0*pt1*pt2;

	int crnr[8];						/* Indices of box corners. */

	/* Don't indent to save space. */
	for (ii=0; ii<pt0; ii++)	{
		
		iisc=ii*ibk;
		ipsc=((ii+1)%pt0)*ibk;

	for (jj=0; jj<pt1; jj++)	{

		jjsc=jj*pt2;
		jpsc=((jj+1)%pt1)*pt2;

	for (kk=0; kk<pt2; kk++)	{

		kp=(kk+1)%pt2;

		/* Now we have the point. */
		/* Build the boxes that the point will affect. */
		crnr[0]=iisc+jjsc+kk;	crnr[4]=iisc+jjsc+kp;
		crnr[1]=ipsc+jjsc+kk;	crnr[5]=ipsc+jjsc+kp;
		crnr[2]=ipsc+jpsc+kk;	crnr[6]=ipsc+jpsc+kp;
		crnr[3]=iisc+jpsc+kk;	crnr[7]=iisc+jpsc+kp;

		/* Store the normals. */	
		for (nn=0; nn<40*aps; nn+=aps)	{
		
			/* Sum it up. */
			avenns[crnr[0]+nn]=0.0;
			for (pp=0; pp<8; pp++)
				avenns[crnr[0]+nn]+=enns[crnr[pp]+nn];
		}

	}
	}
	}

}


/* "Gateway" routine. */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  	double *enns, *avenns;
	double *pts;			/* I don't like this. */
	int aps;
  
	/* I'm not doing ANY error checking on the inputs. Be careful. */

	/* Assign pointers to input and get the scalar. */
	pts = mxGetPr(prhs[0]);
	enns = mxGetPr(prhs[1]);
	
	aps = mxGetM(prhs[1]);

	/* Set the output pointer to the output matrix. */
	plhs[0] = mxCreateDoubleMatrix(aps, 40, mxREAL);
	avenns = mxGetPr(plhs[0]);

	/* Call the computational subroutine. */
	avnmls(pts, enns, avenns);

}

