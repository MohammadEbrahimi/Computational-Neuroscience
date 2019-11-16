#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <omp.h>



void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	 /* Check for proper number of arguments. */
 	if(nrhs!=2) {
    	mexErrMsgTxt("Two input required.");
  	} else if(nlhs>2) {
    	mexErrMsgTxt("Too many output arguments.");
  	}
	
	/*Find number of rows and columns in input data*/
	int k = mxGetM(prhs[0]);
	int p = mxGetN(prhs[0]);			
	/*Get pointer to input data*/
	double* x = mxGetPr(prhs[0]);
    double* f=mxGetPr(prhs[1]);

	/*create space for output*/
 	 plhs[0] = mxCreateDoubleMatrix(k,p, mxREAL);

	/*Get pointer to output array*/
	double* y = mxGetPr(plhs[0]);
	/*Return sum to output array*/
    int i;
	#pragma omp parallel for shared(x,f,y) private(i)
	for(i = 0;i<p;i++)
	{
        int j;int ind=i*k;
        double temp=x[ind];
        int rec=ind;
        int mnum=1;
        for(j = 0;j<k;j++)
        {
           
         if (f[ind+j]>0){
             temp+=x[ind+j+1];
             mnum+=1;
         }
         else{
             int ii;
             for(ii=rec;ii<rec+mnum;ii++){
                 double fmean=temp/mnum;
                 y[ii]=fmean;
             }
            temp=x[ind+j+1]; 
            mnum=1;
            rec=ind+j+1;
         }
         
        }
	     
	}
    
   
}


