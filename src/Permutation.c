#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Permutation.h"
/*
 *Description: This function computes results for two sample permutation test.
 *Parameters:
 *input : Numeric vector in R representing data set to be permuted
 *factor: Factor vector in R representing assignment of data to groups.
 *Return:
 *A numeric vector containing p value, number of permutations and observed test statistics.
*/
SEXP calculate_pvalue(SEXP input_x,SEXP input_y){

	int ysize=LENGTH(input_y);
	int xsize=LENGTH(input_x);
	double* x=Calloc(xsize,double);
	double* y=Calloc(ysize,double);
	memcpy(x,REAL(input_x),sizeof(double)*xsize);
	memcpy(y,REAL(input_y),sizeof(double)*ysize);
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, 3));
	Result* result=calculatePvalue(x,y,xsize,ysize);
	REAL(ans)[0]=result->pvalue;
	REAL(ans)[1]=result->permutations;
	REAL(ans)[2]=result->statistics;
	Free(x);
	Free(y);
	UNPROTECT(1);
	return ans;
}

