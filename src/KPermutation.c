#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "Permutation.h"
/*
 *Description: This function computes results for K sample permutation test.
 *Parameters:
 *input : Numeric vector in R representing data set to be permuted
 *factor: Factor vector in R representing assignment of data to groups.
 *Return:
 *A numeric vector containing p value, number of permutations and observed test statistics.
*/
SEXP calculate_K_pvalue(SEXP input,SEXP factor){
	int size=LENGTH(input);
	SEXP ans;
	PROTECT(ans=allocVector(REALSXP,3));
	Result* result=calculateKPermPvalue(REAL(input),INTEGER(factor),size);
	REAL(ans)[0]=result->pvalue;
	REAL(ans)[1]=result->permutations;
	REAL(ans)[2]=result->statistics;
	UNPROTECT(1);
	return ans;
}

/*
 *Description: This function computes observed test statistics (F statistics) for a K sample permutation test.
 *Parameters:
 *input : Numeric vector in R representing data set to be permuted
 *factor: Factor vector in R representing assignment of data to groups.
*/
SEXP calculate_FStatistics(SEXP input,SEXP factor){
	int size=LENGTH(input);
	SEXP ans;
	PROTECT(ans=allocVector(REALSXP,1));
	REAL(ans)[0]=calculateFStatistics(REAL(input),INTEGER(factor),size);
	UNPROTECT(1);
	return ans;
}
/*
 *Description: This function computes observed test statistics (reduced statistics) for a K sample permutation test.
 *Parameter:
 *input : Numeric vector in R representing data set to be permuted
 *factor: Factor vector in R representing assignment of data to groups.
*/
SEXP calculate_ReducedStatistics(SEXP input,SEXP factor){
	int size=LENGTH(input);
	SEXP ans;
	PROTECT(ans=allocVector(REALSXP,1));
	REAL(ans)[0]=calculateReducedStatistics(REAL(input),INTEGER(factor),size);
	UNPROTECT(1);
	return ans;
}




