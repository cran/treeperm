#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Permutation.h"
#include <stdio.h>



//=================Global variables=================
/*
 * This variable represents observed test statistics in K sample test
 * It represents mean of permutation distribution in two sample test
*/
double T;
//=================Global variables for K sample permutation=================
//size of data set
int SIZE;
//the memory space for tree expansion.
Node INITIAL;
//sum of current node
double TOTAL;
//precision tolerance
double PRECISION_ROOT;
//=================Global variables two sample permutation=================
//size of first group
int XSIZE;
//size of data set
int TSIZE;
//comparator for test statistics
int (*double_cmp)(const void*,const void*);
//data set
double* data;
//sum of first group (observed test statistics)
double SUMX;
//distance from mean of distribution to observed test statistics
double CRITICAL;

//=================Common functions=================
/*
 *==============description==============
 *Comparator for greater than.
 *==============parameters==============
 *d1: The first value to be compared.
 *d1: The second value to be compared.
 *==============return value==============
 *Result of comparation, 1 for true and -1 for false.
*/
int double_cmp1(const void* d1,const void* d2){
	return *((double*)d1)>*((double*)d2)?1:-1;
}
/*
 *==============description==============
 *Comparator for less than.
 *==============parameters==============
 *d1: The first value to be compared.
 *d1: The second value to be compared.
 *==============return value==============
 *Result of comparation, 1 for true and -1 for false.
*/
int double_cmp2(const void* d1,const void* d2){
	return *((double*)d1)<*((double*)d2)?1:-1;
}
/*
 *==============description==============
 *Allocate memory space for a node.
 *==============parameters==============
 *size: size of this node.
 *==============return value==============
 *The new node.
*/
Node create_node(int size){
	return Calloc(size,int);
}
/*
 *==============description==============
 *Fill a node with a sequence starting from 1.
 *==============parameters==============
 *xsize: size of this node.
 *==============return value==============
 *The new node.
*/
Node retrieve_initial(int xsize){
	Node node=create_node(xsize);
	int i;
	for(i=0;i<xsize;i++){
		node[i]=i;
	}
	return node;
}
/*
 *==============description==============
 * Find absolute value.
 *==============parameters==============
 * x: The value to be processed
 *==============return value==============
 * Absolute value of x.
*/
double Abs(double x){
	if(x<0){
		return -x;
	}
	return x;
}
//=================Two sample permutation functions=================
/*
 *==============description==============
 * Tree algorithm for permuting two groups of data.
 *==============parameters==============
 * node: Which node to be expanded.
 * limit: Start with which position from RHS
 * sum: Cumulativreturn value==============
 * Total number of valid nodes under (and including) this node.
*/
double countExtreme(Node node,int limit,double sum){
	int i=limit;
	double counter;
	//Direction conditioning and early pruning
	if(Abs(sum-T)<CRITICAL||(double_cmp==double_cmp1&&sum>=T)||(double_cmp==double_cmp2&&sum<T)){
		return 0;
	}else{
		counter=1;
	}
	//Expand right most elment.
	if(node[i]<TSIZE-1){
		if(limit==XSIZE-1||node[i]+1!=node[i+1]){
			//Memory swapping
			node[i]+=1;
			counter+=countExtreme(node,limit,sum-data[node[i]-1]+data[node[i]]);
			//Memory recover
			node[i]-=1;
		}

	}
	//Expand all permissible positions from right to left
	while(i--){
		if(node[i]+1!=node[i+1]){
			//Memory swapping
			node[i]+=1;
			counter+=countExtreme(node,i,sum-data[node[i]-1]+data[node[i]]);
			//Memory recover
			node[i]-=1;
		}
	}
	return counter;
}
/*
 *==============description==============
 *This function calculates possible combinations for selecting elements from a group.
 *Code has been optimised.
 *==============parameters==============
 * total: Total number of elements
 * select: How many elements selected from the group.
 *==============return value==============
 * total C select
*/
double combinations(double total,double select){
	//Check if inputs are valid
	if(total<0||select<0||total<select){
		return 0;
	}
	if(total==select){
		return 1;
	}
	double start=total-select;
	if(start<select){
		start=select;
	}
	double i=2;
	double loss=1;
	double result=++start;
	while(start++<total){
		result*=start;
		loss*=i;
		i++;
	}
	return result/loss;
}
/*
 *==============description==============
 * Main function for performing a two sample permutations test. Code has been optimised.
 *==============parameters==============
 * x: First sample
 * y: second sample
 * xsize: Size of x
 * ysize: Size of y
 *==============return value==============
 * Result of permutation: # of permutations, p value, observed test statistics
*/
Result* calculatePvalue(double* x,double* y,int xsize,int ysize){
	//Precision tolerance
	PRECISION_ROOT=sqrt(DOUBLE_EPS);
	//Memory allocations
	XSIZE=xsize;
	TSIZE=xsize+ysize;
	data=Calloc(TSIZE,double);
	double* pos=data;
	int i=XSIZE;
	SUMX=0;
	double total=0;
	while(i--){
		SUMX+=*x;
		*pos=*x;
		x++;
		pos++;
	}
	i=ysize;
	total=SUMX;
	while(i--){
		total+=*y;
		*pos=*y;
		y++;
		pos++;
	}
	//Determine observed test statistics
	T=total*XSIZE/TSIZE;
	//Determine critical bounds
	CRITICAL=Abs(T-SUMX)-PRECISION_ROOT;
	//Sort data set
	qsort(data,TSIZE,sizeof(double),double_cmp1);
	//Construct the first node, this is the only piece of memory required by by tree expansions.
	Node initial=retrieve_initial(XSIZE);
	//The first cumulative sum
	double minimum=0;
	for(i=0;i<XSIZE;i++){
		minimum+=data[i];
	}
	double_cmp=&double_cmp1;
	double counter=countExtreme(initial,XSIZE-1,minimum);

	//Start the same process from RHS of distribution for a two-way test
	qsort(data,TSIZE,sizeof(double),double_cmp2);
	minimum=0;
	for(i=0;i<XSIZE;i++){
		minimum+=data[i];
	}
	double_cmp=&double_cmp2;
	counter+=countExtreme(initial,XSIZE-1,minimum);

	Result* result=Calloc(1,Result);
	result->permutations=combinations(TSIZE,XSIZE);
	result->pvalue=counter/result->permutations;
	result->statistics=SUMX/XSIZE;
	Free(data);
	Free(initial);
	return result;
}

//=========================K sample permutation functions=========================

/*
 *==============description==============
 *Comparator for greater than.
 *==============parameters==============
 *i1: The first value to be compared.
 *i1: The second value to be compared.
 *==============return value==============
 *Result of comparation, 1 for true and 0 for false.
*/
int int_cmp_greater(const void* i1,const void* i2){
	return *((int*)i1)>*((int*)i2)?1:0;
}
/*
 *==============description==============
 *Comparator for less than.
 *==============parameters==============
 *i1: The first value to be compared.
 *i1: The second value to be compared.
 *==============return value==============
 *Result of comparation, 1 for true and 0 for false.
*/
int int_cmp_less(const void* i1,const void* i2){
	return *((int*)i1)<*((int*)i2)?1:0;
}
/*
 *==============description==============
 * Comparator for less than.
 *==============parameters==============
 * i1: The first value to be compared.
 * i1: The second value to be compared.
 *==============return value==============
 * Result of comparation, 1 for true and -1 for false.
*/
int int_cmp_sort(const void* i1,const void* i2){
	return *((int*)i1)<*((int*)i2)?1:-1;
}

/*
 *==============description==============
 * Create a group structure based on information implied in index (factor) set.
 *==============parameters==============
 * indices: index set to be inferred.
 * size: size of indices
 *==============return value==============
 * A group structure containing information of indices.
*/
Group* fillGroups(int* indices,int size){
	int i,j=0;
	int ngroup=1;
	int* copy=Calloc(size,int);
	//Make a copy
	memcpy(copy,indices,sizeof(int)*size);
	//Sort the copy
	qsort(copy,size,sizeof(int),int_cmp_sort);
	//Count group number
	for(i=1;i<size;i++){
		if(copy[i]!=copy[i-1]){
			ngroup++;
		}
	}
	//Create elements for group structure
	int* group=Calloc(ngroup,int);
	int* factor=Calloc(ngroup,int);

	group[j]=1;
	factor[j]=copy[0];
	//Fill those elements
	for(i=1;i<size;i++){
		if(copy[i]==copy[i-1]){
			//Count number of elements in a group
			group[j]++;
		}else{
			//Enter next group
			j++;
			group[j]=1;
			factor[j]=copy[i];
		}
	}
	//Wrap all elements with a group structure
	Group* result=Calloc(1,Group);
	result->count=group;
	result->factor=factor;
	result->size=ngroup;
	return result;
}

/*
 *==============description==============
 * Clean memory space for a group structure.
 *==============parameters==============
 * group: which group to be cleaned.
*/
void cleanGroup(Group* group){
	Free(group->count);
	Free(group->factor);
	Free(group);
}
/*
 *==============description==============
 * Calculate sum of variance within each group
 *==============parameters==============
 * data: data set
 * indices: group assignment information of data
 * group: group information of indices
 * size: size of data
 *==============return value==============
 * sum of variance within each group
*/
double calculateCurrentWithin(double* data,int* indices,Group* group,int size){
	int i,j;
	double sum,sumsquare,result=0;
	for(i=0;i<group->size;i++){
		sum=0;
		sumsquare=0;
		for(j=0;j<size;j++){
			if(indices[j]==(group->factor)[i]){
				sum+=data[j];
				sumsquare+=(data[j]*data[j]);
			}
		}
		result+=(sumsquare-sum*sum/((double)(group->count)[i]));
	}
	return result;
}

/*
 *==============description==============
 * Calculate observed test statistics in reduced form
 *==============parameters==============
 * data: data set
 * indices: group assignment information of data
 * group: group information of indices
 * size: size of data
 *==============return value==============
 * Observed test statistics in reduced form
*/

double calculateStatistics(double* data,int* indices,Group* group,int size){
	int i,j;
	double sum,result=0;
	for(i=0;i<group->size;i++){
		sum=0;
		for(j=0;j<size;j++){
			if(indices[j]==(group->factor)[i]){
				sum+=data[j];
			}
		}
		result+=sum*sum/((double)(group->count)[i]);
	}
	return result;
}

/*
 *==============description==============
 * Calculate variance of a node
 *==============parameters==============
 * data: data set
 * node: which node used
 * size: size of data
 *==============return value==============
 * Variance of current node
*/

double variance(double* data,Node node,int size){
	int i;
	double sum=0;
	double squaresum=0;
	for(i=0;i<size;i++){
		squaresum+=(data[node[i]]*data[node[i]]);
		sum+=data[node[i]];
	}
	return (squaresum-sum*sum/((double)size));
}

/*
 *==============description==============
 * Calculate sum of squared elements in a node
 *==============parameters==============
 * data: data set
 * node: which node used
 * size: size of data
 *==============return value==============
 * sum of squared elements in a node
*/
double squareSum(double* data,Node node,int size){
	double result=0;
	int i;
	for(i=0;i<size;i++){
		result+=data[node[i]];
	}
	return result*result;
}
/*
 *==============description==============
 * Calculate the factorial of a number
 *==============parameters==============
 * value: which value used
 *==============return value==============
 * value!
*/
double factorial(double value){
	double result=1;
	while(value>1){
		result*=value;
		value--;
	}
	return result;
}

/*
 *==============description==============
 * Calculate distinct permutations from k groups of elements
 *==============parameters==============
 * count: how large each group is
 *==============return value==============
 * distinct permutations from k groups of elements
*/

double multipleCombinations(int* count,int count_size){
	double result=1;
	double size=0;
	int i;
	for(i=0;i<count_size;i++){
		result*=factorial(count[i]);
		size+=count[i];
	}
	result=factorial(size)/result;
	return result;
}

/*
 *==============description==============
 * Calculate sum of data represented by current node.
 *==============parameters==============
 * data: dataset
 * node: which node used
 * size: size of node
 *==============return value==============
 * sum of data represented by current node
*/
double groupSum(double* data,Node node,int size){
	int j;
	double result=0;
	for(j=0;j<size;j++){
		result+=data[node[j]];
	}
	return result;
}
/*
 *==============description==============
 * Reverse all indices in current node
 *==============parameters==============
 * node: which node used
 * size: size of node
*/
void reverse(Node node,int size){
	int i,temp;
	for(i=0;i<size/2;i++){
		temp=node[i];
		node[i]=node[size-1-i];
		node[size-1-i]=temp;
		
	}
}
/*
 *==============description==============
 * This is the tree algorithm that can count the number of K sample permutation that fail into the critical region of its test using reduced statistics.
 * The tree structure follows from two sample case, but in multiple dimension by separating trees into layers.
 * The permutation tree structure is generated systematically using depth-first search on the indices.
 * The only tree prune happens at the second last layer. The property of F statistics (or reduced statistics) rejects the advantage 
 * of cutting the process in sequence as there exists valid nodes even if its ancestor is invalid. 
 *==============parameters==============
 * data: Array holding the actual value of data
 * node: Current group to be expanded
 * group: The group information describing how data is divided
 * limit: The maximum index of node that is allowed to be permuted.
 * space: The remaining length of data to be permuted.
 * cumulativeT: Current sum mean squares for all layers above
 * sum: The cumulative sum of for all layers above
 * int_cmp: Comparator that decide which way the permutation is going to progress.
 * loc: The index of the node to be swapped with index limit in horizontal tree expansion. When a new layer is reached, its value is set to be the index at the start of remaining groups.
*/

double divideGroup(double* data,Node node,Group* group,int next,int limit,int space,double cumulativeT,double sum,int (*int_cmp)(const void*,const void*),int previous,double group_mean,double perm_mean,int loc){
	//Append current layer of permutation to cumulative test statistics. Check if this is a new layer
	if(previous!=next){
		group_mean=groupSum(data,node,(group->count)[next]);
		perm_mean=(TOTAL-sum)*(group->count)[next]/space;
	}
	double current=cumulativeT+group_mean*group_mean/((double)(group->count)[next]);
	//If last layer is reached
	if(next==group->size-1){
		//Check is the current permutation is more extreme than observed.
		if(T-current>PRECISION_ROOT){
			return 0;
		}else{			
			return 1;
		}
	}
	//This is the direction conditioning that forces a tree search to reverse or stops if its direction when the the group reach another side.
	if((int_cmp==int_cmp_less&&group_mean>=perm_mean)||(int_cmp==int_cmp_greater&&group_mean<perm_mean)){
		return 0;
	}

	double result=0;

	//Prunning at different layers

	if(T-current>PRECISION_ROOT){
		//Progress to another layer
		result=divideGroup(data,node+(group->count)[next],group,next+1,(group->count)[next+1]-1,space-(group->count)[next],current,sum+group_mean,int_cmp,next,0,0,(group->count)[next+1]);
		//Two sided test only when the layer is above or equal to three, layer two is a special case.
		if(next<=group->size-3){
			//reverse the data by reversing the indices, switch comparator direction.
			if(int_cmp==int_cmp_less){
				int_cmp=int_cmp_greater;
			}else{
				int_cmp=int_cmp_less;
			}
			reverse(node+(group->count)[next],space-(group->count)[next]);
			result+=divideGroup(data,node+(group->count)[next],group,next+1,(group->count)[next+1]-1,space-(group->count)[next],current,sum+group_mean,int_cmp,next,0,0,(group->count)[next+1]);
			//reverse back
			reverse(node+(group->count)[next],space-(group->count)[next]);
			if(int_cmp==int_cmp_less){
				int_cmp=int_cmp_greater;
			}else{
				int_cmp=int_cmp_less;
			}
		}
	}else{
		//If the current sum square is bigger than test statistics, no need to try the rest of layers. Return their total number immediately using formula
		result=multipleCombinations((group->count)+next+1,(group->size)-next-1);
	}

	//Tree pruning at layer two.
	//It is impossible for this scheme to cut the tree in any other layers above the second last one unless a better test statistics can be suggested.
	if(result==0&&next==group->size-2){
		return 0;
	}

	int temp;
	//Exchange group element with the remaining data at current layer.
	if(loc<space){
		//Maximum limit case should check if the next shift is out of boundary
		if(limit==(group->count)[next]-1){
			temp=node[limit];
			node[limit]=node[loc];
			node[loc]=temp;
			result+=divideGroup(data,node,group,next,limit,space,cumulativeT,sum,int_cmp,next,group_mean+data[node[limit]]-data[node[loc]],perm_mean,loc+1);
			temp=node[limit];
			node[limit]=node[loc];
			node[loc]=temp;
		}else{
			if(int_cmp(node+loc,node+limit+1)&&int_cmp(node+limit,node+loc)){
				temp=node[limit];
				node[limit]=node[loc];
				node[loc]=temp;
				result+=divideGroup(data,node,group,next,limit,space,cumulativeT,sum,int_cmp,next,group_mean+data[node[limit]]-data[node[loc]],perm_mean,loc+1);
				temp=node[limit];
				node[limit]=node[loc];
				node[loc]=temp;
			}
		}
	}
	//Expand the permissible indices in the current node from RHS to LHS
	while(limit--){
		loc=(group->count)[next];
		if(int_cmp(node+loc,node+limit+1)&&int_cmp(node+limit,node+loc)){
			temp=node[limit];
			node[limit]=node[loc];
			node[loc]=temp;
			result+=divideGroup(data,node,group,next,limit,space,cumulativeT,sum,int_cmp,next,group_mean+data[node[limit]]-data[node[loc]],perm_mean,loc+1);
			temp=node[limit];
			node[limit]=node[loc];
			node[loc]=temp;
		}
	}
	return result;
}
/*
 *==============description==============
 * Calculate F test statistics for observed permutation
 *==============parameters==============
 * input: data set
 * factor: group assignments of input.
 *==============return value==============
 * F test statistics for observed permutation
*/
double calculateFStatistics(double* input,int* factor,int size){
	int* indices=Calloc(size,int);
	memcpy(indices,factor,sizeof(int)*size);
	Group* group=fillGroups(indices,size);
	double *data=Calloc(size,double);
	memcpy(data,input,sizeof(double)*size);
	Node initial=retrieve_initial(size);
	double var=variance(data,initial,size);
	if(var==0){
		return 0;
	}
	double ssw=calculateCurrentWithin(data,indices,group,size);
	double result= ((var-ssw)/((double)group->size-1))/(ssw/((double)(size-group->size)));
	cleanGroup(group);
	Free(initial);
	Free(data);
	Free(indices);
	return result;
}
/*
 *==============description==============
 * Calculate reduced test statistics for observed permutation
 *==============parameters==============
 * input: data set
 * factor: group assignments of input.
 *==============return value==============
 * reduced test statistics for observed permutation
*/
double calculateReducedStatistics(double* input,int* factor,int size){
	int* indices=Calloc(size,int);
	memcpy(indices,factor,sizeof(int)*size);
	Group* group=fillGroups(indices,size);
	double *data=Calloc(size,double);
	memcpy(data,input,sizeof(double)*size);
	
	double result= calculateStatistics(data,indices,group,size);
	cleanGroup(group);
	Free(data);
	Free(indices);
	return result;
}

/*
 *==============description==============
 * Main function for performing a K sample permutations test. Code has been optimised.
 *==============parameters==============
 * input: data set
 * factor: group assignments of input.
 * size: size of input
 *==============return value==============
 * permutation test results: p value, number of permutations, and observed test statistics (F test statistics)
*/
Result* calculateKPermPvalue(double* input,int* factor,int size){
	//Precision tolerance
	PRECISION_ROOT=sqrt(DOUBLE_EPS);
	//Copy data
	Result* result=Calloc(1,Result);
	int* indices=Calloc(size,int);
	memcpy(indices,factor,sizeof(int)*size);
	Group* group=fillGroups(indices,size);
	if(size==group->size){
		return 0;
	}
	double *data=Calloc(size,double);
	memcpy(data,input,sizeof(double)*size);
	//Create initial node to be expanded
	Node initial=retrieve_initial(size);
	INITIAL=initial;
	//Compute observed test statistics
	T=calculateStatistics(data,indices,group,size);
	double var=variance(data,initial,size);
	if(var==0){
		result->statistics=0;
	}else{
		double ssw=calculateCurrentWithin(data,indices,group,size);
		result->statistics=((var-ssw)/((double)group->size-1))/(ssw/((double)(size-group->size)));
	}
	//Sort data set
	qsort(data,size,sizeof(double),double_cmp1);
	//Global variable assignments
	SIZE=size;
	TOTAL=groupSum(data,initial,size);
	//LHS tree algorithm
	double count=divideGroup(data,initial,group,0,(group->count)[0]-1,size,0,0,int_cmp_less,-1,0,0,(group->count)[0]);
	double total=((double)multipleCombinations(group->count,group->size));
	double pvalue;
	reverse(initial,size);
	//RHS tree algorithm
	count+=divideGroup(data,initial,group,0,(group->count)[0]-1,size,0,0,int_cmp_greater,-1,0,0,(group->count)[0]);
	pvalue=count/total;
	result->pvalue=pvalue;
	result->permutations=total;
	//Release memory space
	cleanGroup(group);
	Free(initial);
	Free(data);
	Free(indices);
	return result;
}

