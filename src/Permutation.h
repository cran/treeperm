
#include "Rheaders.h"

/*
 *Headers for C debugging mode
#define Calloc(a,b) ((b*)malloc(sizeof(b)*a))
#define Free free
#define DOUBLE_EPS .000001
*/
#define Node int*
/*
 * Result of a permutation test
*/
typedef struct Result{
	//P value
	double pvalue;
	//Total number of permutations
	double permutations;
	//Observed test statistics
	double statistics;
	
} Result;
/*
 *This structure represents assignment of data into groups.
*/
typedef struct Group{
	//How many elements each group has
	int* count;
	//An array of group indices
	int* factor;
	//How many groups
	int size;
} Group;

/*
extern double T;
//Global variables for K sample permutation
extern double VAR;
extern int SIZE;
extern Node INITIAL;
extern double TOTAL;
extern double PRECISION_ROOT;
//Global variables two permutation
extern int XSIZE;
extern int TSIZE;
extern int (*double_cmp)(const void*,const void*);
extern double* data;
extern double SUMX;
extern double CRITICAL;
*/
//================Common functions================
/*
extern int double_cmp1(const void* d1,const void* d2);
extern int double_cmp2(const void* d1,const void* d2);
extern Node create_node(int size);
extern Node retrieve_initial(int xsize);
extern double Abs(double x);
*/
//================Two sample special================
/*
extern Node retrieve_x(double* data,double* x);
extern double countExtreme(Node node,int limit,double sum);
extern double combinations(double total,double select);
*/
extern Result* calculatePvalue(double* x,double* y,int xsize,int ysize);

//================K sample special=================
/*
extern int int_cmp_greater(const void* i1,const void* i2);
extern int int_cmp_less(const void* i1,const void* i2);
extern int int_cmp_sort(const void* i1,const void* i2);
extern Group* fillGroups(int* indices,int size);
extern void cleanGroup(Group* group);
extern double calculateCurrentWithin(double* data,int* indices,Group* group,int size);
extern double calculateStatistics(double* data,int* indices,Group* group,int size);
extern double variance(double* data,Node node,int size);
extern double squareSum(double* data,Node node,int size);
extern double factorial(double value);
extern double multipleCombinations(int* count,int count_size);
extern int nextBetween(Node node,int value1,int value2,int size,int (*int_cmp)(const void*,const void*));
extern double groupMeanSquare(double* data,Node node,int size);
extern double groupSum(double* data,Node node,int size);
extern void reverse(Node node,int size);
double divideGroup(double* data,Node node,Group* group,int next,int limit,int space,double cumulativeT,double sum,int (*int_cmp)(const void*,const void*),int previous,double group_mean,double perm_mean,int loc);
*/
extern double calculateFStatistics(double* input,int* factor,int size);
extern double calculateReducedStatistics(double* input,int* factor,int size);
extern Result* calculateKPermPvalue(double* input,int* factor,int size);
