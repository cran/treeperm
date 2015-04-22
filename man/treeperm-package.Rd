\name{treeperm-package}
\alias{treeperm-package}
\docType{package}
\title{
	Exact and Asymptotic K Sample Permutation Test
}
\description{
	A package that performs exact K sample permutation test of data locations
}
\details{
\tabular{ll}{
Package: \tab treeperm\cr
Type: \tab Package\cr
Version: \tab 1.6\cr
Date: \tab 2015-04-22\cr
License: \tab GPL-2\cr
}
An implementation of permutation tests in R, supporting both exact and asymptotic K sample test of data locations. The p value of exact tests is found using tree algorithms. Tree algorithms treat permutations of input data as tree nodes and perform constraint depth-first searches for permutations that fall into the critical region of a test systematically. Pruning of tree search and optimisations at C level enable exact tests for certain large data sets.
}
\author{
	Qiao Kang

Maintainer: Qiao <qk@st-andrews.ac.uk>

}
\keyword{ Permutation treeperm }
\seealso{
Permutation KPermutation plot.treeperm treeperm
}
\examples{
value<-c(0,190,0,0,10,0,0,0,0,0,0,110,0,0,52,0,8,0,50,0,0,137,965,110)
label<-as.factor(c(1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2))
data<-data.frame(label,value)
exact<-treeperm(value~label,frame=data,type="exact")
exact
summary(exact)
approximate<-treeperm(value~label,frame=data,type="approximate",size=999)
approximate
summary(approximate)
plot(exact,size=999)
}
