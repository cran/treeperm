\name{Permutation}
\alias{Permutation}
\alias{print.Permutation}
\title{
Two Sample Exact Permutation Test
}
\description{
	This function performs a two sample permutation test of data locations by returning a p value. It is a wrapper function for C level computations. Users are not expected to use this function directly because there is no error checking mechanisms. However, for purpose of performance testing and other demands, this function is visible in this package.
}
\usage{
Permutation(data, factor)
\method{print}{Permutation}(x,\dots)
}
\arguments{
  \item{x}{
The Permutation class to be printed.
}
  \item{data}{
	The set of all data to be permuted
}
  \item{factor}{
	The set of indices, indicating group assignments of values in vector data. This parameter should be in factor type..
}
  \item{\dots}{
Reserved for future.
}
}
\details{
	Details about implementations are located in Permutation.c.
}
\value{
\item{Statistics }{The observed test statistics for input data. It is the sum of all elements in first group.}
\item{Permutations }{The total number of permutations can be performed on input groups.}
\item{pvalue }{The p value of this two sample permutation test.}
}
\references{
Fisher, R. A. 1925, 'Statistical Methods for Research Workers', Oliver and Boyd, Edinburgh
}
\author{
Qiao Kang
}
\note{
Although the return value uses a group mean value as test statistics, the actual computation uses the group sum to speed up computation.
}
\examples{
value<-c(0,190,0,0,10,0,0,0,0,0,0,110,0,0,52,0,8,0,50,0,0,137,965,110)
label<-as.factor(c(1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2))
print(Permutation(value,label))
}
