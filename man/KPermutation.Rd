\name{KPermutation}
\alias{KPermutation}
\alias{print.KPermutation}
\title{
K Sample Exact Permutation Test
}
\description{
	This function performs a K sample permutation test of data locations by returning a p value. It is a wrapper function for C level computations. Although it supports two sample permutation test, users should really use the Permutation function because it is much faster in handling two sample case. Users are not expected to use this function directly because there is no error checking mechanisms. However, for purpose of performance testing and other demands, this function is visible in this package.
}
\usage{
KPermutation(data, factor)
\method{print}{KPermutation}(x,\dots)
}

\arguments{
  \item{x}{
The KPermutation class to be printed.
}
  \item{data}{
	The set of all data to be permuted
}
  \item{factor}{
	The set of indices, indicating group assignments of values in vector data. This parameter should be in factor type.
}
  \item{\dots}{
Reserved for future.
}
}
\details{
	Details about implementations are located in KPermutation.c.
}
\value{
\item{Statistics }{The observed test statistics for input data. F statistics is used in this case.}
\item{Permutations }{The total number of permutations can be performed on input groups.}
\item{pvalue }{The p value of this two sample permutation test.}
}
\references{
Ernst M. D 2004, 'Permutation Methods: A Basis for Exact Inference', Statistical Science, Vol. 19, No. 4, p 676-685
}
\author{
Qiao Kang
}
\note{
Although the return value contains F statistics, the actual computation uses a reduced statistics to speed up computation.
}
\examples{
value<-c(0,190,0,0,10,0,0,0,0,0,0,110,0,0,52,0,8,0,50,0,0,137,965,110)
label<-as.factor(c(1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2))
print(KPermutation(value,label))
function (data, factor) 
{
    factor <- as.factor(factor)
    p <- .Call("calculate_K_pvalue", data, as.integer(factor))
    result <- list(Fstatistics = p[3], Permutations = p[2], pvalue = p[1])
    class(result) <- "KPermutation"
    result
  }
}
