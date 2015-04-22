\name{treeperm}
\alias{treeperm}
\alias{treeperm.formula}
\alias{treeperm.default}
\alias{print.treeperm}
\alias{summary.treeperm}
\alias{print.summary.treeperm}
\title{
Permutation Tests Using Tree Algorithms
}
\description{
This function will perform a K samnple permutation test of data locations. Users may choose either asymptotic test or exact test.
}
\usage{
\method{treeperm}{formula}(formula, frame = list(), type, size, \dots)
\method{treeperm}{default}(x,data,factor, type, size, \dots)
\method{print}{treeperm}(x,\dots)
\method{summary}{treeperm}(object,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Generic argument, used only in print.treeperm to represent the treeperm class.
}
  \item{formula}{
The formula to perform permutation test. The form should be value~label, where value is the set of data and label classifies value into groups.
}
  \item{frame}{
The dataframe storing vectors used in formula.
}
  \item{data}{
The data set to be permuted
}
  \item{factor}{
The index set that classifies parameter data into groups.
}
  \item{type}{
Which type of permutation the user want to use, either 'exact' or 'approximate'.
}
  \item{size}{
If the type is 'approximate', user must specify the size of permutations they want.
}
  \item{object}{
The treeperm class to be summarised.
}
  \item{\dots}{
Reserved for future.
}
}
\details{
This function calls KPermutation or Permutation depending on the levels of factor.
}
\value{
\item{result}{A permutation object returned by a permutation function, either class Permutation or class KPermutation depending the levels of argument factor.}
\item{call}{User's calling script to this function at highest level.}
\item{data}{The input data set.}
\item{factor}{The index set that classifies data into groups.}
}
\author{
Qiao
}
\seealso{
Permutation
KPermutation
}
\examples{
value<-c(0,190,0,0,10,0,0,0,0,0,0,110,0,0,52,0,8,0,50,0,0,137,965,110)
label<-as.factor(c(1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2))
data<-data.frame(label,value)
exact<-treeperm(value~label,frame=data,type="exact")
summary(exact)
approximate<-treeperm(value~label,frame=data,type="approximate",size=999)
summary(approximate)
}
