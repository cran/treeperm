\name{plot.treeperm}
\alias{plot.treeperm}
\title{
Plot Permutations
}
\description{
This function use large sample approximation to visualize a permutation class after user performed a treeperm test.
}
\usage{
\method{plot}{treeperm}(x, size, \dots)
}
\arguments{
  \item{x}{
The treeperm class to be visualized
}
  \item{size}{
The size of simulated distribution
}
  \item{\dots}{
Not necessarily required, this is reserved for future extensions.
}
}
\details{
This function is independent of the type of permutation (either exact or approximate) contained in the treeperm object. It will use the stored data and factor in treeperm create an asymptotic permutation distribution and plot that distribution.
}
\value{
A histogram for visualization.
%% ...
}
\author{
Qiao Kang
}
\note{
This function simulates a distribution for a given distribution. Please set the random seed in advance.
}

\seealso{
treeperm
}
\examples{
value<-c(0,190,0,0,10,0,0,0,0,0,0,110,0,0,52,0,8,0,50,0,0,137,965,110)
label<-as.factor(c(1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2))
data<-data.frame(label,value)
exact<-treeperm(value~label,frame=data,type="exact")
plot(exact,size=999)
function (x, size, ...) 
{
    ran <- GetDistribution(x, size)
    r <- hist(ran, breaks = 50, main = "Permutations estimated by Monto carlo method", 
        xlab = "F  statistics", ylab = "Frequency", col = "lightblue")
    top <- max(r$counts)
    points(x$Fstatistics, -top/150, type = "p", pch = 17, col = "red", 
        lwd = 3)
    text(x$Fstatistics, -top/50, cex = 0.65, "Observed F statistics")
  }
}
