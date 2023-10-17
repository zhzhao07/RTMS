\name{RTMS}
\alias{RTMS}
\title{Recommendation and Trimmed Mean Score}
\description{
This package combines multiple screening methods by retaining a certain number of variables from each method and then taking a trimmed average of the rankings from the rest methods, resulting in a more informative variable combination.}
\usage{
RTMS(x,y,rp,method,CS,mc)
}
\arguments{
  \item{x}{Matrix of predictors.}

  \item{y}{Response variable.}

  \item{rp}{The number of recommendate predictors. The default is 3.}

  \item{method}{The parameter can take one of four values: "3S," "4I," "7M," and "XM." 
- "3S" corresponds to Recursive Feature Elimination with the following three variable selection methods: ISIS, SIRS, and DC-SIS.
- "4I" corresponds to Recursive Feature Elimination with the following four variable importance ranking methods: SOIL-AIC, SOIL-BIC, SOIL-ARM, and RF.
- "7M" corresponds to Recursive Feature Elimination with all seven variable selection methods.
- "XM" corresponds to Recursive Feature Elimination with at least two of the seven variable selection methods. The default is "3S."}

  \item{CS}{A selection of at least two characters from the set c("DC-SIS", "SIRS", "ISIS", "ARM", "AIC", "BIC", "RF"). This vector indicates the variable importance ranking methods used when the "method" parameter is set to "XM" for constructing RTMS.}

  \item{mc}{The number of lower-ranking methods trimmed.}

}
}

\details{
See the paper provided in Reference section.
}

\value{
A \code{'RTMS'} object is retured. The components are:

\item{RTMS}{The combines rankings.}
\item{DCSIS}{The DCSIS rankings.}
\item{SIRS}{The SIRS rankings.}
\item{ISIS}{The ISIS rankings.}

}



\examples{

### Example usage:
Below is an example of using the function RTMS:
#generate simulation data
library(RTMS)
n<-100
p<-1000
beta<-c(3,1.5,0,0,2,rep(0,p-5))
b0<-1
x<-matrix(rnorm(n*p,0,1),nrow=n,ncol=p)
e<-rnorm(n,0,3)
y<-x\%*\%beta+b0+e

R3  <- RTMS(x, y, rp=3, method="3S", mc=1)$RTMS
#the Recommendation and Trimmed Mean Score
R3  <- R3$RTMS
#the DCSIS Ranking
RDCSIS  <- R3$DCSIS

#the SIRS Ranking
RSIRS  <- R3$SIRS

#the ISIS Ranking
RISIS  <- R3$ISIS
}


