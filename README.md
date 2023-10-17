# RTMS
This package combines multiple screening methods by retaining a certain number of variables from each method and then taking a trimmed average of the rankings from the remaining methods, resulting in a more informative variable combination.   
### Authors: 
Zhihao Zhao `<zhzhao@cueb.edu.cn>`, Li Wen `<wlwendy1008@163.com>`, Yuhong Yang `<yyang@stat.umn.edu>`.  
### Installation
To install this package in R, run the following commands:  
`library(devtools)`  
`devtools::install_github("zhzhao07/RTMS")`  



### Example usage:
Below is an example of using the function RTMR:  
```#generate simulation data  
library(RTMS)  
n    <- 50  
p    <- 8  
beta <- c(3,1.5,0,0,2,0,0,0)  
b0   <- 1  
x    <- matrix(rnorm(n*p,0,1),nrow=n,ncol=p)  
e    <- rnorm(n,0,3)  
y    <- x%*%beta+b0+e  

#the Recommendation and Trimmed Mean Ranking
R3  <- RTMS(x, y, 3)$RTMS
#the DCSIS Ranking
RDCSIS  <- R3$DCSIS

#the SIRS Ranking
RSIRS  <- R3$SIRS

#the ISIS Ranking
RISIS  <-R3$ISIS

