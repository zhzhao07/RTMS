RTMS<- function(x,y,rp,method="3S",CS=NULL,mc=1)
{
  p=ncol(x); n=length(y);x=as.matrix(x);y=c(y)
  colnames(x)=c(paste("X",1:p,"",sep =""))
if((method %in% c("3S","4I","7M","XM"))==FALSE)
{stop("The specified 'method' is wrong.")}
if(rp>=p)
{stop("The number of recommended variables cannot exceed the total number of variables.")}

 if(method=="XM"){
if(mc>=length(CS))
{stop("The number of trimmed methods exceeds the total number of candidate methods.")}
if(all(CS %in% c("DC-SIS","SIRS","ISIS","ARM","AIC","BIC","RF"))==FALSE)
{stop("The names in the specified variable screening methods contain some errors.")}
}
  if(method=="3S"||method=="7M"||method=="XM")
  {
  R1=screenIID(x, y, method = "DC-SIS")$rank
  R2=screenIID(x, y, method = "SIRS")$rank
  R3=SIS(x, y, family='gaussian', tune='bic',nsis = p)$ix0

  a1=paste("X",1:length(R1),sep="")
  nr3=1:length(R3)
  names(nr3)= paste("X",R3,sep="")
  R30=nr3[a1]
  a=cbind(R1,R2,as.numeric(R30) )
  a2=rank((rowSums(a)-apply(a,1,max))/2)
  a32=a2

  names(R1)=names(R2)=names(a2)=a1
 b=names(sort(a2))
 b1=names(sort(R1))
 b2=names(sort(R2))
 b3=names(sort(R30))
 b33=paste("X",as.numeric(SIS(x, y, family='gaussian', tune='bic',nsis =rp)$ix0),sep="")

 a23=cbind(b1[1:rp],b2[1:rp],b33)
 A=unique(names(sort(a2[a23])))
 a32[1:length(A)]=A
 a32= unique(c(a23,b))
 RTMS3S=a32
 object=list(RTMS3S=RTMS3S,DCSIS=b1,SIRS=b2,ISIS=b3)
 }
  
 if(method=="4I"||method=="7M"||method=="XM"){
 Rs1=SOIL(x=x,y=y,family='gaussian',weight_type = "ARM",no_rep = 100)$importance
 Rs2=SOIL(x=x,y=y,family='gaussian',weight_type = 'AIC')$importance
 Rs3=SOIL(x=x,y=y,family='gaussian',weight_type = 'BIC')$importance


 S1=match(Rs1,Rs1[order(c(Rs1),decreasing =TRUE)])
 S2=match(Rs2,Rs2[order(c(Rs2),decreasing =TRUE)])
 S3=match(Rs3,Rs3[order(c(Rs3),decreasing =TRUE)])
       names(S1)=paste("X",1:length(S1),sep="")
       names(S2)=paste("X",1:length(S2),sep="")
       names(S3)=paste("X",1:length(S3),sep="")
 b4=names(sort(S1))  #ARM
 b5=names(sort(S2))  #AIC
 b6=names(sort(S3))  #BIC

data_rf=data.frame(y=y,x=x)
Rf=randomForest(y~.,data=data_rf,importance=TRUE)
     rfim=Rf$importance[,2]
     F1=match(rfim,rfim[order(c(rfim),decreasing =TRUE)])
     RF=F1 #1:M0
     names(RF)=paste("X",1:length(RF),sep="")
 b71=names(sort(RF))  #RF
 
 ab=cbind(S1,S2,S3,F1)
 ab2=rank((rowSums(ab)-apply(ab,1,max))/3)
 
 aa3=cbind(b4[1:rp],b5[1:rp],b6[1:rp],b71[1:rp])
 AB=unique(names(sort(ab2[aa3])))
 aa32=ab2
 aa32[1:length(AB)]=AB
 aa32= unique(c(aa3,names(sort(ab2))))
 RTMS4I=aa32
 object=list(RTMS4I=RTMS4I,ARM=b4,AIC=b5,BIC=b6,RF=b71)
 }
 
 if(method=="7M"){
 RRR=cbind(a23,cbind(b4,b5,b6,b71)[1:rp,])
 d=cbind(R1,R2,as.numeric(R30),S1,S2,S3,F1)
 d2=(sort((rowSums(d)-apply(d,1,max))/6))

     B=unique(names(sort(d2[RRR])))
     d32=1:p
     d32[1:length(B)]=B
     d23= unique(c(B,names(d2)))
     D=d23
     b7=D  #RF
object=list(RTMS7M=b7,DCSIS=b1,SIRS=b2,ISIS=b3,ARM=b4,AIC=b5,BIC=b6,RF=b71)
}
     
     
  if(method=="XM"){
 XR=cbind(b1[1:rp],b2[1:rp],b33[1:rp],b4[1:rp],b5[1:rp],b6[1:rp],b71[1:rp])
 d=cbind(R1,R2,as.numeric(R30),S1,S2,S3,F1)
colnames(XR)=colnames(d)=c("DC-SIS","SIRS","ISIS","ARM","AIC","BIC","RF")
d=d[,CS];
d2=sort((rowSums(d)-rowSums((t(apply(d,1,function(x) sort(x, decreasing = TRUE))))[,1:mc]))/(ncol(d)-mc))

#CS=CS
#d2= d2[,CS]

     B=unique(names(sort(d2[XR])))
     d32=1:p
     d32[1:length(B)]=B
     d23= unique(c(B,names(d2)))
     D=d23
     b7=D  #RF
object=list(RTMSXM=b7,d=cbind(b1[1:rp],b2[1:rp],b33[1:rp],b4[1:rp],b5[1:rp],b6[1:rp],b71[1:rp]))
}

 class(object)    = "RTMS"
 object
}
