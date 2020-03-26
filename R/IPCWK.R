#' @title Kendall's Tau Partial Corr. for Survival Trait and Biomarkers
#' @description The function proposes the inverse probability-of-censoring weighted (IPCW) Kendall's tau to measure the association of the survival trait with biomarkers and Kendall's partial correlation to reflect the relationship of the survival trait with interaction variable conditional on main effects, as described in Wang and Chen (2020) <doi:10.1093/bioinformatics/btaa017>. 
#' @param data The \eqn{N \times P} matrix of data. There are \eqn{N} individuals in matrix, with one individual in each row. The \eqn{P} columns orderly included the observable times which are time-to-event or censoring times and without ties at the event times, the status is a binary variable with 1 indicating the event has occured and 0 indicating (right) censoring, and the (\eqn{P-2}) main predictors. Note that the missing values of predictors in the data are not allowed in this version.
#' @param standarize Setting to "FALSE", point out the original gene espressions profiles are not standarized. We are going to standarize the gene features automatically by IPCWK function. Setting to "TRUE", point out the original gene espressions profiles have been standarized. We maintain the original gene features matrix to do the following analysis.
#' @return Returns a list with components
#' @return \item{CORR}{The \eqn{3 \times K} proposed correlation matrix, where K is the dimension of all biomarkers include main and interaction effects with one method in each row. The first row means the IPCW Kendall's tau correlation; the second row means the original Kendall's tau correlation without considering censoring scheme and the final row means the proposed Kendall's partial correlation. Note that the column names of "CORR" matrix is gene index, in which "Ga" means the ath main gene and "Ga&Gb" means the ath main gene interact with the bth main gene.}
#' @references Wang JH, and Chen YH* (2020) Interaction Screening by Kendall's Partial Correlation for Ultrahigh-dimensional Data with Survival Trait. published in  Bioinformatics <doi:10.1093/bioinformatics/btaa017>.
#' @import MASS survival stats utils
#' @export
#' @examples
#' set.seed(123)
#' library(MASS)
#' library(survival)
#'
#' numbeta=50
#' N=100 
#' beta0=matrix(0, numbeta, 1)
#' sigma1=0.5^abs(outer(1:numbeta, 1:numbeta, "-"))
#' W1=mvrnorm(N, beta0, sigma1, tol=1e-8, empirical=FALSE)
#' Z1=W1^2-1
#'
#' ### produce quadratic and two-way interaction effects ###
#' tempZZ1=matrix(Z1^2, N, numbeta)
#' tempZZ2=model.matrix(~(.)^2 - . -1, data=data.frame(Z1))
#' AZ=cbind(Z1, cbind(tempZZ1, tempZZ2))
#'
#' ### identify true predictors: G1, G10, G1&G1, G10&G10, G1&G10, G10&G20 ###
#' a=3
#' beta=matrix(0, dim(AZ)[2], 1)
#' beta[1,]=-0.8*a
#' beta[10,]=a
#' beta[51,]=1.2*a
#' beta[60,]=a
#' beta[109,]=-1.2*a
#' beta[515,]=a
#'
#' ### simulate survival time follows linear transformation model ###
#' C=matrix(runif(N,0,1), N, 1)
#' ST=X=S=matrix(0, N, 1)
#' temp=rexp(N)
#' ST=as.numeric(0.5*log(2*temp*exp(-AZ%*%beta)+1))
#' X=pmin(ST, C) 
#' S=(ST==X)
#' survdata=cbind(X, S, Z1)
#'
#' ### perform IPCWK function ###
#' test=IPCWK(data=survdata, standarize="FALSE")
#' true=which(beta!=0)
#' sum(order(-abs(test$CORR[1,]))[1:20] %in% true) ### IPCW-tau
#' sum(order(-abs(test$CORR[2,]))[1:20] %in% true) ### Kendall's tau
#' sum(order(-abs(test$CORR[3,]))[1:20] %in% true) ### PC-IPCW-tau

                ##################################################################
                #                   CREATE IPCWK FUNCTION                        #
                ##################################################################

IPCWK=function(data,standarize){

#################################################################
#1
tt=function(x){xx=x[,1]*x[,2]}

#2
PAS=function(z) {
g=paste("G", z, sep = "")
return(g)
}

#3
PAS1=function(z) {
g=paste("G", z[1], "&" ,"G", z[1],  sep = "")
return(g)
}

#4
PAS2=function(z) {
g=paste("G", z[1], "&" ,"G", z[2],  sep = "")
return(g)
}
#################################################################

DATA=data; N=dim(DATA)[1]; numbeta=dim(DATA)[2]-2; Z1=matrix(DATA[,3:(numbeta+2)],N,numbeta)

###################################################
Cen=1-DATA[,2]
km=survfit(Surv(DATA[,1],Cen)~1)
kmT=summary(km)$time ; kmS=summary(km)$surv
tempkms=c(rep(0,N))
for (i in 1:N){
temp=sum(DATA[i,1]>=kmT)
if (temp==0) tempkms[i]=1 else tempkms[i]=kmS[temp]
}
kms=tempkms
####################################################
####################################################
if (standarize){
Z1=Z1
} else {
Z1=scale(Z1)
}
colnames(Z1)=PAS(c(1:numbeta))

tempZZ1=matrix(0,N,numbeta); GGnames=c(rep(0,numbeta))
tempZZ1=Z1^2
for (k in 1:numbeta){
GGnames[k]=PAS1(k)
}
colnames(tempZZ1)=GGnames

tempcombn=combn(c(1:numbeta),2)
numss=dim(tempcombn)[2]; GGnames=c(rep(0,numss))
tempZZ2=matrix(0,N,numss)
for (k in 1:numss){
temp=tempcombn[,k]
tempZZ2[,k]=as.matrix(tt(Z1[,temp]))
GGnames[k]=PAS2(temp)
}
colnames(tempZZ2)=GGnames

tempZZ=cbind(tempZZ1,tempZZ2)
tempZZ=scale(tempZZ)

AZ=cbind(Z1,tempZZ); X=DATA[,1]; S=DATA[,2]
survdata=cbind(X, S, AZ); DATA=data.frame(survdata)
##############################################################

w=c(DATA[,2]/kms^2); w[is.infinite(w)]=0; w[is.na(w)]=0
v=c(DATA[,1])
W=matrix(c(w),N,N)
V=matrix(c(v),N,N)
f=N^2-N

#############################################################
#          START  compute IPCW Kendall's Tau                #
#############################################################

Cor=matrix(0,3,dim(AZ)[2])
for (j in 1:dim(AZ)[2]){

#1, > >
X=matrix(c(DATA[,j+2]),N,N)
XX=sign(t(X)-X)
VV=sign(t(V)-V)

XX[XX==-1]=0
VV[VV==-1]=0
A1W=sum(W*XX*VV)/f
A1=sum(XX*VV)/f

#3, < >
XX=sign(t(X)-X)
VV=sign(t(V)-V)
XX[XX==1]=0; XX[XX==-1]=1
VV[VV==-1]=0
B1W=sum(W*XX*VV)/f
B1=sum(XX*VV)/f

############################################################
Cor[1,j]=2*A1W-2*B1W #### modified
Cor[2,j]=2*A1-2*B1   #### no-weighted
}
#############################################################
#                END compute IPCW Kendall's Tau             #
#############################################################

############################################################# 
#         START compute Kendall's partial correlation       #
#############################################################
######################################## PART 2 start
PC1=c(rep(0,numbeta))
for (k in 1:numbeta){
temp=cor(DATA[,2+(k+numbeta)],DATA[,2+k],method="kendall")
PC1[k]=(Cor[1,k+numbeta]-Cor[1,k]*temp)/((sqrt(1-Cor[1,k]^2))*(sqrt(1-temp^2)))
}
######################################## PART 2 end
######################################## PART 3 start
int=combn(c(1:numbeta),2); numss=dim(int)[2]
PC2=c(rep(0,numss))
for (k in 1:numss){
temp1=cor(DATA[,(2+numbeta*2+k)],DATA[,(2+int[1,k])],method="kendall")
aPC2=(Cor[1,k+(numbeta*2)]-Cor[1,int[1,k]]*temp1)/((sqrt(1-Cor[1,int[1,k]]^2))*(sqrt(1-temp1^2)))

a=int[2,k]; b=int[1,k]
temp2=cor(DATA[,(2+a)],DATA[,(2+b)],method="kendall")
temp3=cor(DATA[,(2+(numbeta*2+k))],DATA[,(2+a)],method="kendall")
temp4=cor(DATA[,(2+(numbeta*2+k))],DATA[,(2+b)],method="kendall")
bPC2=(Cor[1,a]-Cor[1,b]*temp2)/(sqrt(1-Cor[1,b]^2)*sqrt(1-temp2^2))
cPC2=(temp3-temp4*temp2)/(sqrt(1-temp4^2)*sqrt(1-temp2^2))

PC2[k]=(aPC2-bPC2*cPC2)/(sqrt(1-bPC2^2)*sqrt(1-cPC2^2))
}
Cor[3,]=c(Cor[1,1:numbeta],PC1,PC2)
Cor[is.na(Cor)]=0
colnames(Cor)=colnames(AZ)
rownames(Cor)=c("IPCW Kendall's tau","Original Kendall's tau","Kendlla's partial correlation")
####################################### PART 3 end
############################################################# 
#           END compute Kendall's partial correlation       #
#############################################################
IPCWKlist=list(CORR=Cor)
return(IPCWKlist)

} ### end function
