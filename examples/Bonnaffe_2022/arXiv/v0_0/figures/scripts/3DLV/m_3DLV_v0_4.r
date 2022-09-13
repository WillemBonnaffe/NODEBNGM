#############
## RSCRIPT ##
#############

## goal: perform assessment of quality of inference of NODE by gradient matching 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update:
## 08-04_2022 - created v0_0
##            - simplified code by generating and storing ground truth
## 05-05-2022 - created v0_1
##            - implemented nonlinear 3D Lotka-Volterra
## 13-05-2022 - created v0_2
##            - kept only one model (model 4)
##            - renamed file m_3DLV_v0_0.r
## 25-05-2022 - created v0_1
##            - changed model to make sure nonlinear dynamics are identifiable
##            - created v0_2
##            - tested a model where the prey gains resitance to the predator at high density
## 31-05-2022 - created v0_4

##############
## INITIATE ##
##############

## 
library(deSolve)

#
###

#################
## 3D LV MODEL ##
#################

## goal:

## define model
t_true     = seq(0,100,0.1)
Y_true_0   = c(1.0,1.5,0.25)
# p          = runif(8,0,2) # c(1.0,1.0,2.5,1.5,1.0,1.0,1.0,1.0)
p          = c(1.0,3.0,2.5,1.5,1.0,1.0,1.0,1.0) 
# load("out/parameters_slow_oscillations.RData")
#
ddt.Y_true = function(t,x,p) list(c((p[1]*(1-x[1]/p[2]) - p[3]*1/(1+p[8]*x[1])*x[2] - p[4]*1/(1+p[8]*x[1])*x[3]) * x[1],
                                    (p[3]*1/(1+p[8]*x[1])*x[1] - p[5]*x[3] - p[6])*x[2],
                                    (p[4]*1/(1+p[8]*x[1])*x[1] + p[5]*x[2] - p[7])*x[3]))
#
ddx.r_true = function(t,x,p) list(c( - p[1]/p[2] - p[3]*x[2]*(-1)*p[8]/((1+p[8]*x[1])^2) - p[4]*x[3]*(-1)*p[8]/((1+p[8]*x[1])^2),
                                     - p[3]*1/(1+p[8]*x[1]),
                                     - p[4]*1/(1+p[8]*x[1]),
                                     p[3]*(1/(1+p[8]*x[1])+x[1]*(-1)*p[8]/((1+p[8]*x[1])^2)),
                                     0,
                                     - p[5],
                                     p[4]*(1/(1+p[8]*x[1])+x[1]*(-1)*p[8]/((1+p[8]*x[1])^2)),
                                     p[5],
                                     0
                                     ))
## run ode
Yhat_true  = ode(y=Y_true_0,times=t_true,func=ddt.Y_true,parms=p)

## ground truth
Yhat_true     = Yhat_true[,2:4]
ddt.Yhat_true = t(apply(Yhat_true,1,function(x)unlist(ddt.Y_true(NULL,x,p))))
ddx.rhat_true = t(apply(Yhat_true,1,function(x)unlist(ddx.r_true(NULL,x,p))))

## visualise results
pdf("out/TS_3DLV.pdf")
#
par(cex.lab=1.25)
plot(t_true,Yhat_true[,1],col="black",ylim=c(0,max(Yhat_true)*1.25),cex=0,xlab="Time",ylab="Density")
lines(t_true,Yhat_true[,1],col="green",type="l")
lines(t_true,Yhat_true[,2],col="blue",type="l")
lines(t_true,Yhat_true[,3],col="red",type="l")
legend("topleft",legend=c("Algae","Flagellate","Rotifer"),lty=1,col=c("green","blue","red"),bty="n")
par(cex.lab=1)
#
## create TS 
s            = seq(0,100*10,10) + 1
TS           = cbind(t_true[s],Yhat_true[s,])
colnames(TS) = c("t","G","B","R")
#
## visualise TS
plot(data.frame(TS))
#
## visualise ground truth
layout(mat=matrix(1:9,ncol=3,byrow=F))
for(i in 1:3)
{
    plot(t_true,Yhat_true[,i])
    plot(t_true,ddt.Yhat_true[,i])
    plot(-1:1,xlim=c(min(t_true),max(t_true)),ylim=c(-3,3),cex=0)
    for(j in 1:3)
    {
        lines(t_true,ddx.rhat_true[,(i-1)*3+j])
    }
}
#
dev.off()

## save results
write.table(TS,file="out/TS_3DLV.csv",sep=";",row.names=FALSE)
LV_GT = list("t_true"=t_true,"Yhat_true"=Yhat_true,"ddt.Yhat_true"=ddt.Yhat_true,"ddx.rhat_true"=ddx.rhat_true) 
save(LV_GT,file="out/GT_3DLV.RData")
save(p,file="out/parameters.RData")

#
###
