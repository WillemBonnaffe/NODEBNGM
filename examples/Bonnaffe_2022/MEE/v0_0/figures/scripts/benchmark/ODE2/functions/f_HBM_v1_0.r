#############
## f_HBM.r ##
#############

## goal: function for bayesian analysis

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update: 
## 13-03-2019 - cleaned
## 28-03-2019 - added extra functions (e.g. chainList.write(), ...) + normalized function names
## 08-04-2019 - added extra autocor function
## 20-04-2019 - added/corrected the unlist function
## 21-04-2019 - added a polygon.gamma/truncnorm functions for plotting purposes (might move it to a graphical display thingy)
## 27-03-2022 - created v0_5 
##            - introduced function argmax
## 21-06-2022 - created v0_6
##            - introduced significance in summary tab
## 22-08-2022 - created v0_7
##            - cleaned code
## 11-12-2022 - added na removal in quantiles function

#####################
## chainList.apply ##
#####################

## goal: apply function to chain list

## args:
## chainList - list - list of chains
## f         - f    - function
chainList.apply = function(chainList,f)
{
    chain = chainList.unlist(chainList)
    f_ensemble  = t(apply(chain,1,f))
    f_q0.05 = apply(f_ensemble,2,quantile,p=c(0.05),na.rm=T)
    f_q0.95 = apply(f_ensemble,2,quantile,p=c(0.95),na.rm=T)
    f_mean  = apply(f_ensemble,2,mean,na.rm=T)
    return(list("f_mean"     = f_mean,
    			"f_q0.05"    = f_q0.05,
    			"f_q0.95"    = f_q0.95,
    			"f_ensemble" = f_ensemble))
}

#
###

##########################
## chainList.argmaxPost ##
##########################

## goal: return maximum a posteriori

## args:
## chainList - list - list of chains
chainList.argmaxPost = function(chainList)
{
    chainTab = chainList.unlist(chainList)
    MaP      = chainTab[which.max(chainTab[,1]),-1]
    return(MaP)
}

#
###

######################
## chainList.unlist ##
######################

## goal: convert chain list into a matrix

## args:
## chainList - list - list of chains
chainList.unlist = function(chainList)
{
    chainTab = NULL
    for(i in 1:length(chainList))
    {
        chainTab = rbind(chainTab, chainList[[i]])
    }
    colnames(chainTab) = colnames(chainList[[1]])  
    return(chainTab)
}

#
###

#####################
## chainList.write ##
#####################

## goal: store a chainList object in a folder

## args:
## chainList    - list   - list of chains
## pathToFolder - string - path to folder containing chains folders
chainList.write = function(chainList, pathToFolder)
{
    ## create folder
    system(paste("mkdir ", pathToFolder, "/chains", sep=""))
    
    ## path to chains folder
    pathToChainsFolder = paste(pathToFolder,"/chains",sep="")
        
    ## check if there exists chains already
    chainIndexVect = 1:length(chainList) + length(list.files(pathToChainsFolder))
    
    ## for each chains
    for(i in 1:length(chainIndexVect))
    {
        ##
        system(paste("mkdir ",pathToChainsFolder,"/chain_", chainIndexVect[i], sep=""))
        write.table(chainList[[i]], paste(pathToChainsFolder,"/chain_", chainIndexVect[i],"/chain_", chainIndexVect[i], ".csv", sep=""), sep=";")
    }         
}

#
###

####################
## chainList.read ##
####################

## goal: read a list of chains stored in chain folder

## args:
## pathToFolder - string - path to the folder containing the chain folder
chainList.read = function(pathToFolder)
{
    ## path to chains folder
    pathToChainsFolder = paste(pathToFolder,"/chains",sep="")
    
    ## initiate
    chainList = list()
    
    ## discard outlier chains denoted by #
    if(length(grep("#",list.files(pathToChainsFolder)))>0)
    {

        fileVect = list.files(pathToChainsFolder)[-grep("#",list.files(pathToChainsFolder))]
        
    }else
    {
        fileVect = list.files(pathToChainsFolder)
    } 
    
    ## read chains
    for(i in 1:length(fileVect))
    {
        chainList[[i]] = read.table(paste(pathToChainsFolder,"/",fileVect[i],"/",fileVect[i],".csv",sep=""), sep=";", header=T)
    }
    
    return(chainList)
}

#
###

####################
## chainList.burn ##
####################

## goal: removes initial iterations of chain list

## args:
## chainList - list - list of chains
## burnin    - int  - number iterations to remove
chainList.burn = function(chainList, burnin)
{
    for(l in 1:length(chainList)){chainList[[l]] = chainList[[l]][-burnin,]}
    return(chainList)
}

#
###

####################
## chainList.thin ##
####################

## goal: sample chain list

## args:
## chainList - list - list of chains
chainList.thin = function(chainList)
{
    thinin = seq(1,nrow(chainList[[1]]),nrow(chainList[[1]])/1000); 
    for(l in 1:length(chainList)){chainList[[l]] = chainList[[l]][thinin,]}
    return(chainList)
}

#
###

##########################
## chainList.summaryTab ##
##########################

## goal: compute summary statistics (including significance level)

## args:
## chainList - list - list of chains 
chainList.summaryTab = function(chainList,returnSignif=F)
{
    ## cut chains in half
    chainList_final = list() ; j = 1
    for(i in 1:length(chainList))
    {
        chainList_final[[j]] = chainList[[i]][1:(nrow(chainList[[i]])/2),]; j = j+1
        chainList_final[[j]] = chainList[[i]][(nrow(chainList[[i]])/2):nrow(chainList[[i]]),]; j = j+1
    }
    chainList = chainList_final
    
    ## unlist chains
    chains = NULL ; for(i in 1:length(chainList)){chains = rbind(chains, chainList[[i]])}
    
    K = nrow(chainList[[1]])
    w = apply(X=matrix(data=unlist(lapply(X=chainList, FUN=function(x){apply(X=x, MARGIN=2, FUN=var)})), byrow=T, ncol=ncol(chainList[[1]])), MARGIN=2, FUN=mean)
    b = K * apply(X=matrix(data=unlist(lapply(X=chainList, FUN=function(x){apply(X=x, MARGIN=2, FUN=mean)})), byrow=T, ncol=ncol(chainList[[1]])), MARGIN=2, FUN=var)
    r_hat = sqrt(1 + 1/K * (b/w - 1))
    # p = apply(chains,2,FUN=function(X){length(X[X>0])/length(X[X<0])})
	
	## estimates table
	estimatesTab = cbind(
        MaP  = as.vector(t(chains[which.max(chains[,1]),])),
        mean = apply(X=chains, FUN=mean, MARGIN=2),
        sd   = apply(X=chains, FUN=sd, MARGIN=2),
        t(apply(X=chains, FUN=quantile, MARGIN=2, probs=c(0.05,0.5,0.95))),
        r_hat)
    
    ## significant level
    if (returnSignif == T)
    {
	    signif       = gsub("1",x=gsub(pattern="0",x=t(apply(estimatesTab[,c(4,6)]>0,1,diff)),replacement="*"),"ns")
	    estimatesTab = data.frame(cbind(round(estimatesTab,4),signif))
    }

    ## terminate
    return(list(estimates=estimatesTab,covmat=cov(chains[,-1])))
}

#
###

#########################
## chainList.bayesPlot ##
#########################

## goal: visualise parameter confidence intervals 

## args: 
## chainList    - list - list of chains
## logTransform - bool - whether parameters should be log transformed
chainList.bayesPlot = function(chainList, logTransform=F,labels=NULL,main=NULL)
{
    ##
    estimatesTab = chainList.summary(chainList)[["estimates"]]
    estimatesTab = estimatesTab[-1,]
    estimatesTab = estimatesTab[order(estimatesTab[,"mean"]),]
    
    ## plotting parameters
    par(mar=c(5.1,6.1,4.1,4.1))

    ##
    if(is.null(labels)) labels = rownames(estimatesTab)
  
    ##
    if(logTransform == T)
    {
        ##
        estimatesTab = log(estimatesTab)
        
        ##
        plot(0:nrow(estimatesTab),cex=0, bty="n",xlim=c(min(c(0,estimatesTab)),max(c(0,estimatesTab))), xlab="Estimates", ylab="", yaxt="n",xaxt="n",main=main)
        
        ##
        for(i in 1:nrow(estimatesTab))
        {
            points(estimatesTab[i,"mean"],i, pch=1); # lines(log(c(exp(estimatesTab[i,"mean"])-2*exp(estimatesTab[i,"sd"]),exp(estimatesTab[i,"mean"])+2*exp(estimatesTab[i,"sd"]))),c(i,i),lwd=2, col="red"); 
            lines(estimatesTab[i,c("5%","95%")],c(i,i),lwd=2,col="red"); # lines(estimatesTab[i,c("0%","100%")],c(i,i))
        }
        
        ##
        axis(side = 2, at = 1:nrow(estimatesTab), labels = rownames(estimatesTab), las=1)
        
        ## 
        order_max = round(log(exp(max(estimatesTab)))/log(10))
        order_min = round(log(exp(min(estimatesTab)))/log(10))
        
        for(o in (order_min+1):order_max)
        {
            X = seq(0.1,1,1/10)*10^o
            axis(side=1, at=log(X), labels=rep("",length(X)))
            axis(side=1, at=log(max(X)), labels=max(X))
        }
        #
    } else
    {
        ##
        plot(0:nrow(estimatesTab),cex=0, bty="n",xlim=c(min(c(0,estimatesTab)),max(c(0,estimatesTab))), xlab="Estimates", ylab="", yaxt="n",main=main)
	    lines(c(0,0),c(1,nrow(estimatesTab)),lty=3)
        
        ##
        for(i in 1:nrow(estimatesTab))
        {
            points(estimatesTab[i,"mean"],i); # lines(c(estimatesTab[i,"mean"]-2*estimatesTab[i,"sd"],estimatesTab[i,"mean"]+2*estimatesTab[i,"sd"]),c(i,i),lwd=2, col="red"); 
            lines(estimatesTab[i,c("5%","95%")],c(i,i),col="red", lwd=2); # lines(estimatesTab[i,c("0%","100%")],c(i,i))
        }
        
        ##
        axis(side = 2, at = 1:nrow(estimatesTab), labels = labels, las=1)
    }
    ## plotting parameters
    par(mar=c(5.1,5.1,4.1,4.1))
}

#
###

########################
## chainList.postPlot ##
########################

## goal: plot posterior distrubtion

## args:
## CMChainList - list - list of chains
## nPoints     - int  - number of points per plot
chainList.postPlot = function(chainList, nPoints)
{
    chains = NULL ; for(i in 1:length(chainList)){chains = rbind(chains, chainList[[i]])}
    par(mfrow=c(ncol(chains),ncol(chains)), mar=c(0.5,0.5,0.5,0.5), oma=c(2,2,2,2))
    chains = chains[order(chains[,1]),]
    
    ## labels
    if(is.null(colnames(chains))){colnames(chains) = paste("X",1:ncol(chains))}
    
    for(i in 1:ncol(chains))
    {
        for(j in 1:ncol(chains))
        {
            
            if(i==ncol(chains)){xaxt = "s"}else{xaxt = "n"} ; if(j==1){yaxt = "s"}else{yaxt = "n"}
            
            if(i==j)
            {
                plot(1:2, cex=0, xlab="", ylab="", xlim=c(min(chains[,j]),max(chains[,j])), ylim=c(min(chains[,i]),max(chains[,i])), xaxt=xaxt, yaxt=yaxt)
                text(x=(min(chains[,j])+max(chains[,j]))/2, y=(min(chains[,i])+max(chains[,i]))/2 ,labels = paste(colnames(chains)[j]), cex=1.5)
            }
            
            if(j==i+1)
            {
                hist(chains[,j], freq=F, main="", xaxt="n", yaxt="n") # barplot(height=density(chains[,j])$y, main="", xaxt="n", yaxt="n", horiz = F)
                lines(density(chains[,j]), col="red")
            }
            
            if(i>j)
            {
                plot(1:2, cex=0, xlab="", ylab="", xlim=c(min(chains[,j]),max(chains[,j])), ylim=c(min(chains[,i]),max(chains[,i])), xaxt=xaxt, yaxt=yaxt)                        
                subset = seq(1,nrow(chains),nrow(chains)/nPoints)
                points(chains[subset,j], chains[subset,i], col=grey(level=0.0+0.9*((chains[subset,1]-min(chains[subset,1]))/(max(chains[subset,1])-min(chains[subset,1]))), alpha=0.75), pch=16)
                points(x=mean(chains[,j]), y=mean(chains[,i]), pch=1, col="red")
                points(x=chains[which.max(chains[,1]),j], y=chains[which.max(chains[,1]),i], pch=8, col="red")
            }
            
            if(j>i+1)
            {
                plot(1:2, cex=0, xlab="", ylab="", xlim=c(min(chains[,j]),max(chains[,j])), ylim=c(min(chains[,i]),max(chains[,i])), xaxt="n", yaxt="n", bty="n")                        
            }
        }
    }
    par(mfrow=c(1,1), mar=c(4,4,1,1))
}

#
###

#########################
## chainList.tracePlot ##
#########################

## goal: plot the traces of a chain list

## arg: 
## chainList - list - list of chains 
chainList.tracePlot = function(chainList)
{
    ##
    par_mar_old = par()$mar
    par_oma_old = par()$oma
    par(mfrow=c(4,2), mar=c(5,5,0,0), oma=c(0,0,1,1), cex.lab=1.5)
    
    ## assemble chains together
    chain = NULL
    for(j in 1:length(chainList))
    {
        chain = rbind(chain, cbind(nChain=j, chainList[[j]]),deparse.level=0)
    }
    
    ## visualize chains
    for(i in 2:dim(chain)[2])
    {
        ## plots
        plot(chain[,i],ylab=paste(colnames(chain)[i]), xlim=c(0,nrow(chain)/length(chainList)), type="l", col="white")
        for(j in 1:length(chainList)){lines(1:length(chain[chain[,"nChain"]==j ,i]), chain[chain[,"nChain"]==j,i], col=rainbow(length(chainList),start=0.5,end=0.75,alpha=0.5)[j])}

        ## hist
        hist(chain[,i],freq=F,main="", xlab=paste(colnames(chain)[i]))
        for(j in 1:length(chainList)){lines(density(chain[chain[,"nChain"]==j,i]), col=rainbow(length(chainList),start=0.5,end=0.75,alpha=0.5)[j])}
    }
    par(mfrow=c(1,1),mar=par_mar_old,oma=par_oma_old)
}

#
###

############
## acplot ##
############

## goal: plot autocorrelation with increasing order

## args:
## X     - vector - data vector
## order - int    - order up to which ac are calculated
## xlab  - string - label of x axis
## ylab  - string - label of y axis
acPlot = function(X,order,xlab="",ylab="",add=F,col="grey")
{
    orderVect = corVect = NULL
    for(i in 1:order)
    {
        corVect_local = NULL; X_local = X
        for(j in 1:10)
        {
            corVect_local = c(corVect_local,cor(X_local[-c((length(X_local)-i+1):length(X_local))],X_local[-c(1:i)]))
            X_local = X_local[-c(1:(length(X)/10))]
        }
        corVect = c(corVect, mean(corVect_local))
    }
    ##
    barplot(height = corVect, ylim=c(-1,1),xlab=xlab,ylab=ylab,cex.lab=1.5,col=col,add=add,border=adjustcolor("black",alpha=0))
}

#
###

######################
## chainList.acPlot ##
######################

## goal: plot acplots for each elements of a chain list

## args:
## chainList - list - list of chains
chainList.acPlot = function(chainList,legend=NULL)
{
    k = 1
    if(is.null(legend) ) mainVect = c("a.","b.","c.","d.","e.","f.","g.","h.","i.","j.","k.","l.")
    for(l in 1:length(chainList))
    {
        for(i in 2:dim(chainList[[l]])[2])
        {
            acPlot(X=chainList[[l]][,i],order=50,xlab=if(l == length(chainList))"Order"else"",ylab=if(i==2)"Autocorrelation"else"")
            legend("topright",legend=mainVect[k],bty="n",cex=1.5)
            k = k + 1
        }    
    }
}

#
###

###################
## chainList.ESS ##
###################

## goal: compute the effective sample size of a chain list

## args:
## chainList - list - list of chains
chainList.ESS = function(chainList)
{
    ESSList = list()
    l = 1
    for(l in 1:length(chainList))
    {
        ESSVect = NULL
        i = 2
        for(i in 2:dim(chainList[[l]])[2])
        {
            X = chainList[[l]][,i]
            orderVect = corVect = NULL
            o = 1
            for(o in 1:50)
            {
                corVect_local = NULL; X_local = X
                for(j in 1:10)
                {
                    corVect_local = c(corVect_local,cor(X_local[-c((length(X_local)-o+1):length(X_local))],X_local[-c(1:o)]))
                    X_local = X_local[-c(1:(length(X)/10))]
                }
                corVect = c(corVect, mean(corVect_local))
            }
            ESSVect = c(ESSVect,dim(chainList[[l]])[1]/(1+2*sum(abs(corVect))))
        }
        ESSList[[l]] = ESSVect
    }
    return(ESSList)
}

#
###


