##############
## f_DEMC.r ##
##############

## goal: Rcpp implementation of a single chain differential evolution monte carlo algorithm

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

##############
## INITIATE ##
##############

#
###

##########
## DEMC ##
##########

## goal: R implementation of the DE-MC algorithm

## args: 
# dTarget - function - target function to explore
# Theta_0 - vector   - vector of initial parameter values
# gamma   - float    - scaling parameter
# epsilon - matrix   - covariance matrix for the target distribution
# nIt     - int      - number of iterations
# msg     - bool     - whether messages should be printed at each iterations
DEMC = function(dTarget, Theta_0, gamma, epsilon, nIt, msg=TRUE)
{

    ####################
    ## INITIATE CHAIN ##
    
    ## initiate param vector
    Theta_local = Theta_0
    d = length(Theta_local)
    
    ## evalutate the local value of target at theta_0
    target_local = dTarget(paramVect=Theta_local)
    
    ## initiate chain
    Theta_chain = as.matrix(t(c(dLogPost=target_local, Theta_local)))
    
    ################
    ## ITERATIONS ##
    
    i = 1 
    k = accepted = count = 0
    for(i in 1:nIt)
    {
        
        k = k+1
        
        ## draw mean and dispersion parameters
        Theta_newLocal = Theta_local + gamma*(diff(Theta_chain[round(runif(2,max(1,k/2),k)),-1])) + epsilon*runif(d,-1,1)  
        
        ## calculate new target density
        target_newLocal = dTarget(paramVect=Theta_newLocal)
        
        ## calculate metropolis-hasting ratio
        r = exp(target_newLocal - target_local); if(is.nan(r)){r = 0}
        
        ## test to accept or reject MH ratio
        if(runif(1,0,1) < r) 
        {
            Theta_local = Theta_newLocal
            target_local = target_newLocal
            accepted = accepted + 1
        }
        
        ## update chain
        Theta_chain = rbind(Theta_chain, c(target_local, Theta_local), deparse.level=0)
        
        # print current state of chain and acceptance rate
        if(msg==TRUE&count==round((nIt/10)))
        {
            message(paste("It: ", k,
                          " | p: ", round(accepted/k, 2),
                          " | dPostLocal: ", round(target_local, 2),
                          " | dPostTarget: ", round(target_newLocal, 2),
                          " | Theta: ", 
                          round(Theta_local[1], 4),
                          round(Theta_local[2], 4),
                          round(Theta_local[3], 4),
                          round(Theta_local[4], 4)
            ))
            count = 0
        }
        count = count + 1
        
    }
    
    ## print acceptance rate
    message(paste("p = ",round(accepted/k, 2),sep=""))
    
    ## TERMINATION ##
    return(Theta_chain)
}

#
###
