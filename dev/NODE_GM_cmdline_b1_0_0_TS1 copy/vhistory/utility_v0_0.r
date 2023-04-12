###############
## f_utility ##
###############

## goal: utility function repository


###############
## FUNCTIONS ##
###############

## .plot.DIN
##
## goal: plot the dynamical interaction network of the system
##
## input:
## effectsMat - matrix - matrix of pairwise effects between system variables (e.g. row 1 col 2 is the effect of variable 2 on variable 1)
## weightsMat - matrix - matrix of pairwise weights of the effects between system variables (e.g. row 1 col 2 corresponds to the contribution of variable 2 on variable 1)
## labels     - vector - vector of the names of the variables in the matrix
##
## output:
.plot.DIN = function(effectsMat,weightsMat,labels)
{
	## dimensions
	N = dim(effectsMat)[1]

	## scale effects and contributions
	effectsMat = (effectsMat>0)*1

	## angles
	theta = seq(0,2*pi,(2*pi)/N)
	x = cos(theta)
	y = sin(theta)
	x_ = cos(theta+(2*pi*0.05))
	y_ = sin(theta+(2*pi*0.05))
	x__ = 1.25*cos(theta+(2*pi*0.025))
	y__ = 1.25*sin(theta+(2*pi*0.025))

	## plot interactions
	plot(x=c(-1:1)*2,y=c(-1:1)*2,cex=0,bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
	for(i in 1:N)
	{
		for(j in 1:N)
		{
			color_ = if(effectsMat[i,j]>0){"green"}else{"red"}
			# points(x__[i],y__[i],cex=30/N)
			text(x__[i],y__[i],labels=labels[i])
			if(weightsMat[i,j]*10>0)
			{
				arrows(x0=x[j],x1=x_[i],y0=y[j],y1=y_[i],lwd=weightsMat[i,j]*10,col=color_,length=0.1)
			}
		}
	}
}

## .plot.diamondplot
##
## goal: plot distributionns of the variables in a table
##
## input:
## df - matrix - matrix containing variables to be plotted by columns
##
## output:
## none
.plot.diamond = function(df,y_lim=NULL,colVect=NULL,y_lab="")
{
	## initiate graphical parameters
	alpha = 0.2 # scaling factor for the plot
	n_col = ncol(df)
	n_row = nrow(df)
	x = 1:n_col
	x_lim = c(min(x)-alpha,max(x)+alpha)
	tmp = as.vector(apply(df,2,function(x)density(x)$x))
	if(is.null(y_lim)) y_lim = c(min(tmp),max(tmp));
	if(is.null(colVect)) colVect = rainbow(n_col,alpha=0.5) 

	## plot
	plot(rep(1,n_row),df[,1],xlim=x_lim,ylim=y_lim,cex=0,xlab="",ylab=y_lab)
	for(j in 1:n_col)
	{
		## coordindates
		x_ = rep(j,n_row)
		y_ = df[,j]
		x_density = density(df[,j])$x
		y_density = density(df[,j])$y
		y_density = y_density/max(y_density)*alpha # normalise

		## draw distributions
		polygon(x=c(j-y_density,rev(j+y_density)),y=c(x_density,rev(x_density)),col=colVect[j],border=NA)
	}
	lines(0:(n_col+1),rep(0,n_col+2),lty=2)
}

## anchorSampling
##
## goal: function to perform anchor sampling of a given posterior distribution
##
## input:
## n         - int - number of samples to take
## fn_prior  - fn  - function to sample from prior distribution of parameters 
## fn_argmax - fn  - function that optimises posterior distribution given parameters sampled from the priors
##
## output:
## none
anchorSample = function(n,fn_prior,fn_argmax)
{
	ensemble = NULL
	for(i in 1:n)
	{
		message(paste(i,"/",n,sep=""))
		Omega_0 = fn_prior()
		Omega_0 = fn_argmax(Omega_0)
		ensemble = rbind(ensemble,Omega_0)
	}
	return(ensemble)
}



#
###
