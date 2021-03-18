rm(list = ls())
graphics.off()
#R versions <4.0.0 convert strings to factors, specify default behaviour
formals(data.frame)$stringsAsFactors <- FALSE
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster            DATE: 2019 06 07
#     MODIFIED:	James Foster            DATE: 2021 03 18
#
#  DESCRIPTION: Fits a maximum-likelihood quadrimodal distribution by constructing
#               the probability density function in the same way as a bimodal
#               "mixed von Mises" distribution. The log likelihood is then compared
#               with a circular uniform distribution. The likelihood is compared
#               using a likelihood ratio test, for which the degrees of freedom 
#               are the minimum number of parameters to describe a quadrimodal   
#               von Mises distribution. The test is demonstrated here with
#               simulated data.
#               
#               
#      OUTPUTS: Plots.
#
#	   CHANGES: -description
#             -
#
#   REFERENCES: Fitak R. & Johnsen S. (2017).
#               Bringing the analysis of animal orientation data full circle:
#               model-based approaches with maximum likelihood.
#               Journal of Experimental Biology 220: 3878-3882
#               https://doi.org/10.1242/jeb.167056
#
#               Schnute J. T. & Groot K. (1992) 
#               Statistical analysis of animal orientation data.
#               Animal Behaviour 43, 15–33 (1992).
#               
#               Oliveriusová, L., Němec, P., Pavelková, Z. et al. (2014)
#               Spontaneous expression of magnetic compass orientation in an 
#               epigeic rodent: the bank vole, Clethrionomys glareolus. 
#               Naturwissenschaften 101, 557–563.
#               https://doi.org/10.1007/s00114-014-1192-0
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Comment  
#- Test run  
#- Plots          


# Load relevant packages --------------------------------------------------
require(circular)
require(CircMLE)


# Set up probability density functions ------------------------------------

# for reference, this is the density function for a "mixed" von Mises
# (two distributions combined, e.g. bimodal) from the "circular" package

# DmixedvonmisesRad <- function(x, mu1, mu2, kappa1, kappa2, prop) {
  # vm <- prop/(2 * pi * besselI(x=kappa1, nu=0, expon.scaled = TRUE)) * (exp(cos(x - mu1) - 1))^kappa1 + (1 - prop)/(2 * pi * besselI(x=kappa2, nu=0, expon.scaled = TRUE)) * (exp(cos(x - mu2) - 1))^kappa2
  # return(vm)
# }

# this is the function from circular for extracting random draws from that function
# N.B. Draws are extracted from two separate von Mises distributions, and combined
# using the proportion defined by the user

# RmixedvonmisesRad <- function(n, mu1, mu2, kappa1, kappa2, prop) {
  # result <- rep(NA, n)
  # test <- runif(n)
  # n1 <- sum(test < prop)
  # n2 <- n - n1
  # res1 <- RvonmisesRad(n1, mu1, kappa1)
  # res2 <- RvonmisesRad(n2, mu2, kappa2)
  # result[test < prop] <- res1
  # result[test >= prop] <- res2
  # return(result)
# }

# Set up the random quadrimodal function in the same way
# N.B. for simplicity it is assumed that:
# - there are equal numbers of samples from the distributions around each mean
# (i.e. equal proportions, 1/4 of the distribution)
# - each mean has the same mean vector length as all others
# (i.e. kappa1 = kappa2 = kappa3 = kappa4)
# - the means are equally-spaced (i.e. separated by exactly 90°)
RquadvonmisesRad <- function(n, mu, kappa) {
  result <- rep(NA, n)#pre-assign the output vector
  test <- runif(n)#make a uniform distribution (0,1) of the correct length
  n1 <- sum(test < 0.25)#mean 1 has as many entries as the bottom quarter
  n2 <- sum(test >= 0.25 & test < 0.5 )#mean 2 has as many entries as the 2nd quarter
  n3 <- sum(test >= 0.5 & test < 0.75 )#mean 3 has as many entries as the 3rd quarter
  n4 <- sum(test >= 0.75)#mean 4 has as many entries as the 4th quarter
  #create entries for each distribution
  res1 <- rvonmises(n1, mu, kappa)
  res2 <- rvonmises(n2, mu + pi/2, kappa)
  res3 <- rvonmises(n3, mu + pi, kappa)
  res4 <- rvonmises(n4, mu + 3*pi/2, kappa)
  #assign entries for each distribution based on the quarters of the random uniform
  result[test < 0.25] <- res1
  result[test >= 0.25 & test < 0.5 ] <- res2
  result[test >= 0.5 & test < 0.75 ] <- res3
  result[test >= 0.75] <- res4
  return(circular(result))#output in circular format
}

# Set up the probability density function, sharing total available probability (100%)
# between each of the means equally (25% each).
# At a given angle, probability density is the sum of the probability of an observation
# from the distribution centred on _any_ of the four means.
DquadvonmisesRad <- function(x, mu, kappa) {
  vm <- 0.25/(2 * pi * besselI(x=kappa, nu=0, expon.scaled = TRUE)) * (exp(cos(x - mu) - 1))^kappa +
	   0.25/(2 * pi * besselI(x=kappa, nu=0, expon.scaled = TRUE)) * (exp(cos(x - (mu+pi/2)) - 1))^kappa +
	   0.25/(2 * pi * besselI(x=kappa, nu=0, expon.scaled = TRUE)) * (exp(cos(x - (mu+pi)) - 1))^kappa +
	   0.25/(2 * pi * besselI(x=kappa, nu=0, expon.scaled = TRUE)) * (exp(cos(x - (mu+3*pi/2)) - 1))^kappa
  return(vm)
}

# Fit the quadrimodal model using the same methods Fitak & Johnsen use in CircMLE
# to fit bimodal (or other mixed) von Mises distributions
QuadMod<- function(data, BadStart, nchains, method, niter){
    if (missing(BadStart)) 
        BadStart = 10^9
    else BadStart = BadStart
    if (BadStart < 0) 
        stop("The value for starting parameters outside the preset limits must be >0")
    if (missing(nchains)) 
        nchains = 5
    else nchains = nchains
    if (nchains < 1) 
        stop("Must set the number of chains to an integer >=1")
    if (missing(niter)) 
        niter = 5000
    else niter = niter
    if (niter < 1000) 
        warning("At least 1000 iterations is recommended but not required. Check ?optim for details.")
    if (missing(method)) 
        method = "BFGS"
    else method = method
    if (method != "Nelder-Mead" & BadStart == Inf) 
        stop("Except for Nelder-Mead, all other optimization algorithms require finite starting parameters")
    # Set up a log likelihood function for the optimiser
    mQuad = function(params) {
        if (params[1] < 0 | params[1] > 2 * pi | params[2] <= 
            0 | params[2] > 227) {
            R = BadStart
            return(R)
        }
        else {
            R = DquadvonmisesRad(data, mu = circular(params[1]),  
                kappa = params[2])
            R = -sum(log(R))#negative log likelihood is output
            # this quantity will be minimised by the optimiser
            # minimum -log likelihood and maximum likelihood occur at the same
            # parameter values
            return(R)
        }
    }
    # start looking for the mean at random angles
    q1 = as.numeric(rcircularuniform(nchains, control.circular = list(modulo = "2pi")))
    # start looking for kappa somewhere between 1 and 5
    k1 = as.numeric(sample(1:5, nchains, replace = T))
    mQuad.out = list()
    #run the optimiser for each chain
    for (i in 1:nchains) {
        chain.out = suppressWarnings(stats::optim(c(q1[i], k1[i]), 
            fn = mQuad, method = method, control = list(maxit = niter)))
        names(chain.out)[2] = "lik"
        mQuad.out[[i]] = chain.out
    }
    #choose the parameters from the chain with the maximum likelihood
    min = which.min(sapply(mQuad.out, function(x) x[2]))
    return(mQuad.out[[min]])
}

# Simulate a quadrimodal distribution -------------------------------------

xx <- RquadvonmisesRad(n = 40,#reasonable sample size
                       mu = circular(pi/3),#mean at 60°
                       kappa = exp(4)#mean vector length of almost 1
                       )
xx <- (xx-pi)*180/pi #now in degrees
#Fit a quadrimodal model
qm <- QuadMod(xx*pi/180)#now in radians


# Plot the data -----------------------------------------------------------
# raw data
plot(circular(xx, type = 'angles', units = 'degrees', zero = pi/2, rotation = 'clock'),
	stack = T, sep = 0.03, shrink = 1.5)
# Fitted mean vectors
arrows.circular(circular(-qm$par[1]+c(0,pi/2,pi,3*pi/2)),
			shrink = A1(qm$par[2]),
			col = 'darkred', len = 0.1)


# Test against a uniform distribution -------------------------------------


uf <- M1(xx)#uniform
un <- M2A(xx)#unimodal
ax <- M3A(xx)#axial (equal variances)
npar <- 2#this symmetric quad has 2 parameters, primary mean and concentration
#calculate Akaike Information Criterion,
# assuming we know that, if oriented, the distribution would be:
# - quadrimodal
# - symmetric
# - equally spaced
qm.aic <- 2*npar + 2* qm$lik
#tell the user about the outcome of the test
# Unless there is a lot of very concentrated data,
# either the uniform distribution
# or the axial distribution has lower AIC 
if(2*uf$lik > qm.aic & ax$lik > qm$lik & un$lik > qm$lik )
{
  message('Quadrimodal!')
  message('primary mean = ', round(qm$par[1]*180/pi,1), '°')
  message('mean vector length = ', round(A1(qm$par[2]),5))
}else
  {
  if(ax$lik < qm$lik)
    {
    message('Axial')
    message('primary mean = ', round(ax$par[1]*180/pi,1), '°')
    message('mean vector length = ', round(A1(ax$par[2]),5))
  }else
  {
    if(un$lik < qm$lik)
    {
      message('Unimodal')
      message('mean angle = ', round(un$par[1]*180/pi,1), '°')
      message('mean vector length = ', round(A1(un$par[2]),5))
    }else
   {message('Uniform')}
  }
}
message('AIC uniform = ', round(2*uf$lik,2))
message('AIC unimodal = ', round(2*un$lik+2*length(ax$par),2))
message('AIC axial = ', round(2*ax$lik+2*length(ax$par),2))
message('AIC quadrimodal = ', round(qm.aic,2))
