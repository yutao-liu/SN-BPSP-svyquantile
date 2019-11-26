###########################################################
# Title: Simulate survey population
#
# Author: Yutao Liu
#
# Date: 08/21/2017
#
# R version: 3.3.2
###########################################################

# if ( !require(pacman) ) install.packages(pacman)
# pacman::p_load(rstan, dplyr, loo, car)
pacman::p_load(sn)
pacman::p_load(pps)

## N = population size, n = PPS sample size
## distr = c('SkewN', 'Gamma', 'LogNorm')
## trans = c('none', 'log', 'sqrt')

SimSurvyPopu <- function (N, n, size.ty, distr, assoc, xi = NULL, alpha = NULL, sigma = NULL, gamma = NULL, trans = NULL,
                          doFigure = FALSE, seed) {
  
  set.seed(seed)
  
  if (size.ty == 1) size = rgamma(n = N, shape = 1.5, rate = .001);
  if (size.ty == 2) size = rgamma(n = N, shape = 5, rate = .0025);
  
  # ------- probability of selection -----------------------
  std.size = (size/sum(size))*n  # standardized size
  
  if (distr == "SkewN") {
    
    if (assoc == FALSE) {
      if ( is.null(xi) ) {
        xis <- 150 + mean(std.size)
      } else {
        xis = xi
      }
    } else {
      xis <- 150 + 100*std.size
    }
    omegas = sqrt( (alpha^2 + 1) * sigma^2 * std.size^(2*gamma) )
    y = rsn(n = N, xi = xis, omega = omegas)
  }
  
  # if (distr == "MixN") {
  #   
  #   n1 = ceiling(.9 * N)
  #   y1 = sapply(1:n1, function (i) {
  #     return( rnorm(n = 1, mean = xis[i], sd = sigma) )
  #   }) 
  #   
  #   y2 = sapply( (n1 + 1):N, function (i) {
  #     return( rnorm(n = 1, mean = xis[i] + discr * sigma, sd = 2 * sigma) )
  #   })
  #   
  #   y = c(y1, y2)
  #   
  #   y1 = y1 - min(y) + 5; y2 = y2 -min(y) + 5;   yplot = y - min(y) + 5   # avoid negative observations
  #   
  #   if (doFigure) {
  #     
  #     pdf(file = paste(distr, size.ty, assoc, ".pdf", sep = "_"), height = 6, width = 6.5)
  #     par( mfrow = c(1, 1) )
  #     plot(density(yplot), main = "Density of survey variable", col = "red", ylim = c( 0, max(density(yplot)$y ) * 1.3))
  #     lines(density(y1))
  #     lines(density(y2)); rug(y2)
  #     graphics.off()  
  #   }
  # }
  
  if (distr == "Gamma") {
    
    if ( assoc == FALSE) {
      shape <- 0.3*log(size)
      scale <- 5e2/log(size)
    } else {
      shape = 0.3*log(size) 
      scale = sqrt(size)
    }
  
    y <- rgamma(n = N, shape = shape, scale = scale)
  }
  
  if (distr == "LogNorm") {
    
    if ( assoc == FALSE ) {
      mu <- -5 + mean(std.size)
    } else {
      mu <- -2 + 5*std.size + log(1e3)
    }
    
    y <- rlnorm(n = N, meanlog = mu, sdlog = sigma^2 * std.size^(2*gamma))
  }
  
  if ( !is.null(trans) ) {
    if( any( !trans %in% c('none', 'log', 'sqrt') ) ) {
      stop("The transformation specified is not available!")
    }
  
    if (trans[1] == 'log') std.size <- exp(std.size)
    if (trans[1] == 'sqrt') std.size <- std.size^2
    
    if (trans[2] == 'log') y <- exp(y)
    if (trans[2] == 'sqrt') y <- y^2
  }
  
  if ( min(y) < 0 ) {
    cat('Location shift to avoid mininal negative value.')
    y = y - min(y) + 5   # avoid negative observations
  }
  
  probs <- std.size/sum(std.size)*n
  
  popu = data.frame(sub = 1:N, std.size, probs, y) 
  
  return(popu)
  
}


# ------- PPS sampling ----------------------------------

ppss.smpl <- function (n, popu, seed = 2017) {
  set.seed(seed)
  N <- nrow(popu)
  if (!( 'probs' %in% names(popu) ) ) {
    popu <- mutate(popu, probs = std.size/sum(std.size)*n)
  }
  permutation = sample(1:N)
  current = popu[permutation, ]  # permute population
  ppss.sample = ppss(current$probs, n)
  smpl = current[ppss.sample, ]  
  nonsmpl = current[-ppss.sample, ]
  return( list(smpl = smpl, nonsmpl = nonsmpl) )
}



