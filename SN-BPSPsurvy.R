###########################################################
# Title: Skew-Normal Bayesian penalized spline model-based 
#        inference with survey data 
#
# Author: Yutao Liu
#
# R version: 3.3.2
###########################################################

pacman::p_load(rstan, loo, car)
pacman::p_load(grid)
pacman::p_load(gridExtra)

# ------- Transform input data -------------------------

## transfm[1] = transformation for probability of selection
## transfm[2] = transformation for survey variable

transFn <- function(inputdat, transfm) {
  
  if( !transfm[1] %in% c('none', 'log', 'sqrt', 'logit') ) {
    stop("The transformation specified for probability of selection is not available!")
  }
  
  if( !transfm[2] %in% c('none', 'log', 'sqrt') ) {
    stop("The transformation specified for survey variable is not available!")
  }
  
  if (transfm[1] == 'log') inputdat$probs <- inputdat$probs %>% log(.)
  
  if (transfm[1] == 'sqrt') inputdat$probs <- inputdat$probs %>% sqrt(.)
  
  if (transfm[1] == 'logit') inputdat$probs <- inputdat$probs %>% logit(.)

  if (transfm[2] == 'log') inputdat$y <- inputdat$y %>% log(.)
  
  if (transfm[2] == 'sqrt') inputdat$y <- inputdat$y %>% sqrt(.)
  
  return( inputdat )
}

# ------- Bayesian inference with survey data ----------

SNsurvyBayes <- function (dat, transfm = c('none', 'none'), mod, smpl.n = 50, num.knots = 15, n.iter = 1e3, n.chains = 3, seed = c(1, 1), percentiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), doSample = TRUE, doFigure = FALSE) {
  
  print(smpl.n)
  # ------- probability of selection ----------------
  if (!( 'probs' %in% names(dat) ) ) {
    dat <- mutate(dat, probs = size/sum(size)*smpl.n)
  }
  
  # ------- sample from population ------------------
  # set.seed(seed[1])
  if (doSample) {
    origdat <- ppss.smpl (n = smpl.n, popu = dat, seed = seed[1])
    
  } else {
    origdat <- dat
  }
  
  Pi <- lapply(origdat, function (inputdat) {
    return (inputdat$probs)
  } )
  
  moddat <- lapply(origdat, transFn, transfm)
  
  if (doSample) {
    plotdat <- moddat$smpl
  } else {
    plotdat <- moddat
  }
  
  # ------- plot sample ------------------
  
  if (doFigure & doSample) {
    
    # pdf(file = paste(distr, size.ty, assoc, transfm[1], transfm[2], "smpl.pdf", sep = "_"), height = 4.5, width = 8)
    # 
    # p1 <- ggplot(data = dat, aes(probs) ) + geom_histogram(aes(y=..density..), bins = 100) + geom_density(alpha = .2, fill="#FF6666") + labs(title = "Distribution of Probability of Selection", x = "Probability of Selection") 
    # p2 <- ggplot(dat, aes(probs, y)) + geom_point(colour = "black", alpha = .3) + geom_point(data = origdat$smpl, colour = "red") + labs(title = "Scatter Plot of Survey Pupulation", x = "Probability of Selection", y = 'Survey Variable')
    # p3 <- ggplot(data = plotdat, aes(probs) ) + geom_histogram(aes(y=..density..), bins = 20) + geom_density(alpha = .2, fill="#FF6666") + labs(title = "Distribution of (Transformed) Probability of Selection", x = "(Transformed) Probability of Selection") 
    # p4 <- ggplot(plotdat, aes(probs, y) ) + geom_point(colour = "red") + labs(title = 'PPS Sample', x = '(Transformed) Probability of Selection', y = '(Transformed) Survey Variable');
    # 
    # grid.arrange(p1, p2, p3, p4, ncol = 2)
    # graphics.off()
    
    pdf(file = paste(distr, size.ty, assoc, transfm[1], transfm[2], "smpl.pdf", sep = "_"), height = 3, width = 8)
    p1 <- ggplot(dat, aes(probs, y) ) + geom_point(colour = "black", alpha = .3) +
           geom_point(data = origdat$smpl, colour = "red") +
           labs(title = "Survey Pupulation", x = "Probability of Selection", y = 'Survey Variable')
    p2 <- transFn(dat, transfm) %>% ggplot(., aes(probs, y) ) +
           geom_point(colour = "black", alpha = .3) +
           geom_point(data = plotdat, aes(probs, y), colour = "red") +
           labs(title = 'Survey Pupulation', x = '(Transformed) Probability of Selection', y = '(Transformed) Survey Variable');
    grid.arrange(p1, p2, ncol = 2)
    graphics.off()
    
  }
  
  summary(moddat$smpl) %>% print(.)
  
  if (doSample) {
    
    # ----- design matrices: polynomial and spline -------------
    knots = with(moddat$smpl, quantile(probs, probs = seq(0, 1, length = (num.knots+2))[ -c(1, (num.knots+2) ) ] ) )
    
    X = with(moddat$smpl, cbind(1, probs) )
    Z = with(moddat$smpl, outer(probs, knots, "-") )  # truncated 
    Z = Z * (Z > 0)
    
    X.pred = with(moddat$nonsmpl, cbind(1, probs) )
    Z.pred = with(moddat$nonsmpl, outer(probs, knots, "-") )  # truncated 
    Z.pred = Z.pred * (Z.pred > 0)
    
    n.nonsmpl <- nrow(moddat$nonsmpl)
    
    if (mod == 'SN-B2PSP') {
      moddat.stan <- list(n1 = nrow(moddat$smpl), n2 = n.nonsmpl, p = ncol(X), K = num.knots, y = moddat$smpl$y, X = X, Z = Z, predX = X.pred, predZ = Z.pred)
      fit <- stan(file = "SN-B2PSP-prediction.stan", data = moddat.stan, iter = 100, chains = n.chains, seed = seed[2])
      fit <- stan(fit = fit, data = moddat.stan, iter = n.iter, chains = n.chains)
      print(fit, c('beta', 'b', 'lambda', 'nu', 'alpha', 'lp__') )
    }
    
    if (mod == 'SN-B2PSP-cauchy') {
      moddat.stan <- list(n1 = nrow(moddat$smpl), n2 = n.nonsmpl, p = ncol(X), K = num.knots, y = moddat$smpl$y, X = X, Z = Z, predX = X.pred, predZ = Z.pred)
      fit <- stan(file = "SN-B2PSP-cauchy.stan", data = moddat.stan, iter = 100, chains = n.chains, seed = seed[2])
      fit <- stan(fit = fit, data = moddat.stan, iter = n.iter, chains = n.chains)
      print(fit, c('beta', 'b', 'lambda', 'nu', 'alpha', 'lp__') )
    }
    
	if (mod == 'SN-B2PSP-uniform') {
      moddat.stan <- list(n1 = nrow(moddat$smpl), n2 = n.nonsmpl, p = ncol(X), K = num.knots, y = moddat$smpl$y, X = X, Z = Z, predX = X.pred, predZ = Z.pred)
      fit <- stan(file = "SN-B2PSP-uniform.stan", data = moddat.stan, iter = 100, chains = n.chains, seed = seed[2])
      fit <- stan(fit = fit, data = moddat.stan, iter = n.iter, chains = n.chains)
      print(fit, c('beta', 'b', 'lambda', 'nu', 'alpha', 'lp__') )
    }
	
    if (mod == 'SN-B2PSP-omega') {
      moddat.stan <- list(n1 = nrow(moddat$smpl), n2 = n.nonsmpl, p = ncol(X), K = num.knots, y = moddat$smpl$y, X = X, Z = Z, predX = X.pred, predZ = Z.pred)
      fit <- stan(file = "SN-B2PSP-omega.stan", data = moddat.stan, iter = 100, chains = n.chains, seed = seed[2])
      fit <- stan(fit = fit, data = moddat.stan, iter = n.iter, chains = n.chains)
      print(fit, c('beta', 'b', 'lambda', 'nu', 'alpha', 'lp__') )
    }
    
    if (mod == 'SN-poly') {
      moddat.stan <- list(n1 = nrow(moddat$smpl), n2 = n.nonsmpl, p = ncol(X), K = num.knots, y = moddat$smpl$y, X = X, Z = Z, predX = X.pred, predZ = Z.pred, Pi = origdat$smpl$probs, predPi = origdat$nonsmpl$probs)
      fit <- stan(file = "SN-poly-prediction.stan", data = moddat.stan, iter = 50, chains = n.chains, seed = seed[2])
      fit <- stan(fit = fit, data = moddat.stan, iter = n.iter, chains = n.chains)
      print(fit, c('beta', 'b', 'sigma2', 'gamma', 'alpha', 'lp__') )
    }
    
    if (mod == 'SN-poly-cauchy') {
      moddat.stan <- list(n1 = nrow(moddat$smpl), n2 = n.nonsmpl, p = ncol(X), K = num.knots, y = moddat$smpl$y, X = X, Z = Z, predX = X.pred, predZ = Z.pred, Pi = origdat$smpl$probs, predPi = origdat$nonsmpl$probs)
      fit <- stan(file = "SN-poly-cauchy.stan", data = moddat.stan, iter = 50, chains = n.chains, seed = seed[2])
      fit <- stan(fit = fit, data = moddat.stan, iter = n.iter, chains = n.chains)
      print(fit, c('beta', 'b', 'sigma', 'gamma', 'alpha', 'lp__') )
    }
	
	if (mod == 'SN-poly-uniform') {
      moddat.stan <- list(n1 = nrow(moddat$smpl), n2 = n.nonsmpl, p = ncol(X), K = num.knots, y = moddat$smpl$y, X = X, Z = Z, predX = X.pred, predZ = Z.pred, Pi = origdat$smpl$probs, predPi = origdat$nonsmpl$probs)
      fit <- stan(file = "SN-poly-uniform.stan", data = moddat.stan, iter = 50, chains = n.chains, seed = seed[2])
      fit <- stan(fit = fit, data = moddat.stan, iter = n.iter, chains = n.chains)
      print(fit, c('beta', 'b', 'sigma', 'gamma', 'alpha', 'lp__') )
    }
	
	if (mod == 'SN-poly-omega') {
      moddat.stan <- list(n1 = nrow(moddat$smpl), n2 = n.nonsmpl, p = ncol(X), K = num.knots, y = moddat$smpl$y, X = X, Z = Z, predX = X.pred, predZ = Z.pred, Pi = origdat$smpl$probs, predPi = origdat$nonsmpl$probs)
      fit <- stan(file = "SN-poly-omega.stan", data = moddat.stan, iter = 50, chains = n.chains, seed = seed[2])
      fit <- stan(fit = fit, data = moddat.stan, iter = n.iter, chains = n.chains)
      print(fit, c('beta', 'b', 'sigma', 'gamma', 'alpha', 'lp__') )
    }
    # print(fit, c('predxi', 'predomega') )
    
    xi.draws <- extract(fit, "xi", permuted = TRUE)[[1]]
    omega.draws <- extract(fit, "omega", permuted = TRUE)[[1]]
    alpha.draws <- extract(fit, "alpha", permuted = TRUE)[[1]]
    predxi.draws <- extract(fit, "predxi", permuted = TRUE)[[1]]
    predomega.draws <- extract(fit, "predomega", permuted = TRUE)[[1]]
    
    # yhats <- extract(fit, "yhat", permuted = TRUE)[[1]]
    
    chain.length <- length(alpha.draws)
    
    # ------- posterior predictive check ---------
    
    yrep <- sapply(1:chain.length, function (t) {
      yrep.iter <- rsn(n = smpl.n, xi = xi.draws[t, ], omega = omega.draws[t, ], alpha = alpha.draws[t])
      return( yrep.iter )
    } )
    
    # if ( transfm[2] == 'none' ) {
    #   yrep <- yrep * (yrep > 0)
    # } else {
    #   
    #   if (transfm[2] == 'log') {
    #     yrep <- exp(yrep)
    #   } 
    #   
    #   if (transfm[2] == 'sqrt') {
    #     yrep <- ( yrep * (yrep > 0) )^2
    #   }
    # }
    
    Trepmat <- ( ( t(yrep) - xi.draws ) / omega.draws )^2
    Trep <- apply(Trepmat, 1, mean)
    
    Tmat <- sweep(xi.draws, 2, moddat$smpl$y)
    Tmat <- ( Tmat/omega.draws )^2
    Tobs <- apply(Tmat, 1, mean)
    
    yrepSD <- apply(yrep, 2, stats::sd)
    yobsSD <- stats::sd(moddat$smpl$y)
    
    pval1 <- mean(yrepSD > yobsSD)
    pval2 <- mean(Trep > Tobs)
    
    # sum(Tmat[2, ] == xi.draws[2, ] - moddat$smpl$y); sum(Tmat[1, ] != xi.draws[1, ] - moddat$smpl$y)
    # cat('posterior predictive p-value with fitted sample:')
    # mean(Trep > Tobs) %>% print(.)
    
    if (doFigure) {
      pdf(file = paste(distr, size.ty, assoc, transfm[1], transfm[2], mod, "ppp_obs.pdf", sep = "_"), height = 7, width = 14)
      F1 <- data.frame(sd = yrepSD) %>% ggplot(., aes(sd) ) + geom_histogram(bins = 50) +
        labs(x = 'standard deviation of y-rep', subtitle = paste('p = ', round(pval1, 2), sep = '')) +
        geom_vline(xintercept = yobsSD, color = 'red')
      F2 <- data.frame(Trep, Tobs) %>% ggplot(., aes(x = Tobs, y = Trep)) + geom_point(colour = "black", alpha = .3) +
        geom_abline(intercept = 0, slope = 1, color = 'red') +
        labs(x = 'T (y, theta)', y = 'T (y-rep, theta)', subtitle = paste('p = ', round(pval2, 2), sep = ''))
      gridExtra::grid.arrange(F1, F2, ncol = 2)
      graphics.off()
    }
    
    # ------- prediction -------------------------
    
    yhats <- sapply(1:chain.length, function (t) {
      yhats.iter <- rsn(n = n.nonsmpl, xi = predxi.draws[t, ], omega = predomega.draws[t, ], alpha = alpha.draws[t])
      return( yhats.iter )
    } )
    
    Thatmat <- ( ( t(yhats) - predxi.draws ) / predomega.draws )^2
    That <- apply(Thatmat, 1, mean)
    
    Tmatunobs <- sweep(predxi.draws, 2, moddat$nonsmpl$y)
    # sum(Tmatunobs[2, ] == predxi.draws[2, ] - moddat$nonsmpl$y); sum(Tmatunobs[1, ] != predxi.draws[1, ] - moddat$nonsmpl$y)
    Tmatunobs <- ( Tmatunobs/predomega.draws )^2
    Tunobs <- apply(Tmatunobs, 1, mean)
    
    cat('posterior predictive p-value with non-sample units:')
    mean(That > Tunobs) %>% print(.)
    
    if (doFigure) {
      pdf(file = paste(distr, size.ty, assoc, transfm[1], transfm[2], mod, "ppp_unobs.pdf", sep = "_"), height = 7, width = 7)
      p <- data.frame(That, Tunobs) %>% ggplot(., aes(Tunobs, That)) + geom_point(colour = "black", alpha = .3) +
        geom_abline(intercept = 0, slope = 1, color = 'red') + coord_equal(ratio = 1)
      print(p)
      graphics.off()
    }
    
    if ( transfm[2] == 'none' ) {
      yhats <- yhats * (yhats > 0)
    } else {
      
      if (transfm[2] == 'log') {
        yhats <- exp(yhats)
      } 
      
      if (transfm[2] == 'sqrt') {
        yhats <- ( yhats * (yhats > 0) )^2
      }
    }
    
    # if (doFigure) {
    #   
    #   # pdf(file = paste(distr, size.ty, assoc, transfm[1], transfm[2], mod, "prediction.pdf", sep = "_"), height = 6, width = 13)
    #   # par( mfrow = c(1, 2) )
    #   # 
    #   # plot( density(origdat$nonsmpl$y), main = "Non-sampled and Prediction -- y", ylim = c( 0, max(density(origdat$nonsmpl$y)$y ) * 1.5), lwd = 2.5 )
    #   # for (i in chain.length:(chain.length - 10) ) {
    #   #   lines(density(yhats[ , i]), col = "#FF0000CC")
    #   # }
    #   # 
    #   # with( origdat$nonsmpl, plot( probs, y, pch = 19, main = "Non-sampled and Prediction", ylim =  c( min(y) - 0.2*( max(y) - min(y) ), max(y) + 0.2*( max(y) - min(y) ) )) )
    #   # with( origdat$nonsmpl, points(probs, yhats[ , chain.length - 10], pch = 19, col = "#FF0000CC") )
    #   # graphics.off()
    # 
    #   # quantlvl <- runif(n.nonsmpl, min = 0, max = 1) %>% sort(.)
    #   # 
    #   # figurequantest <- sapply(1:chain.length, function(i) {
    #   #   return( quantile( yhats[ , i], quantlvl ) )
    #   # })
    #   # 
    #   # figurequantpred <- apply(figurequantest, 1, median)
    #   # figurequantobs <- quantile(origdat$nonsmpl$y, quantlvl)
    #   # n.nonsmpl.quar <- floor(n.nonsmpl/4)
    #   # 
    #   # pdf(file = paste(distr, size.ty, assoc, transfm[1], transfm[2], mod, "prediction.pdf", sep = "_"), height = 8, width = 8)
    #   # p1 <- data.frame(pred = figurequantpred, obs = figurequantobs)[1:n.nonsmpl.quar, ] %>%
    #   #        ggplot(., aes(obs, pred)) + geom_point(colour = "black", alpha = .3) +
    #   #        geom_abline(intercept = 0, slope = 1, color = 'red') + coord_equal(ratio = 1) +
    #   #        labs(title = 'Quantile-Quantile Plot', x = 'Non-Sampled', y = 'Predicted')
    #   # p2 <- data.frame(pred = figurequantpred, obs = figurequantobs)[(n.nonsmpl.quar+1):(2*n.nonsmpl.quar), ] %>%
    #   #        ggplot(., aes(obs, pred)) + geom_point(colour = "black", alpha = .3) +
    #   #        geom_abline(intercept = 0, slope = 1, color = 'red') + coord_equal(ratio = 1) +
    #   #        labs(title = 'Quantile-Quantile Plot', x = 'Non-Sampled', y = 'Predicted')
    #   # p3 <- data.frame(pred = figurequantpred, obs = figurequantobs)[(2*n.nonsmpl.quar+1):(3*n.nonsmpl.quar), ] %>%
    #   #        ggplot(., aes(obs, pred)) + geom_point(colour = "black", alpha = .3) +
    #   #        geom_abline(intercept = 0, slope = 1, color = 'red') + coord_equal(ratio = 1) +
    #   #        labs(title = 'Quantile-Quantile Plot', x = 'Non-Sampled', y = 'Predicted')
    #   # p4 <- data.frame(pred = figurequantpred, obs = figurequantobs)[(3*n.nonsmpl.quar+1):n.nonsmpl, ] %>%
    #   #        ggplot(., aes(obs, pred)) + geom_point(colour = "black", alpha = .3) +
    #   #        geom_abline(intercept = 0, slope = 1, color = 'red') + coord_equal(ratio = 1) +
    #   #        labs(title = 'Quantile-Quantile Plot', x = 'Non-Sampled', y = 'Predicted')
    #   # grid.arrange(p1, p2, p3, p4, ncol = 2)
    #   # graphics.off()
    #   
    #   # nonsmpl.sort <- sort(origdat$nonsmpl$y)
    #   # # pred.sort <- apply(yhats, 1, median) %>% sort(.)
    #   # pred.sort <- apply(yhats[ , 1:10], 2, sort)
    #   # 
    #   # sortplot <- data.frame(nonsmpl.sort, pred.sort)
    #   # 
    #   # trm <- c(.01, .99)*n.nonsmpl %>% as.integer(.)
    #   # sortplot <- sortplot[ trm[1]:trm[2], ]
    #       
    #   # pdf(file = paste(distr, size.ty, assoc, transfm[1], transfm[2], mod, "prediction.pdf", sep = "_"), height = 8, width = 8)
    #   # p <- ggplot(sortplot, aes(nonsmpl.sort, X1) ) +
    #   #   geom_point(colour = "black", alpha = .3) + geom_line(aes(nonsmpl.sort, X1 ), colour = "red", alpha = .3) + 
    #   #   geom_point(aes(nonsmpl.sort, X2 ), colour = "black", alpha = .3) + geom_line(aes(nonsmpl.sort, X2 ), colour = "red", alpha = .3) + 
    #   #   geom_point(aes(nonsmpl.sort, X3 ), colour = "black", alpha = .3) + geom_line(aes(nonsmpl.sort, X3 ), colour = "red", alpha = .3) + 
    #   #   geom_point(aes(nonsmpl.sort, X4 ), colour = "black", alpha = .3) + geom_line(aes(nonsmpl.sort, X4 ), colour = "red", alpha = .3) + 
    #   #   geom_point(aes(nonsmpl.sort, X5 ), colour = "black", alpha = .3) + geom_line(aes(nonsmpl.sort, X5 ), colour = "red", alpha = .3) + 
    #   #   geom_point(aes(nonsmpl.sort, X6 ), colour = "black", alpha = .3) + geom_line(aes(nonsmpl.sort, X6 ), colour = "red", alpha = .3) + 
    #   #   geom_point(aes(nonsmpl.sort, X7 ), colour = "black", alpha = .3) + geom_line(aes(nonsmpl.sort, X7 ), colour = "red", alpha = .3) + 
    #   #   geom_point(aes(nonsmpl.sort, X8 ), colour = "black", alpha = .3) + geom_line(aes(nonsmpl.sort, X8 ), colour = "red", alpha = .3) + 
    #   #   geom_point(aes(nonsmpl.sort, X9 ), colour = "black", alpha = .3) + geom_line(aes(nonsmpl.sort, X9 ), colour = "red", alpha = .3) + 
    #   #   geom_point(aes(nonsmpl.sort, X10), colour = "black", alpha = .3) + geom_line(aes(nonsmpl.sort, X10), colour = "red", alpha = .3) +
    #   #   labs(title = 'Quantile-Quantile Plot', x = 'Non-Sampled', y = 'Predicted') +
    #   #   coord_equal(ratio = 1, xlim = c(min(sortplot), max(sortplot)), ylim = c(min(sortplot), max(sortplot)))
    #   # print(p)
    #   # graphics.off()
    # }
    
    est.percentiles <- sapply(1:chain.length, function(i) {
      return( quantile( c(origdat$smpl$y, yhats[ , i]), percentiles ) )
    })
    
    res.percentiles = rbind( apply(est.percentiles, 1, median), apply(est.percentiles, 1, quantile, c(.025, .975)) )
    
    est.total <- sapply(1:chain.length, function(i) {
      return( sum( c(origdat$smpl$y, yhats[ , i]) ) )
    })
    
    res.total = c( mean = median(est.total), quantile(est.total, c(.025, .975) ) )
    
    res = rbind( res.total, t(res.percentiles) ); row.names(res)[1] = "total"; colnames(res) = c("est", "CI 2.5%", "CI 97.5%")
  }
  
  return (res)
  
}
