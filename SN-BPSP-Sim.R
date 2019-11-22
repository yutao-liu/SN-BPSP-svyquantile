###########################################################
# Title: Skew-Normal Bayesian penalized spline model-based 
#        inference with artifical survey data 
#
# Author: Yutao Liu
#
# R version: 3.3.2
###########################################################

rm( list = ls() )

start.time <- Sys.time()

if ( !require(pacman) ) install.packages(pacman)
pacman::p_load(dplyr)
# pacman::p_load(pps)
pacman::p_load(ggplot2)
pacman::p_load(grid)
pacman::p_load(gridExtra)

source('SimSurvyPopu.R')
source('SN-BPSPsurvy.R')

# ------- Global Options -------------------------------

plot.popu <- TRUE # plot.popu <- FALSE
multitask <- TRUE # multitask <- FALSE
taskid <- 1
nrep <- 10

## ------- Simulation Setting: Survey Population -------

artfcl = FALSE
size.ty = 1
distr = 'ndatss' # distr = 'farm'
assoc = TRUE
trans = c('none', 'none')

# readdir <- 'C:/Dropbox/PhD@Columbia/Dissertation/Project 1 SN-BPSP/Code'

if (artfcl) {
  
  smpl.n <- 200
  
  if (distr == 'SkewN') {
    
    popu <- SimSurvyPopu (N = 2e3, n = 2e2, size.ty, distr, assoc, alpha = 4, sigma = 12, gamma = .8, trans, seed = 1)
    
  } else if (distr == 'Gamma') {
    
    popu <- SimSurvyPopu (N = 2e3, n = 2e2, size.ty, distr, assoc, trans, seed = 1)
    
  } else if (distr == 'LogNorm') {
    
    popu <- SimSurvyPopu (N = 2e3, n = 2e2, size.ty, distr, assoc, sigma = .8, gamma = 0.1, trans, seed = 1)
    
  }
  
} else if (distr == 'farm') {
  
  smpl.n <- 50
  assoc = 'unknown'
  
  farm <- read.csv('farmData.csv')
  
  popu <- select(farm, c(AREA, TCC)) %>% filter(., AREA <= 6000)
  names(popu) <- c('size', 'y')
  popu <- popu[order(popu$size), ]
  popu <- mutate(popu, probs = (size/sum(size))*smpl.n)
  
} else {
  
  smpl.n <- 50
  distr = 'ndatss'
  assoc = 'unknown'
  
  popu <- read.csv('popu_c.csv')
}


summary(popu)

## ------- Simulaton Setting: Transformation and Model Options ---

transfms <- c('sqrt', 'sqrt') # c('none', 'logit', 'sqrt', 'log')
mod <- 'SN-B2PSP-uniform' # c('SN-B2PSP', 'SN-B2PSP-cauchy', 'SN-B2PSP-uniform', 'SN-B2PSP-omega', 'SN-poly', 'SN-poly-cauchy', 'SN-poly-uniform')  

# mod <- 'SN-B2PSP-uniform'

# ------- Miscellaneous ----------------------------------

if (multitask) {
  taskid <- Sys.getenv("SGE_TASK_ID") %>% as.integer(.); taskid
}

if (taskid != 1) {
  plot.popu <- FALSE
}

# ------- Simulate Random Seeds ---------------------------

set.seed(taskid)
seeds <- sample( ( (taskid - 1)*1e4 + 1):(taskid*1e4), nrep, replace = FALSE); print(seeds)

infr <- sapply(1:nrep, function (b) {
  print(b)
  SNsurvyBayes (dat = popu, transfm = transfms, mod, smpl.n, num.knots = 15, n.iter = 6e3, n.chains = 3, seed = c(seeds[b], 1), percentiles = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), doSample = TRUE, doFigure = (b == 1)*plot.popu )
}, simplify = 'array')

size.ty 
distr
assoc 
trans 

transfms  # c('none', 'logit', 'sqrt', 'log')
mod  # c('SN-B2PSP', 'SN-poly')

dim(infr)

save(infr, file = paste(distr, size.ty, assoc, transfms[1], transfms[2], mod, taskid, sep = "_") )

end.time <- Sys.time()

end.time - start.time

print(infr)