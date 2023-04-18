# Prep nimble data, constants, and initial values for Chile Puma data for 
#   gSPIM model following Augustine et al. 2021

library(coda)
library(nimble)
library(abind)

source("D:/Puma/genotype-SPIM-main/build.genos.R")
source("D:/Puma/genotype-SPIM-main/init.data.poisson.sampType.R")
source("D:/Puma/genotype-SPIM-main/map.genos.R")
source("D:/Puma/genotype-SPIM-main/NimbleModel genoSPIM Poisson sampType.R")
source("D:/Puma/genotype-SPIM-main/Nimble Functions genoSPIM Poisson sampType.R")
source("D:/Puma/genotype-SPIM-main/sSampler.R")



################################################################################
# Identify data 
M <- 500
J <- nrow(data$X)
K <- data$K
n.cov <- data$n.cov
n.cov.use <- n.cov
n.levels <- data$n.levels
n.rep <- data$n.rep
n.samples <- data$n.samples
K1D <- data$K1D
samp.levels <- 2 #number of sample type covariates. Each type has it's own genotyping error rates.
samp.type <- data$samp.type

built.genos <- build.genos(unique_genos)
ptype <- built.genos$ptype 
ptypeArray <-  built.genos$ptypeArray

gammaMat <- matrix(0, nrow = n.cov, ncol = max(n.levels))
for(l in 1:n.cov){
	gammaMat[l,1:n.levels[l]] = rep(1/n.levels[l], n.levels[l])
	}


################################################################################
# Initial values for latent variables and other data structures; genotyping error rates 
#  stored in a ragged matrix
p.geno.het.init=matrix(NA, nrow = 2, ncol = 3)
p.geno.hom.init=matrix(NA, nrow = 2, ncol = 2)
p.geno.het.init[1,] = c(0.95, 0.025, 0.025)
p.geno.het.init[2,] = c(0.95, 0.025, 0.025)
p.geno.hom.init[1,] = c(0.95, 0.05)
p.geno.hom.init[2,] = c(0.95, 0.05)

inits <- list(lam0 = 1, 
	sigma = 1,
	gammaMat = gammaMat,
	p.geno.het = p.geno.het.init,
	p.geno.hom =  p.geno.hom.init)
	
nimbuild <- init.data.poisson.sampType(data = data,
	M = M,
	inits = inits)



################################################################################
# Padding the structures for n.rep == 1
if(n.rep == 1){
	n.rep.use = 2
		} else {
			n.rep.use = n.rep
			}


################################################################################
# Bundle initial values
Niminits <- list(z = rep(1,M),
	s = nimbuild$s,
	G.true = nimbuild$G.true,
	ID = nimbuild$ID,
	capcounts = rowSums(nimbuild$y.true),
	y.true = nimbuild$y.true,
	G.latent = nimbuild$G.latent,
	theta = nimbuild$thetaArray,
	#psi = 0,
	lam0 = inits$lam0, # 0.02,   #inits$lam0,
	sigma = inits$sigma) #5)   #inits$sigma)
	
	
# Bundle constants
constants <- list(M = M,
	J = J,
	K = K,
	K1D = K1D,
	n.rep = n.rep.use,
	n.samples = n.samples,
	n.cov = n.cov,
	xlim = nimbuild$xlim,
	ylim = nimbuild$ylim,
	na.ind = nimbuild$G.obs.NA.indicator,
	n.levels = n.levels,
	max.levels = max(n.levels),
	ptype = ptypeArray,
	area_sq_km = data$area_sq_km) 


# Bundle  data
Nimdata <- list(y.true = matrix(NA, nrow = M, ncol = J),
	G.obs = nimbuild$G.obs,
	G.true = matrix(NA, nrow = M, ncol = n.cov),
	ID = rep(NA, n.samples),
	samp.type = samp.type,
	z = rep(NA, M),
	X = as.matrix(data$X),
	capcounts = rep(NA, M))


##############################################################################
##############################################################################