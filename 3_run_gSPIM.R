#  Compile and run gSPIM model in NIMBLE; run as 3 different instances to run and save separate 
#   chains

nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)
nimbleOptions('MCMCjointlySamplePredictiveBranches') 



################################################################################
# MCMC settings
nit <- 300000
nb <- 150000



################################################################################
#  Parameters to monitor
parameters <- c('psi', 'lam0', 'sigma', 
	'N', 'n', 
	'p.geno.het', 'p.geno.hom', 'gammaMat', 
	"Dens_sq_km", "Dens_100sq_km")

parameters2 <- c('ID', "G.true")
nt <- 1
nt2 <- 50 



################################################################################
# Generate model and specify samplers
# Can ignore warnings about 1) ID in constants 2) possible size mismatch for G.obs.
start.time <- Sys.time()
Rmodel <- nimbleModel(code = NimModel, 
	constants = constants, 
	data = Nimdata,
	check = FALSE,
	inits = Niminits)
	
conf <- configureMCMC(Rmodel,
	monitors = parameters, 
	thin = nt, 
	monitors2 = parameters2,
	thin2 = nt2,
	useConjugacy = TRUE) 

conf$removeSampler("G.obs")
conf$removeSampler("y.true")
conf$addSampler(target = paste0("y.true[1:",M,",1:",J,"]"),
	type = 'IDSampler', 
	control = list(M = M, 
		J = J, 
		K1D = K1D, 
		n.cov = n.cov, 
		n.levels = n.levels,
		n.samples = n.samples, 
		n.rep = n.rep, 
		samp.type = samp.type,
		this.j = nimbuild$this.j, 
		G.obs = data$G.obs,
		na.ind = nimbuild$G.obs.NA.indicator),
	silent = TRUE)

conf$removeSampler("G.true")
for(i in 1:M){
	for(m in 1:n.cov){ #don't need to update first cat bc it is mark status
		conf$addSampler(target = paste("G.true[",i,",",m,"]", sep=""),
		type = 'GSampler',
		control = list(i = i, 
			m = m, 
			n.levels = n.levels, 
			n.rep = n.rep,
			samp.type = samp.type,
			na.ind = nimbuild$G.obs.NA.indicator[,m,]), 
		silent = TRUE)
		}
	}

conf$removeSampler(paste("s[1:",M,", 1:2]", sep = ""))
for(i in 1:M){
	conf$addSampler(target = paste("s[",i,", 1:2]", sep = ""),
	type = 'sSampler',
	control = list(i = i, 
		xlim = nimbuild$xlim, 
		ylim = nimbuild$ylim,
		scale = 0.25, adaptInterval = 500),
	silent = TRUE)
	}

conf$removeSampler(c("lam0", "sigma"))
conf$addSampler(target = c(paste("lam0"), paste("sigma")),
	type = 'AF_slice',
	control = list(adaptive = TRUE),
	silent = TRUE)



################################################################################
# Compile 
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel) #, showCompilerOutput = TRUE)



################################################################################
# Run 
#   Can ignore nimble warnings about NA or NaN in ptype and theta
#   Can ignore nimble warnings about G.obs value NA or NaN, 
#   due to padding to keep dimensions constant for nimble
start.time2 <- Sys.time()
Cmcmc$run(nit, nburnin = nb) 
end.time <- Sys.time()
#end.time - start.time  
end.time - start.time2 


################################################################################
#  Extract samples ands save
mvSamples = as.matrix(Cmcmc$mvSamples)
mvSamples2 = as.matrix(Cmcmc$mvSamples2)

# Params of interest
idx <- grep("gammaMat", colnames(mvSamples))
mvSamples <- mvSamples[,-idx]

#save(mvSamples, file = "mvSamples_TDP_chain1_4km_grid_6occ.Rdata")
#save(mvSamples2, file = "mvSamples2_TDP_chain1_4km_grid_6occ.Rdata")
#save(mvSamples, file = "mvSamples_TDP_chain2_4km_grid_6occ.Rdata")
#save(mvSamples2, file = "mvSamples2_TDP_chain2_4km_grid_6occ.Rdata")
#save(mvSamples, file = "mvSamples_TDP_chain3_4km_grid_6occ.Rdata")
#save(mvSamples2, file = "mvSamples2_TDP_chain3_4km_grid_6occ.Rdata")


##############################################################################
##############################################################################