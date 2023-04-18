#  Visualize outputs from NIMBLE gSPIM model; assess differences in varying number of 
#   surveys and grid cell size


library(coda)
library(nimble)
library(abind)
library(ggplot2)
library(bayesplot)
library(MCMCglmm)
library(patchwork)

options(scipen = 999)
options(digits = 6)

setwd("D:/Puma/Torres_del_Paine/gSPIM")



################################################################################
# Main scenario: 4x4km grid cell size and max 6 surveys

load("mvSamples_TDP_chain1_4km_grid_6occ.Rdata")
res1 <- mvSamples
rm(mvSamples)

load("mvSamples_TDP_chain2_4km_grid_6occ.Rdata")
res2 <- mvSamples
rm(mvSamples)

load("mvSamples_TDP_chain3_4km_grid_6occ.Rdata")
res3 <- mvSamples
rm(mvSamples)

# Manual burn in and thin
out <- mcmc.list(mcmc(res1[1:nrow(res1),]), 
	mcmc(res2[1:nrow(res1),]), 
	mcmc(res3[1:nrow(res1),]))

#  Traceplots
mcmc_hist(out)
#ggsave("traceplots.png")

#  Summary
MCMCvis::MCMCsummary(out)

#  Coefficient plots
bayesplot_theme_update(text = element_text(size = 16, family = "sans"))

bayesplot::color_scheme_set("brightblue")	
pn <- bayesplot::mcmc_intervals(out, 
	prob_outer = 0.95,
	point_est = "median",
	prob = 0.5,
	pars = c("n")) +
	xlab("Parameter estimate") +
	#ylab("") +
	scale_y_discrete(labels = c("n" = "Unique individs.\nsampled")) 
pn

pN <- bayesplot::mcmc_intervals(out, 
	prob_outer = 0.95,
	point_est = "median",
	prob = 0.5,
	pars = c("Dens_100sq_km")) +
	xlab("Parameter estimate") +
	#ylab("") +
	scale_y_discrete(labels = c("Dens_100sq_km" = "Density\n(ind. per 100 sq. km.)")) 
pN

bayesplot::color_scheme_set("teal")	
pS <- bayesplot::mcmc_intervals(out, 
	prob_outer = 0.95,
	point_est = "median",
	prob = 0.5,
	pars = c("sigma")) +
	xlab("Parameter estimate") +
	ylab("") +
	scale_y_discrete(labels = 
	c("sigma" = "Home range\nspatial scale")) 
pS

pD <- bayesplot::mcmc_intervals(out, 
	prob_outer = 0.95,
	point_est = "median",
	prob = 0.5,
	pars = c("psi", "lam0")) +
	xlab("Parameter estimate") +
	ylab("") +
	scale_y_discrete(labels = 
	c("psi" = "Individual inclusion\nprobability",
		"lam0" = "Baseline\ndetection rate")) 
pD

pN/pS/pD

#ggsave("parameter_estimates.png")



#################################################################################################
# Explore ID posteriors
load("mvSamples2_OP_chain1_4km_grid_6occ.Rdata")
res1 <- as.matrix(mvSamples2)

load("mvSamples2_OP_chain2_4km_grid_6occ.Rdata")
res2 <- as.matrix(mvSamples2)

load("mvSamples2_OP_chain3_4km_grid_6occ.Rdata")
res3 <- as.matrix(mvSamples2)


plot(mcmc(res1[2:nrow(res1),]))

out2 <- res1 %>%
	rbind(res2) %>%	
	rbind(res3)

IDpost <- posterior.mode(out2[100:nrow(out2),])

#calculate posterior probability of pairwise sample matches
#P(sample x belongs to same individual as sample y)
burnin = 2 #where to start. Don't start at 1, is NA.
n.iter = nrow(out2) - burnin+1
pair.probs = matrix(NA, data$n.samples, data$n.samples)

for(i in 1:data$n.samples){
	for(j in 1:data$n.samples){
	count = 0
	for(iter in burnin:n.iter){
		count=count+1*(out2[iter,j] == out2[iter,i])
		}
		pair.probs[i,j] = count/(n.iter-1)
		}
	}

this.samp = 1 #sample number to look at
pair.probs[this.samp,] #probability this sample is from same individual as all other samples
pair.probs[this.samp,data$ID==data$ID[this.samp]] #for simulated data, these are the other samples truly from same individual

i = 1

i = i + 1
i
max(pair.probs[i,][pair.probs[i,] < 1])
which(pair.probs[i,] == max(pair.probs[i,][pair.probs[i,] < 1]))







#####################################################################################################
#  Other grid size and effort scenarios
#   Vary grid cell size: 8x8, 4x4, 2x2 km

# 8x8 km
load("mvSamples_TDP_chain1_8km_grid_6occ.Rdata")
res1 <- mvSamples
rm(mvSamples)

load("mvSamples_TDP_chain2_8km_grid_6occ.Rdata")
res2 <- mvSamples
rm(mvSamples)

load("mvSamples_TDP_chain3_8km_grid_6occ.Rdata")
res3 <- mvSamples

# Manual burn in and thin
out8 <- mcmc.list(mcmc(res1[1:nrow(res1),]), 
	mcmc(res2[1:nrow(res1),]), 
	mcmc(res3[1:nrow(res1),]))
rm(res1, res2, res3)


# 4x4 km
load("mvSamples_TDP_chain1_4km_grid_6occ.Rdata")
res1 <- mvSamples
rm(mvSamples)

load("mvSamples_TDP_chain2_4km_grid_6occ.Rdata")
res2 <- mvSamples
rm(mvSamples)

load("mvSamples_TDP_chain3_4km_grid_6occ.Rdata")
res3 <- mvSamples

# Manual burn in and thin
out4 <- mcmc.list(mcmc(res1[1:nrow(res1),]), 
	mcmc(res2[1:nrow(res1),]), 
	mcmc(res3[1:nrow(res1),]))
rm(res1, res2, res3)


# 2x2 km
load("mvSamples_TDP_chain1_2km_grid_6occ.Rdata")
res1 <- mvSamples
rm(mvSamples)

load("mvSamples_TDP_chain2_2km_grid_6occ.Rdata")
res2 <- mvSamples
rm(mvSamples)

load("mvSamples_TDP_chain3_2km_grid_6occ.Rdata")
res3 <- mvSamples

# Manual burn in and thin
out2 <- mcmc.list(mcmc(res1[1:nrow(res1),]), 
	mcmc(res2[1:nrow(res1),]), 
	mcmc(res3[1:nrow(res1),]))
rm(res1, res2, res3)


#  Coefficient plots
bayesplot_theme_update(text = element_text(size = 14, family = "sans"))
bayesplot::color_scheme_set("teal")	

pN2 <- bayesplot::mcmc_intervals(out2, 
	prob_outer = 0.99,
	point_est = "median",
	prob = 0.85,
	pars = c("Dens_100sq_km")) +
	xlab(NULL) +
	xlim(0, 22) +
	scale_y_discrete(labels = c("Dens_100sq_km" = "2 km")) 
pN2

pN4 <- bayesplot::mcmc_intervals(out4, 
	prob_outer = 0.99,
	point_est = "median",
	prob = 0.85,
	pars = c("Dens_100sq_km")) +
	xlab(NULL) +
	xlim(0, 22) +
	scale_y_discrete(labels = c("Dens_100sq_km" = "4 km")) 
pN4

pN8 <- bayesplot::mcmc_intervals(out8, 
	prob_outer = 0.99,
	point_est = "median",
	prob = 0.85,
	pars = c("Dens_100sq_km")) +
	xlab("\nDensity (ind. per 100 sq. km.)") +
	xlim(0, 22) +
	scale_y_discrete(labels = c("Dens_100sq_km" = "8 km")) 
pN8


pN2/pN4/pN8 + 
	plot_annotation(title = 'Varying grid cell length; all use 6 (max) occasions of effort')
 
 

#####################################################################################################
#  Other effort scenarios
#   Vary effort discretization: max surveys 3, 6, 10

# Max 3 surveys
load("mvSamples_TDP_chain1_4km_grid_3occ.Rdata")
res1 <- mvSamples
rm(mvSamples)

load("mvSamples_TDP_chain2_4km_grid_3occ.Rdata")
res2 <- mvSamples
rm(mvSamples)

load("mvSamples_TDP_chain3_4km_grid_3occ.Rdata")
res3 <- mvSamples
rm(mvSamples)

# Manual burn in and thin
out3 <- mcmc.list(mcmc(res1[1:nrow(res1),]), 
	mcmc(res2[1:nrow(res1),]), 
	mcmc(res3[1:nrow(res1),]))
rm(res1, res2, res3)


# Max 6 surveys
load("mvSamples_TDP_chain1_4km_grid_6occ.Rdata")
res1 <- mvSamples
rm(mvSamples)

load("mvSamples_TDP_chain2_4km_grid_6occ.Rdata")
res2 <- mvSamples
rm(mvSamples)

load("mvSamples_TDP_chain3_4km_grid_6occ.Rdata")
res3 <- mvSamples
rm(mvSamples)

# Manual burn in and thin
out6 <- mcmc.list(mcmc(res1[1:nrow(res1),]), 
	mcmc(res2[1:nrow(res1),]), 
	mcmc(res3[1:nrow(res1),]))
rm(res1, res2, res3)


# Max 10 surveys
load("mvSamples_TDP_chain1_4km_grid_10occ.Rdata")
res1 <- mvSamples
rm(mvSamples)

load("mvSamples_TDP_chain2_4km_grid_10occ.Rdata")
res2 <- mvSamples
rm(mvSamples)

load("mvSamples_TDP_chain3_4km_grid_10occ.Rdata")
res3 <- mvSamples
rm(mvSamples)

# Manual burn in and thin
out10 <- mcmc.list(mcmc(res1[1:nrow(res1),]), 
	mcmc(res2[1:nrow(res1),]), 
	mcmc(res3[1:nrow(res1),]))
rm(res1, res2, res3)


#  Coefficient plots
bayesplot_theme_update(text = element_text(size = 14, family = "sans"))
bayesplot::color_scheme_set("teal")	

pN3 <- bayesplot::mcmc_intervals(out3, 
	prob_outer = 0.99,
	point_est = "median",
	prob = 0.85,
	pars = c("Dens_100sq_km")) +
	xlab(NULL) +
	xlim(0, 22) +
	scale_y_discrete(labels = c("Dens_100sq_km" = "3 occasions")) 
pN3

pN6 <- bayesplot::mcmc_intervals(out6, 
	prob_outer = 0.99,
	point_est = "median",
	prob = 0.85,
	pars = c("Dens_100sq_km")) +
	xlab(NULL) +
	xlim(0, 22) +
	scale_y_discrete(labels = c("Dens_100sq_km" = "6 occasions")) 
pN6

pN10 <- bayesplot::mcmc_intervals(out10, 
	prob_outer = 0.99,
	point_est = "median",
	prob = 0.85,
	pars = c("Dens_100sq_km")) +
	xlab("\nDensity (ind. per 100 sq. km.)") +
	xlim(0, 22) +
	scale_y_discrete(labels = c("Dens_100sq_km" = "10 occasions")) 
pN10


pN3/pN6/pN10 + 
	plot_annotation(title = 'Varying (max) occasions of effort; all use 4 km grid cell length')
 

##############################################################################
##############################################################################
