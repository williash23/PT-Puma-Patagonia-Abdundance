 #  Formatting Patanogia Puma data for SPIM model following Augustine et al. 2021

library(tidyverse)
library(sf)
library(ggplot2)
library(raster)
library(gdalUtils)
library(wesanderson)
library(terra)

setwd("D:/Puma/Torres_del_Paine/gSPIM")



################################################################################
# Using 10km state space buffer buffer per Proffitt et al. 2015
state_space_buff <- 10

# Grid cell size 
grid_cell_m <- 4000


################################################################################
#  Load partial genotype data
caps_tmp <- read.csv("Data/chile_puma_geno_w_sex_and_samptype.csv") %>%
	dplyr::filter(Lab.ID != "Pan19-S200") #### No location info for this sample

#  Needed variables:
#   sample ID (sample)
#   site ID (site)
#   occasion ID (survey) -- need to get this
#   genotype scores starting at column 18
head(caps_tmp)

#  Number of loci
n_loci <- 19
#n_loci <- 20 # includes a loci for sex


# Keep sample type: 1 (if ID'd by lab to individual), 2 (if poor DNA)
samptype_tmp <- caps_tmp$Samptype

#  Pull out the genotype scores - store in a list with 1 list element per locus
genos <-  caps_tmp[,3:40]
#genos <-  caps_tmp[,3:42] # includes columns for two sex alleles
genos_list <- vector("list",n_loci)
first_geno_cols <- seq(from = 1, to = 38, by = 2)
#first_geno_cols <- seq(from = 1, to = 40, by = 2)

for(i in 1:length(first_geno_cols)){
	genos_list[[i]] <- as.matrix(genos[,first_geno_cols[i]:(first_geno_cols[i]+1)])
	}

#  Items needed for model fit:
# 1: a list of all unique alleles at each locus (alleles)
# 2: a list of all unique genotypes stemming from all possible combinations of alleles at each 
#   locus (unique_genos )
# 3: a list indicating which genotype is a heterozygote or homozygote at each locus
#   (zygotype)
# 4: a list indicating which observations|true genotypes constitute correct classification, 
#   allelic dropout, or false allele events
#   1 is correct
#   0 is not possible (allelic dropout for homozygotes)
#   2 is allelic dropout
#   3 is false allele

alleles <- unique_genos  <- ptype <- zygotype <- vector("list", n_loci)

for(i in 1:n_loci){
	alleles[[i]] <- sort(unique(c(genos_list[[i]])))
	chk_0 <- alleles[[i]][1]
	if(chk_0 == 0){
		alleles[[i]] <- alleles[[i]][-1]
		}else{
			alleles[[i]] <- alleles[[i]]
			}
  
	unique_genos [[i]] <- expand.grid(alleles[[i]],alleles[[i]])
	unique_genos [[i]] <- t(apply(unique_genos [[i]],1,sort))
	dupes <- which(duplicated(apply(unique_genos [[i]],1,paste,collapse = "")))
	unique_genos [[i]] <- unique_genos [[i]][-dupes,]
	
	zygotype[[i]] <- ifelse(apply(unique_genos [[i]],1,function(x){x[1] == x[2]}),1,2)
	# 1 is homozygote, 2 hetero
	ptype[[i]] <- matrix(3,nrow = nrow(unique_genos [[i]]), ncol = nrow(unique_genos [[i]]))
	diag(ptype[[i]]) <- 1
	
	for(j in 1:nrow(unique_genos [[i]])){
		for(k in 1:nrow(unique_genos [[i]])){
			if(j==k) next
			if(length(unique(unique_genos [[i]][k,]))==2){ # heterozygote
			if(length(unique(unique_genos [[i]][j,]))==1){ # homozygote
			if(sum(unique_genos [[i]][j,]==unique_genos [[i]][k,])==1){ # could have been allelic dropout
				ptype[[i]][j,k] <- 2
				}
				}
				}
			}
		}
	}	

#  Examples of the formatted data
#  These are the observed alleles at locus 1
alleles[[1]]
#  These are all possible genotypes at locus 1, given these alleles.
#   This is the order in which the locus-level genotypes are enumerated. 1, ..., 3 in this case
unique_genos [[1]]
#  This indicates which unique genotypes at locus 1 is a homozygote (1) and heterozygote (2)
zygotype[[1]]
#  This indicates which type of observation it is if you observe the genotype X (rows) 
#   given the true genotype Y (columns):
#   1 = correct (on the diagonal)
#   2 = allelic dropout (only possible for zygotype[[1]]=1)
#   3 = false allele
ptype[[1]]



################################################################################
# Scat location data
scat_dat_tmp <- read.csv("Data/scat_locs.csv") %>%
	st_as_sf(coords = c("Long", "Lat"), crs = 4326)
scat_sf <- scat_dat_tmp %>%
	st_transform(32718) %>%
	dplyr::filter(Sample %in% caps_tmp$Lab.ID) %>%
	mutate(tmpID = 1:n()) %>%
	dplyr::select(Sample, Individual, tmpID)

# Confirm IDs the same in the scat location data and the scat genetic data
tst1 <- scat_sf$Sample
tst2 <- caps_tmp$Lab.ID
tst1 %in% tst2


#  Transects
tracks_sf <- st_read("Data/Loops.shp") %>%
	st_zm() %>%
	st_transform(32718)

# Distances between scats and tracks
scat_track_dist <- scat_sf %>%
	st_distance(tracks_sf) %>%
	as.data.frame()
min_dist <- apply(scat_track_dist, 1, FUN = min)
max_min_dist <- max(min_dist)
	
tracks_sf_b <- tracks_sf %>%
	st_buffer(100) %>%
	st_combine(.)
	
	
	
################################################################################	
# Make grid around transects
cell_size <- grid_cell_m

pr <- crs(as(tracks_sf, "Spatial"))

xmin <- as.numeric(st_bbox(tracks_sf_b)[1])
xmax <- as.numeric(st_bbox(tracks_sf_b)[3]) + cell_size  # buffer the xmax by cell size to ensure the grid covers all tracks
ymin <- as.numeric(st_bbox(tracks_sf_b)[2])
ymax <- as.numeric(st_bbox(tracks_sf_b)[4]) + cell_size # buffer the ymax by cell size to ensure the grid covers all tracks
x <- seq(from = xmin, to = xmax, by = cell_size)
y <- seq(from = ymin, to = ymax, by = cell_size)
xy <- expand.grid(x = x, y = y)
grid.pts <- SpatialPointsDataFrame(coords= xy, data=xy, proj4string = pr)
gridded(grid.pts) <- TRUE
grid <- as(grid.pts, "SpatialPolygons")
grid_sf_tmp <- grid %>%
	st_as_sf() %>%
	mutate(site = 1:n())

#  Intersect grid and transects to know which grid cells were sampled
grid_track_int <- st_intersects(grid_sf_tmp, tracks_sf, sparse = FALSE)
grid_track_int2 <- st_intersection(grid_sf_tmp, tracks_sf) %>%
	mutate(len = st_length(.)) %>%
	mutate(len_m = as.numeric(len)) %>%
	as.data.frame() %>%
	dplyr::select(site, len_m) 

kp <- grid_track_int2 %>%
	group_by(site) %>%
	mutate(site_len_m = sum(len_m)) %>%
	slice(1) %>%
	as.data.frame() %>%
	dplyr::filter(site_len_m >= 500)

grid_sf <- grid_sf_tmp %>%
	#dplyr::filter(site %in% kp) %>%
	dplyr::filter(site %in% kp$site) %>%
	mutate(site = 1:n())

# Check plot
p_kp <- ggplot() +
	geom_sf(data = grid_sf_tmp, fill = "grey50", colour = "black", size = 1, alpha = 0.5) +
	geom_sf(data = grid_sf, aes(fill = factor(site)), colour = "black", size = 1.5, alpha = 0.5) +
	geom_sf(data = tracks_sf, colour = "darkred", size = 1.5) +
	theme_bw() 
p_kp



################################################################################
# XY centroids of traps (ie grid cells) for model	
traps <- grid_sf %>%
	st_as_sf() %>%
	st_centroid() %>%
	mutate(utmw.trap = st_coordinates(.)[,1],
		utms.trap = st_coordinates(.)[,2]) %>%
	as.data.frame() %>%
	dplyr::select(site, utmw.trap, utms.trap)

#  Scat locations intersected with grid
grid_scat_int <- st_intersection(scat_sf, grid_sf)

#  Issue when grid cells are too small - scats are not intersecting with the grid.....
tst3 <- grid_scat_int$Sample
which(tst2 %in% tst3 == FALSE)

#  Connect scat location and grid site to the capture genetic data
caps_sf <- grid_scat_int %>%
	left_join(caps_tmp, by = c("Sample" = "Lab.ID"))
caps <- caps_sf %>%
	as.data.frame() %>%
	dplyr::select(-Individual, -geometry)
	
# Confirm IDs the same in the interesected data and the scat genetic data
tst4 <- caps$Sample
tst1 %in% tst4


v <- vect(grid_sf_tmp)
r <- rast(v, nrows = 250, ncols = 250)

tracks_u <- tracks_sf %>%
	st_simplify(preserveTopology = FALSE, 150) %>%
	st_buffer(100) %>%
	st_union(.) 
tracks_v <- vect(tracks_u)

tracks_r <- rasterize(tracks_v, r)
plot(tracks_r)
pts <- crds(tracks_r) %>%
	as.data.frame() %>%
	st_as_sf(coords = c(1,2), crs = 32718) 
l <- pts %>%
	st_cast("MULTILINESTRING")



################################################################################
# Generate an effort covaraite based on length of transect walked within grid cell
#  Discretize length into "surveys"
eff <- st_intersection(grid_sf, pts) %>%
	group_by(site) %>%
	mutate(grid_eff = n()) %>%
	slice(1) %>%
	as.data.frame()
max(eff$grid_eff)
# 84 points in a cell --- divide by 15

eff_check <- eff 	%>%
	mutate(effort_exact = grid_eff/15) %>%
	mutate(effort = round(effort_exact)) %>%
	dplyr::select(site, effort)
eff_check$effort[eff_check$effort == 0] <- 1
eff_grid <- grid_sf %>% left_join(eff_check, by = "site")	
# max 6 effort
summ_grid_eff <- eff_grid %>% 
	group_by(site) %>% 
	mutate(max_eff = max(effort)) %>% 
	slice(1) %>% 
	as.data.frame()

eff_grid_w <- eff_grid %>%
	mutate(transect = 1) %>%
	complete(site, effort) %>% #, fill = list(transect = 0)) %>%
	arrange(site, effort) %>%
	as.data.frame() %>%
	dplyr::select(-geometry) %>%	
	pivot_wider(names_from = effort, values_from = transect) %>%
	as.data.frame() 
colnames(eff_grid_w) <- c("site", "1", "2", "3", "4", "5", "6")
eff_grid_w$"1" <- 1
eff_grid_w $"5"[eff_grid_w $"6" == 1] <- 1
eff_grid_w $"4"[eff_grid_w $"5" == 1] <- 1
eff_grid_w $"3"[eff_grid_w $"4" == 1] <- 1
eff_grid_w $"2"[eff_grid_w $"3" == 1] <- 1

# Randomly select which "transect" the scat was found at
caps_vec <- vector()

for(i in 1:nrow(grid_scat_int)){
	cap_tmp <- eff_grid %>%
		dplyr::filter(site == grid_scat_int$site[i])
	cap_seg <- sample(seq(from = 1, to = cap_tmp$effort), 1)
	caps_vec[i] <- cap_seg
	}



################################################################################
#  Summary: number of scats per individual
scat_n <- scat_sf %>%
	group_by(Individual) %>%
	mutate(n_recaps = n()) %>%
	as.data.frame() %>%
	st_as_sf()
scat_n$n_recaps[scat_n$Individual == "Poor DNA"] <- 1

summ_recaps <- scat_n %>% group_by(Individual) %>% slice(1) %>% as.data.frame()

#  Create a polygonal state space to reduce necessary 
#   level of data augmentation
ss_pl <- grid_sf %>%
	st_union(.) %>%
	st_buffer(10000) 
area_sq_km <- as.numeric(st_area(ss_pl))/1000000



################################################################################
#  Combine processed  data into an object for gSPIM model
X <- traps[,2:3]/1000 # working in KM
J <- nrow(X)
#occ <- rep(1, nrow(caps))
occ <- caps_vec
#K <- 1
K <- 6
#tf <- rep(1, nrow(traps))
capsites <- caps$site
#K1D <- rep(K,J) #trap operation matrix, number of occasions trap j is operable
K1D <- eff_grid$effort #trap operation matrix, number of occasions trap j is operable

samples <- caps$Sample
reps <- 1
allsamps <- unique(samples)
n_samples <- length(allsamps)


#  Construct G_obs
G_obs <- array(NA, dim = c(n_samples, n_loci, reps))
trap2 <- occ2 <- rep(NA, n_samples)

#  Had to adjust this part from Augustine example
for(i in 1:n_samples){
	
	idx = which(samples == allsamps[i])
	occ2[i] = occ[idx][1]
	trap2[i] = caps$site[idx][1]
	
	for(j in 1:n_loci){
		
			#scores = genos_list[[j]][idx,]
			scores = cbind(genos_list[[j]][idx,1], genos_list[[j]][idx,2])

			for(k in 1:nrow(scores)){
				
				#if(sum(genos_list[[j]][idx,][k,]) == 0) next	 
					 #G_obs[i,j,k] = which(genos_list[[j]][idx,][k,1] == unique_genos [[j]][,1] &
                           #genos_list[[j]][idx,][k,2] == unique_genos [[j]][,2])
				if(sum(scores[k,]) == 0) next	
					G_obs[i,j,k] = which(scores[k,1] == unique_genos [[j]][,1] &
						scores[k, 2] == unique_genos [[j]][,2])
					 
				
				}
			}
		}


#  Check G_obs
dim(G_obs) #sample by locus by replication number
G_obs[1,,] 

y_obs <- array(0, dim = c(n_samples, J, K))
for(i in 1:n_samples){
	y_obs[i, trap2[i], occ2[i]] = 1
	}

#  Look at y_obs
#  dim = n_samples, ngrid cells that had effort, nmax transects
dim(y_obs) 
all(rowSums(y_obs) == 1)

# Create sample level covariate based on whether lab assigned ID
samptype <- samptype_tmp

# Get this.j, this.k
#  J = sites
#  K = surveys
tmp <- t(apply(y_obs, 1, function(x){which(x == 1, arr.ind = TRUE)}))
this.j <- tmp[,1]
this.k <- tmp[,2]

# IDlist
IDcovs = vector("list", n_loci)
for(i in 1:n_loci){
	IDcovs[[i]] = 1:nrow(unique_genos [[i]])
	}

#  Number loci-level genotypes per locus
n_levels <- unlist(lapply(unique_genos, nrow)) 
	
	
	
################################################################################
#  Bundle data
IDlist <- list()
IDlist$ncat <- n_loci
IDlist$IDcovs <- IDcovs
IDlist$zygotype <- zygotype
IDlist$ptype <- ptype

data <- list()
data$y.obs <- y_obs
data$G.obs <- G_obs
data$X <- X
data$IDlist <- IDlist
data$samp.type <- as.integer(samptype)

data$K <- K
data$n.rep <- reps
data$K1D <- K1D 
data$this.k <- this.k
data$this.j <- this.j
data$n.samples <- n_samples
data$n.levels <- n_levels
data$n.cov <- n_loci

data$buff <- state_space_buff #(10km)
data$area_sq_km <- area_sq_km 



################################################################################
# Visualizations
pal <- c("grey50", wes_palette(20, name = "Zissou1", type = "continuous"))

p0 <- ggplot() +
	geom_sf(data = grid_sf, colour = "black", fill = "transparent") +
	geom_sf(data = tracks_sf, colour = "dodgerblue4") +
	geom_sf(data = scat_n, alpha = 0.75, size = 4, aes(colour = factor(Individual), shape = factor(n_recaps))) +
	geom_sf(data = ss_pl, colour = "red", fill = "transparent") +
	scale_colour_manual(name = "Individual", values = pal) +
	scale_shape_manual(name = "Num. captures per individual", values = c(16,17)) +
	theme_bw() + 	
	theme(text = element_text(size = 14)) 
	#xlim(c(-73.1, -72.5)) +
	#ylim(c(-51.25, -50.84)) 
p0
#ggsave(file = "grid_statespace_and_captures.jpeg")


p1 <- ggplot() +
	#geom_sf(data = eff_grid, aes(colour = effort, fill = effort)) +
	geom_sf(data = eff_grid, colour = "black", aes(fill = effort)) +
	geom_sf(data = pts, colour = "black", size = 0.75, alpha = 0.5) +
	#scale_colour_gradient(name = "Effort", low = "yellow", high = "red", na.value = NA) +
	scale_fill_gradient("Effort", low = "yellow", high = "red", na.value = NA) +
	theme_bw() +
	theme(text = element_text(size = 14)) 
p1
#ggsave(file = "grid_effort.jpeg")


##############################################################################
##############################################################################
