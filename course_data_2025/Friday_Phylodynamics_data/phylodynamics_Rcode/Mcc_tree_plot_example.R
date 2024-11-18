# Plain R script of Mcc_tree_plot_example
# S J Lycett
# 20 Feb 2024

#######################################################
# STEP 0
#######################################################

# uncomment to install if required
#install.packages("ape")
#install.packages("maps")
#install.packages("mapdata")
#install.packages("mapproj")
#install.packages("conicfit")

# now load required packages

# base trees
library(ape)

# maps
library(maps)
library(mapdata)
library(mapproj)

# needed to fit HPDs to ellipses
library(conicfit)


# load Lycett_phylo utility code as separate functions
source("getEl.R")
source("calcDecimalDate.R")
source("read_MCC_tree.R")
source("get_BEAST_cols.R")
source("custom_map_movie.R")

#######################################################
# STEP 1
#######################################################
# define tree file name
treeName <- "cov_net_sim_mper2_120genomes_TN93G4_strict_skygrid_traits_2_mcc.tre"

# read in tree using custom function (could use ggtree also)
tr <- read_latlon_mcc_tr(treeName)
tr <- addDiscreteTraits(tr)
tr <- addLatLonHPD(tr)
tr <- fit_HPDs_to_standard(tr,ltol=0.005)
tr$ntips <- length(tr$tip.label)

# the tree object
print(tr)

# extra attributes
print(attributes(tr))

# discrete traits added
print(tr$propNames)

# discrete trait 1
propIndex <- 1
table(tr$props[1:tr$ntips,propIndex])

# discrete trait 2
propIndex <- 2
table(tr$props[1:tr$ntips,propIndex])

# min-max of the lat-lons
lltbl <- c(min(tr$latlon[,1]),max(tr$latlon[,1]),min(tr$latlon[,2]),max(tr$latlon[,2]))
names(lltbl) <- c("Min Lat","Max Lat","Min Lon","Max Lon")
lltbl

# raw plot of lat-lons
plot(tr$latlon[,2], tr$latlon[,1], xlab="Longitude", ylab="Latitude")
title("Raw plot of spatial coordinates")

#######################################################
# STEP 2
#######################################################

# plot plain tree
tr <- ladderize(tr)
plot(tr, cex=0.5)
add.scale.bar()

# plot tree with discrete trait 1
propIndex <- 1
plot_discrete_tree(tr, propIndex=propIndex, show.tip.label=FALSE)

# plot tree with discrete trait 2
propIndex <- 2
plot_discrete_tree(tr, propIndex=propIndex, show.tip.label=FALSE)

#######################################################
# STEP 3
#######################################################

# plot tree on map with discrete trait 1
plot_mcc_tree_with_hpds(tr, propIndex=1, xlim=c(-15,15), ylim=c(47,65), legpos="topleft")


# plot tree on map with discrete trait 2
plot_mcc_tree_with_hpds(tr, propIndex=2, xlim=c(-15,15), ylim=c(47,65), legpos="topleft")

#######################################################


