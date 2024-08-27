# Plain R script of Mcc_tree_incursions
# S J Lycett
# 20 Feb 2024

# the start of this is the same as Mcc_tree_plot_example.R

#######################################################
# STEP 0
#######################################################

# uncomment to install if required
#install.packages("ape")
#install.packages("maps")
#install.packages("mapdata")
#install.packages("mapproj")
#install.packages("conicfit")
#install.packages("lattice") - needed for levelplot

# now load required packages

# base trees
library(ape)

# maps - not needed here
#library(maps)
#library(mapdata)
#library(mapproj)

# needed to fit HPDs to ellipses - not needed here
#library(conicfit)

# needed for levelplot
library(lattice)

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
# need to add the decimal dates [not always necessary but need to do this in this example]
tr$decDates	   <- apply(as.matrix(apply(as.matrix(tr$tip.label), 1, getEl, ind=1, sep="\\|", fromEnd=TRUE)),1,calcDecimalDate_from_yymmdd,sep="-")
tr$youngestTip <- max(tr$decDates)
tr <- nodeTimes(tr, youngestTip=tr$youngestTip)

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

#######################################################
# STEP 2
#######################################################

# use the properties (props) to find location from-to

# first get the node indices from the edges of the tree (1=ancestral, 2=child)
fromNode    <- tr$edge[,1]
toNode      <- tr$edge[,2]

# also get the respective times of these edges (i.e. the time of the ancestral node & the time of the child node)
fromTime    <- tr$nodeTimes[fromNode]
toTime      <- tr$nodeTimes[toNode]
midEdgeTime <- (fromTime+toTime)/2

#########################
# Country is propIndex = 1
propIndex   <- 1
fromCountry <- tr$props[fromNode,propIndex]
toCountry   <- tr$props[toNode,propIndex]
# number of transitions from Country to Country over the whole tree
country_tbl <- table(fromCountry,toCountry)
print(country_tbl)

# for display, set diagonal to -1
cols <- c("grey90",topo.colors(100))
diag(country_tbl) <- -1
levelplot(country_tbl,col.regions=cols, main="Number of Transitions from-to Country")

#########################
# Place is propIndex = 2
propIndex   <- 2
fromPlace   <- tr$props[fromNode,propIndex]
toPlace     <- tr$props[toNode,propIndex]
# number of transitions from Place to Place over the whole tree
place_tbl <- table(fromPlace,toPlace)
print(place_tbl)

# for display, set diagonal to 0
cols <- c("grey90",topo.colors(100))
diag(place_tbl) <- -1
levelplot(place_tbl,col.regions=cols, main="Number of Transitions from-to Place")


