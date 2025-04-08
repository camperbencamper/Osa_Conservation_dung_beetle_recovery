library("ggplot2")
library("dplyr")
library("vegan")
library("tidyr")
library(nlme)
library(lme4)
library(car)
library(blmeco)
library(emmeans)
library("stringr")
library('DHARMa')
library("glmmTMB")
library("plyr")

set.seed(1)


#data from all years
#comm_data <- read.csv("TEMPORARY_frog_resto_plot_data_all_years.csv", header=TRUE)
data_dung <- read.delim("data_dung.txt")
comm_data <- read.delim("DB_all.txt")
colnames(comm_data)
#comm_data <- comm_data[,-1]

DB_species <- colnames(comm_data[-(1:5)])

#creating a new comm_data column just by plot treatment
#extracting treatments from plot_id column
comm_data$treatment <- substr(comm_data$site_id, 1, 1)

comm_data$treatment_year <- as.factor(paste(comm_data$treatment, comm_data$year, sep="_"))

comm_data$treatment <- factor(comm_data$treatment, levels= c("N", "L", "M", "H", "P"))
comm_data$year <- factor(comm_data$year, levels= c("2017", "2019", "2021"))
#comm_data$treatment_year <- factor(comm_data$treatment_year, levels= c("N_2017", "L_2017", "M_2017", "H_2017", "N_2019", "L_2019", "M_2019", "H_2019", "N_2021", "L_2021", "M_2021", "H_2021", "P_2017"))


#--------------------------------------------------------------------------------------------------------------------------

### converting community matrix into dissimilarity matrix for restoration analyses

#Bray-Curtis

comm_data_ordi <- comm_data
#dropping rows with species abundances == 0
comm_data_ordi <- comm_data_ordi[!!rowSums(comm_data_ordi[names(comm_data_ordi) %in% DB_species]),] #idk why this works but many other solutions didnt, drops rows with species abundances == 0
#!! is not a typo
row.names(comm_data_ordi) <- paste(comm_data_ordi$site_id, comm_data_ordi$year, sep="_")

p_dissim <- vegdist(comm_data_ordi[,names(comm_data_ordi) %in% DB_species], method="bray")
p_dissim <- data.frame(as.matrix(p_dissim))
#write.csv(p_dissim, "test_matrix.csv")
p_dissim <- p_dissim[,grepl("P", names(p_dissim))] #just retaining columns that are primary forest controls
p_dissim <- p_dissim[!grepl("P", row.names(p_dissim)),] #just removing rows that are primary forest controls
p_sim <- 1 - p_dissim #making dissimilarity values into similarity values
p_sim$similarity <- apply(p_sim,1,mean)
p_sim$similarity <- rowMeans(p_sim)
p_sim$ln_similarity <- log10(p_sim$similarity)
p_sim$plot_year <- row.names(p_sim)
p_sim <- p_sim %>% tidyr::separate(plot_year, c("site_id", "year")) #separating row names into columns, automatically chooses non-alphanumeric character to split by
p_sim$treatment <- stringr::str_extract(p_sim$site_id, "[A-Z]+" ) 
p_sim$treatment <- factor(p_sim$treatment, levels= c("N", "L", "M", "H"))
p_sim$year <- factor(p_sim$year, levels= c("2017", "2019", "2021"))

p_dissim$dissimilarity <- apply(p_dissim,1,mean)
p_dissim$dissimilarity <- rowMeans(p_dissim)
p_dissim$ln_dissimilarity <- log10(p_dissim$dissimilarity)
p_dissim$plot_year <- row.names(p_dissim)
p_dissim <- p_dissim %>% tidyr::separate(plot_year, c("site_id", "year")) #separating row names into columns, automatically chooses non-alphanumeric character to split by
p_dissim$treatment <- stringr::str_extract(p_dissim$site_id, "[A-Z]+" ) 
p_dissim$treatment <- factor(p_dissim$treatment, levels= c("N", "L", "M", "H"))
p_dissim$year <- as.integer(p_dissim$year)

#-------------------------------------------------------------------------------------------------------------------------

#Spatial Autocorrelation
# Install and load required package
#if (!require(leaflet)) install.packages("leaflet")
library(leaflet)

# Create a data frame with your coordinates
coords <- read.csv("DB_analysis_2-10-24/github/metadata/osa_resto_plots_site_coords.csv", header=TRUE)

map <- leaflet(coords) %>%
  addTiles() %>%  # Adds default OpenStreetMap tiles
  addMarkers(
    lng = ~LongD,
    lat = ~LatD,
    popup = ~site_id,   # Displays the site_id when you click on the marker
    label = ~site_id    # Also shows site_id as a label when hovering over the marker
  )

# Print the map
map


#---

# Install and load required package
#if (!require(geosphere)) install.packages("geosphere")
library(geosphere)

# Read the CSV file (adjust the file path if needed)
coords <- read.csv("DB_analysis_2-10-24/github/metadata/osa_resto_plots_site_coords.csv", header=TRUE) #54 objects

row.names(coords) <- coords$site_id
coords$treatment <- substr(coords$site_id, 1, 1)

# Inspect the first few rows to understand the column names
head(coords)

# Assume your CSV has columns named 'site', 'lat', and 'lon'
# If your column names differ, change them accordingly.
coord_matrix <- as.matrix(coords[c("LongD", "LatD")])

# Compute the pairwise distance matrix using the Haversine formula
# This returns distances in meters
dist_matrix <- distm(coord_matrix, fun = distHaversine)

# Assign site names as row and column names
rownames(dist_matrix) <- coords$site_id
colnames(dist_matrix) <- coords$site_id

#-------------

##########################################################################################################################

#----------------------------------------------------------------------------------------------

#MORAN'S I w/ RICHNESS

#SI III: Table S1

library(spdep)
# Assume your data frame is 'plot_data' with columns: richness, treatment, lat, lon.
coords2 <- data.frame(cbind(coords[coords$treatment!="P",]$LongD, coords[coords$treatment!="P",]$LatD))

rownames(coords2) <- rownames(coords)[1:40]

#using distance-based network 

threshold <- 150  # meters
nb_dist <- dnearneigh(coords2, 0, threshold)
lw_dist <- nb2listw(nb_dist, style="W", zero.policy=TRUE)
moran.mc(comm_data[comm_data$year=="2017" & comm_data$treatment!="P",]$richness, lw_dist, nsim = 999, alternative="greater", zero.policy=TRUE)
moran.mc(comm_data[comm_data$year=="2019",]$richness, lw_dist, nsim = 999, alternative="greater", zero.policy=TRUE)
moran.mc(comm_data[comm_data$year=="2021",]$richness, lw_dist, nsim = 999, alternative="greater", zero.policy=TRUE)


#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

#further parsing by year and treatment, Moran's I

coords_N <- cbind(coords[coords$treatment=="N",]$LongD, coords[coords$treatment=="N",]$LatD)
nb_N <- dnearneigh(coords_N, 0, threshold)
list_N <- nb2listw(nb_N, style="W", zero.policy=TRUE)
moran_N2017 <- moran.mc(comm_data[comm_data$year=="2017" & comm_data$treatment=="N",]$richness, list_N, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_N2017)
moran_N2019 <- moran.mc(comm_data[comm_data$year=="2019" & comm_data$treatment=="N",]$richness, list_N, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_N2019)
moran_N2021 <- moran.mc(comm_data[comm_data$year=="2021" & comm_data$treatment=="N",]$richness, list_N, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_N2021)

coords_L <- cbind(coords[coords$treatment=="L",]$LongD, coords[coords$treatment=="L",]$LatD)
nb_L <- dnearneigh(coords_L, 0, threshold)
list_L <- nb2listw(nb_L, style="W", zero.policy=TRUE)
moran_L2017 <- moran.mc(comm_data[comm_data$year=="2017" & comm_data$treatment=="L",]$richness, list_L, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_L2017)
moran_L2019 <- moran.mc(comm_data[comm_data$year=="2019" & comm_data$treatment=="L",]$richness, list_L, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_L2019) 
moran_L2021 <- moran.mc(comm_data[comm_data$year=="2021" & comm_data$treatment=="L",]$richness, list_L, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_L2021)

coords_M <- cbind(coords[coords$treatment=="M",]$LongD, coords[coords$treatment=="M",]$LatD)
nb_M <- dnearneigh(coords_M, 0, threshold)
list_M <- nb2listw(nb_M, style="W", zero.policy=TRUE)
moran_M2017 <- moran.mc(comm_data[comm_data$year=="2017" & comm_data$treatment=="M",]$richness, list_M, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_M2017)
moran_M2019 <- moran.mc(comm_data[comm_data$year=="2019" & comm_data$treatment=="M",]$richness, list_M, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_M2019) 
moran_M2021 <- moran.mc(comm_data[comm_data$year=="2021" & comm_data$treatment=="M",]$richness, list_M, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_M2021)

coords_H <- cbind(coords[coords$treatment=="H",]$LongD, coords[coords$treatment=="H",]$LatD)
nb_H <- dnearneigh(coords_H, 0, threshold)
list_H <- nb2listw(nb_H, style="W", zero.policy=TRUE)
moran_H2017 <- moran.mc(comm_data[comm_data$year=="2017" & comm_data$treatment=="H",]$richness, list_H, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_H2017)
moran_H2019 <- moran.mc(comm_data[comm_data$year=="2019" & comm_data$treatment=="H",]$richness, list_H, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_H2019) 
moran_H2021 <- moran.mc(comm_data[comm_data$year=="2021" & comm_data$treatment=="H",]$richness, list_H, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_H2021)


#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------

# SETTING UP MANTEL TESTS
### BRAy_CURTIS

BCdissim <- vegdist(comm_data_ordi[, names(comm_data_ordi) %in% DB_species], method="bray")
BCdissim <- data.frame(as.matrix(BCdissim))
BCsim <- 1 - BCdissim #making dissimilarity values into similarity values
#BCsim_2017$site_id <- row.names(BCsim_2017)
#BCsim_2017$treatment <- stringr::str_extract(BCsim_2017$site_id, "[A-Z]+" ) 
#BCsim_2017$treatment <- factor(BCsim_2017$treatment, levels= c("N", "L", "M", "H"))


N_2017 <- BCsim[
  grepl("N", rownames(BCsim)) & grepl("2017", rownames(BCsim)), 
  grepl("N", colnames(BCsim)) & grepl("2017", colnames(BCsim))]
#rownames(df) <- sub("_.*", "", rownames(df))
#colnames(df) <- sub("_.*", "", colnames(df))
N_2019 <- BCsim[
  grepl("N", rownames(BCsim)) & grepl("2019", rownames(BCsim)), 
  grepl("N", colnames(BCsim)) & grepl("2019", colnames(BCsim))]
N_2021 <- BCsim[
  grepl("N", rownames(BCsim)) & grepl("2021", rownames(BCsim)), 
  grepl("N", colnames(BCsim)) & grepl("2021", colnames(BCsim))]

L_2017 <- BCsim[
  grepl("L", rownames(BCsim)) & grepl("2017", rownames(BCsim)), 
  grepl("L", colnames(BCsim)) & grepl("2017", colnames(BCsim))]
L_2019 <- BCsim[
  grepl("L", rownames(BCsim)) & grepl("2019", rownames(BCsim)), 
  grepl("L", colnames(BCsim)) & grepl("2019", colnames(BCsim))]
L_2021 <- BCsim[
  grepl("L", rownames(BCsim)) & grepl("2021", rownames(BCsim)), 
  grepl("L", colnames(BCsim)) & grepl("2021", colnames(BCsim))]

M_2017 <- BCsim[
  grepl("M", rownames(BCsim)) & grepl("2017", rownames(BCsim)), 
  grepl("M", colnames(BCsim)) & grepl("2017", colnames(BCsim))]
M_2019 <- BCsim[
  grepl("M", rownames(BCsim)) & grepl("2019", rownames(BCsim)), 
  grepl("M", colnames(BCsim)) & grepl("2019", colnames(BCsim))]
M_2021 <- BCsim[
  grepl("M", rownames(BCsim)) & grepl("2021", rownames(BCsim)), 
  grepl("M", colnames(BCsim)) & grepl("2021", colnames(BCsim))]

H_2017 <- BCsim[
  grepl("H", rownames(BCsim)) & grepl("2017", rownames(BCsim)), 
  grepl("H", colnames(BCsim)) & grepl("2017", colnames(BCsim))]
H_2019 <- BCsim[
  grepl("H", rownames(BCsim)) & grepl("2019", rownames(BCsim)), 
  grepl("H", colnames(BCsim)) & grepl("2019", colnames(BCsim))]
H_2021 <- BCsim[
  grepl("H", rownames(BCsim)) & grepl("2021", rownames(BCsim)), 
  grepl("H", colnames(BCsim)) & grepl("2021", colnames(BCsim))]

P_2017 <- BCsim[
  grepl("P", rownames(BCsim)) & grepl("2017", rownames(BCsim)), 
  grepl("P", colnames(BCsim)) & grepl("2017", colnames(BCsim))]

#-------------------------------------------------------

#-------------------------------------------------------

#SI III: Table S3

# MANTEL TESTS

mantel(N_2017, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(N_2017)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(N_2017))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(N_2019, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(N_2019)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(N_2019))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(N_2021, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(N_2021)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(N_2021))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))

mantel(L_2017, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(L_2017)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(L_2017))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(L_2019, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(L_2019)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(L_2019))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(L_2021, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(L_2021)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(L_2021))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))

mantel(M_2017, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(M_2017)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(M_2017))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(M_2019, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(M_2019)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(M_2019))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(M_2021, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(M_2021)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(M_2021))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))

mantel(H_2017, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(H_2017)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(H_2017))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(H_2019, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(H_2019)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(H_2019))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(H_2021, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(H_2021)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(H_2021))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))

mantel(P_2017, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(P_2017)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(P_2017))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))


#######################################################################################
#######################################################################################

#By treatment and year lumped


y2017 <- BCsim[
  grepl("2017", rownames(BCsim)), 
  grepl("2017", colnames(BCsim))]
y2017 <- y2017[!grepl("P", rownames(y2017)),
               !grepl("P", colnames(y2017))]

y2019 <- BCsim[
  grepl("2019", rownames(BCsim)), 
  grepl("2019", colnames(BCsim))]

y2021 <- BCsim[
  grepl("2021", rownames(BCsim)), 
  grepl("2021", colnames(BCsim))]


mantel(y2017, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(y2017)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(y2017))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(y2019, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(y2019)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(y2019))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(y2021, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(y2021)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(y2021))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))

#nothing significant

#----


### MANTEL CORRELOGRAMS - BRAY-CURTIS


ddd<-c()
#ddd<-c(ddd,0)
ddd<-c(ddd,60)
ddd<-c(ddd,120)
ddd<-c(ddd,180)
ddd<-c(ddd,240)
ddd<-c(ddd,300)
ddd<-c(ddd,400)
ddd<-c(ddd,500)
ddd<-c(ddd,700)
ddd<-c(ddd,900)
ddd<-c(ddd,1100)
ddd<-c(ddd,1300)
ddd<-c(ddd,1500)
ddd<-c(ddd,1700)

#SI III: Figure S1

#Mantel Correlograms
BCcorrelog_2017<-mantel.correlog(y2017, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(y2017)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(y2017))))],
  nperm=9999,cutoff=FALSE, break.pts=ddd)
plot(BCcorrelog_2017, alpha=0.05)

BCcorrelog_2019<-mantel.correlog(y2019, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(y2019)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(y2019))))],
  nperm=9999,cutoff=FALSE, break.pts=ddd)
plot(BCcorrelog_2019, alpha=0.05)

BCcorrelog_2021<-mantel.correlog(y2021, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(y2021)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(y2021))))],
  nperm=9999,cutoff=FALSE, break.pts=ddd)
plot(BCcorrelog_2021, alpha=0.05)



#################################################################################################################
#################################################################################################################

### JACCARD MANTEL TESTS

#################################################################################################################
#################################################################################################################


### Mantel Tests by Year

Jdissim <- vegdist(comm_data_ordi[, names(comm_data_ordi) %in% DB_species], method="jaccard")
Jdissim <- data.frame(as.matrix(Jdissim))
Jsim <- 1 - Jdissim #making dissimilarity values into similarity values



N_2017 <- Jsim[
  grepl("N", rownames(Jsim)) & grepl("2017", rownames(Jsim)), 
  grepl("N", colnames(Jsim)) & grepl("2017", colnames(Jsim))]
#rownames(df) <- sub("_.*", "", rownames(df))
#colnames(df) <- sub("_.*", "", colnames(df))
N_2019 <- Jsim[
  grepl("N", rownames(Jsim)) & grepl("2019", rownames(Jsim)), 
  grepl("N", colnames(Jsim)) & grepl("2019", colnames(Jsim))]
N_2021 <- Jsim[
  grepl("N", rownames(Jsim)) & grepl("2021", rownames(Jsim)), 
  grepl("N", colnames(Jsim)) & grepl("2021", colnames(Jsim))]

L_2017 <- Jsim[
  grepl("L", rownames(Jsim)) & grepl("2017", rownames(Jsim)), 
  grepl("L", colnames(Jsim)) & grepl("2017", colnames(Jsim))]
L_2019 <- Jsim[
  grepl("L", rownames(Jsim)) & grepl("2019", rownames(Jsim)), 
  grepl("L", colnames(Jsim)) & grepl("2019", colnames(Jsim))]
L_2021 <- Jsim[
  grepl("L", rownames(Jsim)) & grepl("2021", rownames(Jsim)), 
  grepl("L", colnames(Jsim)) & grepl("2021", colnames(Jsim))]

M_2017 <- Jsim[
  grepl("M", rownames(Jsim)) & grepl("2017", rownames(Jsim)), 
  grepl("M", colnames(Jsim)) & grepl("2017", colnames(Jsim))]
M_2019 <- Jsim[
  grepl("M", rownames(Jsim)) & grepl("2019", rownames(Jsim)), 
  grepl("M", colnames(Jsim)) & grepl("2019", colnames(Jsim))]
M_2021 <- Jsim[
  grepl("M", rownames(Jsim)) & grepl("2021", rownames(Jsim)), 
  grepl("M", colnames(Jsim)) & grepl("2021", colnames(Jsim))]

H_2017 <- Jsim[
  grepl("H", rownames(Jsim)) & grepl("2017", rownames(Jsim)), 
  grepl("H", colnames(Jsim)) & grepl("2017", colnames(Jsim))]
H_2019 <- Jsim[
  grepl("H", rownames(Jsim)) & grepl("2019", rownames(Jsim)), 
  grepl("H", colnames(Jsim)) & grepl("2019", colnames(Jsim))]
H_2021 <- Jsim[
  grepl("H", rownames(Jsim)) & grepl("2021", rownames(Jsim)), 
  grepl("H", colnames(Jsim)) & grepl("2021", colnames(Jsim))]

P_2017 <- Jsim[
  grepl("P", rownames(Jsim)) & grepl("2017", rownames(Jsim)), 
  grepl("P", colnames(Jsim)) & grepl("2017", colnames(Jsim))]

#-------------------------------------------------------

#-------------------------------------------------------

#SI III: Table S4

# MANTEL TESTS

mantel(N_2017, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(N_2017)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(N_2017))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(N_2019, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(N_2019)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(N_2019))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(N_2021, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(N_2021)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(N_2021))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))

mantel(L_2017, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(L_2017)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(L_2017))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(L_2019, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(L_2019)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(L_2019))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(L_2021, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(L_2021)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(L_2021))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))

mantel(M_2017, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(M_2017)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(M_2017))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(M_2019, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(M_2019)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(M_2019))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(M_2021, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(M_2021)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(M_2021))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))

mantel(H_2017, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(H_2017)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(H_2017))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(H_2019, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(H_2019)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(H_2019))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(H_2021, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(H_2021)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(H_2021))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))

mantel(P_2017, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(P_2017)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(P_2017))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))

### nothing significant

#######################################################################################
#######################################################################################

#By treatment and year lumped


y2017 <- Jsim[
  grepl("2017", rownames(Jsim)), 
  grepl("2017", colnames(Jsim))]
y2017 <- y2017[!grepl("P", rownames(y2017)),
               !grepl("P", colnames(y2017))]

y2019 <- Jsim[
  grepl("2019", rownames(Jsim)), 
  grepl("2019", colnames(Jsim))]

y2021 <- Jsim[
  grepl("2021", rownames(Jsim)), 
  grepl("2021", colnames(Jsim))]


mantel(y2017, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(y2017)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(y2017))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(y2019, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(y2019)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(y2019))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))
mantel(y2021, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(y2021)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(y2021))))], 
  method="pearson", permutations=999, strata = NULL,
  na.rm = FALSE, parallel = getOption("mc.cores"))


### Nothing significant

#-------------------------------------------

### MANTEL CORRELOGRAMS - JACCARD

ddd<-c()
#ddd<-c(ddd,0)
ddd<-c(ddd,60)
ddd<-c(ddd,120)
ddd<-c(ddd,180)
ddd<-c(ddd,240)
ddd<-c(ddd,300)
ddd<-c(ddd,400)
ddd<-c(ddd,500)
ddd<-c(ddd,700)
ddd<-c(ddd,900)
ddd<-c(ddd,1100)
ddd<-c(ddd,1300)
ddd<-c(ddd,1500)
ddd<-c(ddd,1700)

#SI III: Figure S2

#Mantel Correlograms
Jcorrelog_2017<-mantel.correlog(y2017, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(y2017)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(y2017))))],
  nperm=9999,cutoff=FALSE, break.pts=ddd)
plot(Jcorrelog_2017, alpha=0.05)

Jcorrelog_2019<-mantel.correlog(y2019, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(y2019)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(y2019))))],
  nperm=9999,cutoff=FALSE, break.pts=ddd)
plot(Jcorrelog_2019, alpha=0.05)

Jcorrelog_2021<-mantel.correlog(y2021, dist_matrix[
  sapply(rownames(dist_matrix), function(x) any(grepl(x, rownames(y2021)))),
  sapply(colnames(dist_matrix), function(x) any(grepl(x, colnames(y2021))))],
  nperm=9999,cutoff=FALSE, break.pts=ddd)
plot(Jcorrelog_2021, alpha=0.05)
