library("ggplot2")
library("dplyr")
library("vegan")
library(cluster)
#library(picante)
#library(FD)
library(nlme)
library(lme4)
library(car)
library(blmeco)
library(emmeans)
library(ggplot2)
library("stringr")
library("plyr")

set.seed(1)

trait_map <- read.table("mafun.txt", sep="\t", header=TRUE)
trait_map$size_group <- gsub(" ", "", trait_map$size_group)
trait_map$diet_breadth <- gsub(" ", "", trait_map$diet_breadth)
trait_map$burial_strategy <- gsub(" ", "", trait_map$burial_strategy)
trait_map$functional_group <- paste(trait_map$size_group, trait_map$diet_breadth, trait_map$burial_strategy, sep="_")
unique(trait_map$functional_group)

comm_data <- read.delim("DB_all.txt")
row.names(comm_data) <- paste(comm_data$site_id, comm_data$year, sep="_")
DB_species <- colnames(comm_data[-(1:5)])


comm_long <- tidyr::pivot_longer(comm_data, cols=DB_species,
                                 names_to="sp_code",
                                 values_to="frequency")
comm_long <- comm_long[comm_long$frequency!=0,]

comm_long <- left_join(comm_long[!(names(comm_long) %in% c("richness","abundance"))], trait_map[c("sp_code", "functional_group")], by = "sp_code")

trait_long <- aggregate(frequency~site_id+year+treatment+functional_group, comm_long[!(names(comm_long) %in% "sp_code")], sum)

trait_df <- data.frame(pivot_wider(trait_long, names_from="functional_group", values_from="frequency"))
row.names(trait_df) <- paste(trait_df$site_id, trait_df$year, sep="_")

functions <- names(trait_df[-(1:3)])

names <- colnames(trait_df)

trait_df <-rbind.fill(trait_df, comm_data[!(row.names(comm_data) %in% row.names(trait_df)),])

trait_df <- trait_df[names(trait_df) %in% names]
trait_df[is.na(trait_df)] <- 0

trait_df$treatment_year <- paste(trait_df$treatment, trait_df$year, sep="_")
row.names(trait_df) <- paste(trait_df$site_id, trait_df$year, sep="_")


rich_sum <- function(x) {sum(x>0)}
trait_df$richness <- apply(trait_df[names(trait_df) %in% functions],1,rich_sum)
trait_df$abundance <- apply(trait_df[names(trait_df) %in% functions],1,sum)

trait_df$treatment <- factor(trait_df$treatment, levels= c("N", "L", "M", "H", "P"))
trait_df$year <- factor(trait_df$year, levels= c("2017", "2019", "2021"))


#---------------------------------------------------------------------------------------------------------

f_traits <- c("tunneler", "roller", "unknown", "dweller", "small", "medium", "large", "generalist", "coprophagous", "necrophagous") 
metadata <- c("treatment", "year", "treatment_year", "site_id", "site_year")

#------------------------------------------------------------------------------------------------------------------------------------------------


trait_df_ordi <- trait_df[!!rowSums(trait_df[names(trait_df) %in% functions]),] #idk why this works but many other solutions didnt, drops rows with species abundances == 0
#!! is not a typo
row.names(trait_df_ordi) <- paste(trait_df_ordi$site_id, trait_df_ordi$year, sep="_")

p_dissim <- vegdist(trait_df_ordi[,names(trait_df_ordi) %in% functions], method="bray")
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
p_sim$treatment <- factor(p_sim$treatment, levels = c("N", "L", "M", "H"))
p_sim$year <- factor(p_sim$year, levels = c("2017", "2019", "2021"))

### creating jaccard similarity matrix
j_dissim <- vegdist(trait_df_ordi[,names(trait_df_ordi) %in% functions], method="jaccard")
j_dissim <- data.frame(as.matrix(j_dissim))
#write.csv(j_dissim, "test_matrix.csv")
j_dissim <- j_dissim[,grepl("P", names(j_dissim))] #just retaining columns that are primary forest controls
j_dissim <- j_dissim[!grepl("P", row.names(j_dissim)),] #just removing rows that are primary forest controls
j_sim <- 1 - j_dissim #making dissimilarity values into similarity values
j_sim$similarity <- apply(j_sim,1,mean)
j_sim$similarity <- rowMeans(j_sim)
j_sim$ln_similarity <- log10(j_sim$similarity)
j_sim$plot_year <- row.names(j_sim)
j_sim <- j_sim %>% tidyr::separate(plot_year, c("site_id", "year")) #separating row names into columns, automatically chooses non-alphanumeric character to split by
j_sim$treatment <- stringr::str_extract(j_sim$site_id, "[A-Z]+" ) 
j_sim$treatment <- factor(j_sim$treatment, levels= c("N", "L", "M", "H"))
j_sim$year <- factor(j_sim$year, levels= c("2017", "2019", "2021"))

p_dissim$dissimilarity <- apply(p_dissim,1,mean)
p_dissim$dissimilarity <- rowMeans(p_dissim)
p_dissim$ln_dissimilarity <- log10(p_dissim$dissimilarity)
p_dissim$plot_year <- row.names(p_dissim)
p_dissim <- p_dissim %>% tidyr::separate(plot_year, c("site_id", "year")) #separating row names into columns, automatically chooses non-alphanumeric character to split by
p_dissim$treatment <- stringr::str_extract(p_dissim$site_id, "[A-Z]+" ) 


#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################

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
coords <- read.csv("DB_analysis_2-10-24/github/metadata/osa_resto_plots_site_coords.csv", header=TRUE)

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

#SI III: Table S2


library(spdep)
# Assume your data frame is 'plot_data' with columns: richness, treatment, lat, lon.
coords2 <- cbind(coords[coords$treatment!="P",]$LongD, coords[coords$treatment!="P",]$LatD)

rownames(coords2) <- rownames(coords)[1:40]

#using distance-based network instead

threshold <- 150  # meters
nb_dist <- dnearneigh(coords2, 0, threshold)
lw_dist <- nb2listw(nb_dist, style="W", zero.policy=TRUE)
moran.mc(trait_df[trait_df$year=="2017" & trait_df$treatment!="P",]$richness, lw_dist, nsim = 999, alternative="greater", zero.policy=TRUE)
moran.mc(trait_df[trait_df$year=="2019",]$richness, lw_dist, nsim = 999, alternative="greater", zero.policy=TRUE)
moran.mc(trait_df[trait_df$year=="2021",]$richness, lw_dist, nsim = 999, alternative="greater", zero.policy=TRUE)


#further parsing by year and treatment, Moran's I

coords_N <- cbind(coords[coords$treatment=="N",]$LongD, coords[coords$treatment=="N",]$LatD)
nb_N <- dnearneigh(coords_N, 0, threshold)
list_N <- nb2listw(nb_N, style="W", zero.policy=TRUE)
moran_N2017 <- moran.mc(trait_df[trait_df$year=="2017" & trait_df$treatment=="N",]$richness, list_N, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_N2017)
moran_N2019 <- moran.mc(trait_df[trait_df$year=="2019" & trait_df$treatment=="N",]$richness, list_N, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_N2019) ### SIGNIFICANT
moran_N2021 <- moran.mc(trait_df[trait_df$year=="2021" & trait_df$treatment=="N",]$richness, list_N, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_N2021)

coords_L <- cbind(coords[coords$treatment=="L",]$LongD, coords[coords$treatment=="L",]$LatD)
nb_L <- dnearneigh(coords_L, 0, threshold)
list_L <- nb2listw(nb_L, style="W", zero.policy=TRUE)
moran_L2017 <- moran.mc(trait_df[trait_df$year=="2017" & trait_df$treatment=="L",]$richness, list_L, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_L2017)
moran_L2019 <- moran.mc(trait_df[trait_df$year=="2019" & trait_df$treatment=="L",]$richness, list_L, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_L2019) ### SIGLIFICALT
moran_L2021 <- moran.mc(trait_df[trait_df$year=="2021" & trait_df$treatment=="L",]$richness, list_L, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_L2021)

coords_M <- cbind(coords[coords$treatment=="M",]$LongD, coords[coords$treatment=="M",]$LatD)
nb_M <- dnearneigh(coords_M, 0, threshold)
list_M <- nb2listw(nb_M, style="W", zero.policy=TRUE)
moran_M2017 <- moran.mc(trait_df[trait_df$year=="2017" & trait_df$treatment=="M",]$richness, list_M, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_M2017)
moran_M2019 <- moran.mc(trait_df[trait_df$year=="2019" & trait_df$treatment=="M",]$richness, list_M, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_M2019) ### SIGMIFICAMT
moran_M2021 <- moran.mc(trait_df[trait_df$year=="2021" & trait_df$treatment=="M",]$richness, list_M, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_M2021)

coords_H <- cbind(coords[coords$treatment=="H",]$LongD, coords[coords$treatment=="H",]$LatD)
nb_H <- dnearneigh(coords_H, 0, threshold)
list_H <- nb2listw(nb_H, style="W", zero.policy=TRUE)
moran_H2017 <- moran.mc(trait_df[trait_df$year=="2017" & trait_df$treatment=="H",]$richness, list_H, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_H2017)
moran_H2019 <- moran.mc(trait_df[trait_df$year=="2019" & trait_df$treatment=="H",]$richness, list_H, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_H2019) ### SIGHIFICAHT
moran_H2021 <- moran.mc(trait_df[trait_df$year=="2021" & trait_df$treatment=="H",]$richness, list_H, nsim = 999, alternative="greater", zero.policy=TRUE)
print(moran_H2021)


#==============================================================================================================================================



# SETTING UP MANTEL TESTS
### BRAy_CURTIS


BCdissim <- vegdist(trait_df_ordi[,names(trait_df_ordi) %in% functions], method="bray")

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

# MANTEL TESTS

#SI III: Table S5

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

#nothing significant

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

#SI III: Figure S3

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

Jdissim <- vegdist(trait_df_ordi[,names(trait_df_ordi) %in% functions], method="jaccard")
Jdissim <- data.frame(as.matrix(Jdissim))
Jsim <- 1 - Jdissim #making dissimilarity values into similarity values
#BCsim_2017$site_id <- row.names(BCsim_2017)
#BCsim_2017$treatment <- stringr::str_extract(BCsim_2017$site_id, "[A-Z]+" ) 
#BCsim_2017$treatment <- factor(BCsim_2017$treatment, levels= c("N", "L", "M", "H"))


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


#----------------------------------------------

#----------------------------------------------

# MANTEL TESTS

#SI III: Table S6

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

#SI III: Figure S4

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

