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
library("performance")


#-----------------------------------------------------------------------------------------------------------

set.seed(1)

#-----------------------------------------------------------------------------------------------------------

veg <- read.csv("vegetation_structure_2-24-24.csv", header=TRUE)
veg$plot_year <- paste(veg$site_ID, veg$year, sep="_")
veg$lealit <- as.numeric(veg$lealit)
veg$leadep <- as.numeric(veg$leadep)

soil <- read.csv("soil_data_2-24-24.csv", header=TRUE)
soil$plot_year <- paste(soil$Plot, soil$Year, sep="_")
soil$CE <- as.numeric(soil$CE)
soil$P <- as.numeric(soil$P)
#"Copy of Monitoreo de crecimiento.csv"

#----------------------------------------------------------------------------------------------------------

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
comm_data$plot_year <- paste(comm_data$site_id, comm_data$year, sep="_")

#----------------------------------------------------------------------------------------------------------

#generating similarity data
comm_data_ordi <- comm_data[!!rowSums(comm_data[names(comm_data) %in% DB_species]),] #idk why this works but many other solutions didnt, drops rows with species abundances == 0
#!! is not a typo
row.names(comm_data_ordi) <- paste(comm_data_ordi$site_id, comm_data_ordi$year, sep="_")

pt_dissim <- vegdist(comm_data_ordi[,names(comm_data_ordi) %in% DB_species], method="bray")
pt_dissim <- data.frame(as.matrix(pt_dissim))
#write.csv(pt_dissim, "test_matrix.csv")
pt_dissim <- pt_dissim[,grepl("P", names(pt_dissim))] #just retaining columns that are primary forest controls
pt_dissim <- pt_dissim[!grepl("P", row.names(pt_dissim)),] #just removing rows that are primary forest controls
pt_sim <- 1 - pt_dissim #making dissimilarity values into similarity values
pt_sim$similarity <- apply(pt_sim,1,mean)
pt_sim$similarity <- rowMeans(pt_sim)
pt_sim$ln_similarity <- log10(pt_sim$similarity)
pt_sim$plot_year <- row.names(pt_sim)
pt_sim <- pt_sim %>% tidyr::separate(plot_year, c("site_id", "year")) #separating row names into columns, automatically chooses non-alphanumeric character to split by
pt_sim$treatment <- stringr::str_extract(pt_sim$site_id, "[A-Z]+" ) 
pt_sim$treatment <- factor(pt_sim$treatment, levels= c("N", "L", "M", "H"))
pt_sim$year <- factor(pt_sim$year, levels= c("2017", "2019", "2021"))
pt_sim$plot_year <- row.names(pt_sim)

### creating jaccard similarity matrix
jt_dissim <- vegdist(comm_data_ordi[,names(comm_data_ordi) %in% DB_species], method="jaccard")
jt_dissim <- data.frame(as.matrix(jt_dissim))
#write.csv(jt_dissim, "test_matrix.csv")
jt_dissim <- jt_dissim[,grepl("P", names(jt_dissim))] #just retaining columns that are primary forest controls
jt_dissim <- jt_dissim[!grepl("P", row.names(jt_dissim)),] #just removing rows that are primary forest controls
jt_sim <- 1 - jt_dissim #making dissimilarity values into similarity values
jt_sim$jt_similarity <- apply(jt_sim,1,mean)
jt_sim$jt_similarity <- rowMeans(jt_sim)
jt_sim$ln_similarity <- log10(jt_sim$jt_similarity)
jt_sim$plot_year <- row.names(jt_sim)
jt_sim <- jt_sim %>% tidyr::separate(plot_year, c("site_id", "year")) #separating row names into columns, automatically chooses non-alphanumeric character to split by
jt_sim$treatment <- stringr::str_extract(jt_sim$site_id, "[A-Z]+" ) 
jt_sim$treatment <- factor(jt_sim$treatment, levels= c("N", "L", "M", "H"))
jt_sim$year <- factor(jt_sim$year, levels=c("2017","2019","2021"))
jt_sim$plot_year <- row.names(jt_sim)

#----------------------------------------------------------------------------------------------------------

### putting data together

#attaching vegetation to similarity data
comm_veg <- left_join(comm_data, veg[!(names(veg) %in% c("treatment", "year", "site_ID", "Abund", "R"))], by="plot_year")

comm_soil <- left_join(comm_data, soil[!(names(soil) %in% c("Treatment", "Year", "Plot"))], by="plot_year")

comm_env <- left_join(comm_veg, comm_soil[-(1:35)], by="plot_year")
colnames(comm_env)

comm_var <- c(37:45, 48:63) #selecting predictor variables

comm_env[comm_var] <- lapply(comm_env[comm_var], as.numeric)

zcomm_env <- comm_env #copying df

#rescaling (z-score) predictors, retaining dataframe structure
zcomm_env[comm_var] <- as.data.frame(lapply(zcomm_env[comm_var], function(x) as.numeric(scale(x))))


#attaching vegetation to similarity data
pt_sim <- left_join(pt_sim, jt_sim[names(jt_sim) %in% c('jt_similarity', 'plot_year')], by="plot_year")

pt_sim_veg <- left_join(pt_sim, veg[!(names(veg) %in% c("treatment", "year", "site_ID", "Abund", "R"))], by="plot_year")

pt_sim_soil <- left_join(pt_sim, soil[!(names(soil) %in% c("Treatment", "Year", "Plot"))], by="plot_year")

pt_sim_env <- left_join(pt_sim_veg, pt_sim_soil[-c(1:19, 21)], by="plot_year")
colnames(pt_sim_env)

pt_var <- c(22:30, 33:48) #selecting predictor variables

pt_sim_env[pt_var] <- lapply(pt_sim_env[pt_var], as.numeric)

ptz_sim_env <- pt_sim_env #copying df

#rescaling (z-score) predictors, retaining dataframe structure
ptz_sim_env[pt_var] <- as.data.frame(lapply(ptz_sim_env[pt_var], function(x) as.numeric(scale(x))))

test <- na.omit(ptz_sim_env)
#rerun everything w/o na.omit???

#-----------------------------------------------------------------------------------------------------------

# BRAY-CURTIS SIMILARITY

#CHECKING COLLINEARITY

all_bct <- glmmTMB(similarity ~ pH+Acidity+Ca+Mg+K+P+SA_percent+Zn+Cu+Fe+Mn+percent_C+percent_N+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf,
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = ptz_sim_env)
check_collinearity(all_bct) #had to remove CICE

all_bct <- glmmTMB(similarity ~ pH+Ca+Mg+K+P+SA_percent+Mn+percent_N+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf,
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = ptz_sim_env)
check_collinearity(all_bct) #had to remove CICE
#removed %C, Acidity, Zn, Fe, Cu

#--------------------------------

###   step AIC magic: Bray-Curtis similarity

all_bct <- glmmTMB(similarity ~ pH+Ca+Mg+K+P+SA_percent+Mn+percent_N+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf,
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = na.omit(ptz_sim_env))
stepAIC(all_bct, direction="backward")

#--------------------------------

#best stepAIC model
best_bct <- glmmTMB(similarity ~  Ca + P + Mn + CN_ratio + canhei + grass + leadep + grahei + (1 |site_id), 
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = ptz_sim_env)

#best model!
best_bct2 <- glmmTMB(similarity ~  Ca + P + Mn + CN_ratio + canhei + grass + leadep + grahei, 
                    family = ordbeta(),
                    start=list(psi = c(-1, 1)),
                    data = ptz_sim_env)
#SI I: Table S18
summary(best_bct2)
histogram(residuals(best_bct2))
plot(residuals(best_bct2))
testZeroInflation(best_bct2)
testDispersion(best_bct2)
testOverdispersion(best_bct2)
plot(simulateResiduals(best_bct2, n=500))

AIC(best_bct, best_bct2)

#-----------------------------------------------------------------------------------------------------------

# JACCARD SIMILARITY


#CHECKING COLLINEARITY

all_jt <- glmmTMB(jt_similarity ~ pH+Acidity+Ca+Mg+P+K+SA_percent+CICE+Zn+Cu+Fe+Mn+percent_C+percent_N+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf,
                 family = ordbeta(),
                 start=list(psi = c(-1, 1)),
                 data = ptz_sim_env)
check_collinearity(all_jt)

all_jt <- glmmTMB(jt_similarity ~ pH+Ca+Mg+P+K+SA_percent+Mn+percent_N+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf,
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = ptz_sim_env)
check_collinearity(all_jt)
#removed: CICE, %C, Acidity, Zn, Fe, Cu


#--------------------------------

###   step AIC magic: Jaccard similarity

all_jt <- glmmTMB(jt_similarity ~ pH+Ca+Mg+P+K+SA_percent+Mn+percent_N+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf,
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = na.omit(ptz_sim_env))
stepAIC(all_jt, direction="backward")

#---------------------------------

#best stepAIC model
best_jt <- glmmTMB(jt_similarity ~  Ca + P + Mn + CN_ratio + canhei + grass + leadep +  grahei + (1 |site_id), 
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = ptz_sim_env)#,
                  #control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))) #optimizer argument to avoid warning

best_jt2 <- glmmTMB(jt_similarity ~ Ca + P + Mn + CN_ratio + canhei + grass + leadep +  grahei, #BEST MODEL FIT
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = ptz_sim_env)
#SI I: Table S18
summary(best_jt2)
histogram(residuals(best_jt2))
plot(residuals(best_jt2))
testZeroInflation(best_jt2)
testDispersion(best_jt2)
testOverdispersion(best_jt2)
plot(simulateResiduals(best_jt2, n=500))

AIC(best_jt, best_jt2)


#-----------------------------------------------------------------------------------------------------------

# ABSOLUTE ABUNDANCE

#zcomm_env <- zcomm_env[!(names(zcomm_env) %in% c("CE", "X", "X.1"))]

#checking collinearity
#already removed CE, P, X, X.1 (missing data)

all_abun <- glmmTMB(abundance ~ pH+Acidity+Ca+Mg+K+CICE+P+SA_percent+Zn+Cu+Fe+Mn+percent_C+percent_N+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf,
                    family = genpois(link="log"),
                    data = zcomm_env[zcomm_env$treatment!="P",])
check_collinearity(all_abun)

all_abun <- glmmTMB(abundance ~ Ca+Mg+K+P+SA_percent+Mn+percent_N+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf,
                    family = genpois(link="log"),
                    data = zcomm_env[zcomm_env$treatment!="P",])
check_collinearity(all_abun)
#removed CICE, %C, Acidity, Zn, Fe, pH, Cu

#--------------------------------

###   step AIC magic: Absolute Abundance


all_abun <- glmmTMB(abundance ~ Ca+Mg+K+P+SA_percent+Mn+percent_N+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf,
                    family = genpois(link="log"),
                    data = na.omit(zcomm_env[zcomm_env$treatment!="P",]))
stepAIC(all_abun, direction="backward")

#--------------------------------

#best stepAIC model
all_abun <- glmmTMB(abundance ~ Mg + P + CN_ratio + canhei + bargro + lealit + leadep + (1 | site_id), #2nd best model
                    family = genpois(link="log"),
                    data = zcomm_env[zcomm_env$treatment!="P",])
summary(all_abun)

all_abun2 <- glmmTMB(abundance ~ Mg + P + CN_ratio + canhei + bargro + lealit + leadep, #BEST FIT MODEL!!! kind of, not by 2 AIC
                     family = genpois(link="log"),
                     data = zcomm_env[zcomm_env$treatment!="P",])
#SI I: Table S17
summary(all_abun2)
histogram(residuals(all_abun2))
plot(residuals(all_abun2))
testZeroInflation(all_abun2)
testDispersion(all_abun2)
testOverdispersion(all_abun2)
plot(simulateResiduals(all_abun2, n=500))

all_abun3 <- glmmTMB(abundance ~ Mg + P + CN_ratio + canhei + bargro + lealit + leadep + distf + (1 | site_id),
                    family = genpois(link="log"),
                    data = zcomm_env[zcomm_env$treatment!="P",])

all_abun4 <- glmmTMB(abundance ~ Mg + P + CN_ratio + canhei + bargro + lealit + leadep + distf, #3rd best model
                    family = genpois(link="log"),
                    data = zcomm_env[zcomm_env$treatment!="P",])
summary(all_abun4)

AIC(all_abun, all_abun2, all_abun3, all_abun4)


#-----------------------------------------------------------------------------------------------------------


# ABSOLUTE RICHNESS

#zcomm_env <- zcomm_env[!(names(zcomm_env) %in% c("CE", "X", "X.1"))]

#checking collinearity

all_richt <- glmmTMB(richness ~ pH+Acidity+Ca+Mg+P+K+CICE+SA_percent+Zn+Cu+Fe+Mn+percent_C+percent_N+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf,
                     family = genpois(link="log"),
                     data = zcomm_env[zcomm_env$treatment!="P",])
check_collinearity(all_richt)

all_richt <- glmmTMB(richness ~ pH+Ca+P+K+SA_percent+Mn+percent_N+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf,
                     family = genpois(link="log"),
                     data = zcomm_env[zcomm_env$treatment!="P",])
check_collinearity(all_richt)
#removed CICE, %C, Acidity, Zn, Cu, Mg

#--------------------------------

###   step AIC magic: Absolute Richness

all_richt <- glmmTMB(richness ~ pH+Ca+P+K+SA_percent+Mn+percent_N+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf,
                     family = genpois(link="log"),
                     data = na.omit(zcomm_env[zcomm_env$treatment!="P",]))
stepAIC(all_richt, direction="backward") 

#--------------------------------

#best stepAIC model
best_richt <- glmmTMB(richness ~  P + Mn + percent_N + CN_ratio + canhei + (1 |site_id),
                     family = genpois(link="log"),
                     data = zcomm_env[zcomm_env$treatment!="P ",],
                     control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))) #optimizer argument to avoid warning

best_richt2 <- glmmTMB(richness ~ P + Mn + percent_N + CN_ratio + canhei, #BEST MODEL *** 
                      family = genpois(link="log"),
                      data = zcomm_env[zcomm_env$treatment!="P",])
#SI I: Table S17
summary(best_richt2)
histogram(residuals(best_richt2))
plot(residuals(best_richt2))
testZeroInflation(best_richt2)
testDispersion(best_richt2)
plot(simulateResiduals(best_richt2, n=500))

AIC(best_richt, best_richt2) 


#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------


############################# FUNCTIONAL TRAIT DATA PREP


#creating functional trait data

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
trait_df$plot_year <- paste(trait_df$site_id, trait_df$year, sep="_")


#----------------------------------------------------------------------------------------------------------

trait_df_ordi <- trait_df[!!rowSums(trait_df[names(trait_df) %in% functions]),] #idk why this works but many other solutions didnt, drops rows with species abundances == 0
#!! is not a typo
row.names(trait_df_ordi) <- paste(trait_df_ordi$site_id, trait_df_ordi$year, sep="_")

pf_dissim <- vegdist(trait_df_ordi[,names(trait_df_ordi) %in% functions], method="bray")
pf_dissim <- data.frame(as.matrix(pf_dissim))
#write.csv(pf_dissim, "test_matrix.csv")
pf_dissim <- pf_dissim[,grepl("P", names(pf_dissim))] #just retaining columns that are primary forest controls
pf_dissim <- pf_dissim[!grepl("P", row.names(pf_dissim)),] #just removing rows that are primary forest controls
pf_sim <- 1 - pf_dissim #making dissimilarity values into similarity values
pf_sim$similarity <- apply(pf_sim,1,mean)
pf_sim$similarity <- rowMeans(pf_sim)
pf_sim$ln_similarity <- log10(pf_sim$similarity)
pf_sim$plot_year <- row.names(pf_sim)
pf_sim <- pf_sim %>% tidyr::separate(plot_year, c("site_id", "year")) #separating row names into columns, automatically chooses non-alphanumeric character to split by
pf_sim$treatment <- stringr::str_extract(pf_sim$site_id, "[A-Z]+" ) 
pf_sim$treatment <- factor(pf_sim$treatment, levels = c("N", "L", "M", "H"))
pf_sim$year <- factor(pf_sim$year, levels = c("2017", "2019", "2021"))
pf_sim$plot_year <- row.names(pf_sim)

### creating jaccard similarity matrix
jf_dissim <- vegdist(trait_df_ordi[,names(trait_df_ordi) %in% functions], method="jaccard")
jf_dissim <- data.frame(as.matrix(jf_dissim))
#write.csv(jf_dissim, "test_matrix.csv")
jf_dissim <- jf_dissim[,grepl("P", names(jf_dissim))] #just retaining columns that are primary forest controls
jf_dissim <- jf_dissim[!grepl("P", row.names(jf_dissim)),] #just removing rows that are primary forest controls
jf_sim <- 1 - jf_dissim #making dissimilarity values into similarity values
jf_sim$jf_similarity <- apply(jf_sim,1,mean)
jf_sim$jf_similarity <- rowMeans(jf_sim)
jf_sim$ln_similarity <- log10(jf_sim$jf_similarity)
jf_sim$plot_year <- row.names(jf_sim)
jf_sim <- jf_sim %>% tidyr::separate(plot_year, c("site_id", "year")) #separating row names into columns, automatically chooses non-alphanumeric character to split by
jf_sim$treatment <- stringr::str_extract(jf_sim$site_id, "[A-Z]+" ) 
jf_sim$treatment <- factor(jf_sim$treatment, levels= c("N", "L", "M", "H"))
jf_sim$year <- factor(jf_sim$year, levels= c("2017", "2019", "2021"))
jf_sim$plot_year <- row.names(jf_sim)

#----------------------------------------------------------------------------------------------------------

env_var <- c("cancov","canhei","bargro","lealit","grass","veg","leadep","grahei","distf", 
             "pH","Acidity","Ca","Mg","K","CICE","SA_percent","P","Zn","Cu","Fe","Mn","CE","percent_C","percent_N","CN_ratio")

### putting datasets together!

#attaching vegetation to similarity data
trait_veg <- left_join(trait_df, veg[!(names(veg) %in% c("treatment", "year", "site_ID", "Abund", "R"))], by="plot_year")

trait_soil <- left_join(trait_df, soil[!(names(soil) %in% c("Treatment", "Year", "Plot"))], by="plot_year")

trait_env <- left_join(trait_veg, trait_soil[-(1:17)], by="plot_year")
colnames(trait_env)

trait_var <- c(19:27, 30:45) #selecting predictor variables

trait_env[trait_var] <- lapply(trait_env[trait_var], as.numeric)

ztrait_env <- trait_env #copying df

#rescaling (z-score) predictors, retaining dataframe format
ztrait_env[trait_var] <- as.data.frame(lapply(ztrait_env[trait_var], function(x) as.numeric(scale(x))))



#attaching vegetation to similarity data
pf_sim <- left_join(pf_sim, jf_sim[names(jf_sim) %in% c('jf_similarity', 'plot_year')], by="plot_year")

pf_sim_veg <- left_join(pf_sim, veg[!(names(veg) %in% c("treatment", "year", "site_ID", "Abund", "R"))], by="plot_year")

pf_sim_soil <- left_join(pf_sim, soil[!(names(soil) %in% c("Treatment", "Year", "Plot"))], by="plot_year")

pf_sim_env <- left_join(pf_sim_veg, pf_sim_soil[-c(1:19, 21)], by="plot_year")
colnames(pf_sim_env)

pf_var <- c(22:30, 33:48) #selecting predictor variables

pf_sim_env[pf_var] <- lapply(pf_sim_env[pf_var], as.numeric)

pfz_sim_env <- pf_sim_env #copying df

#rescaling (z-score) predictors, retaining dataframe format
pfz_sim_env[pf_var] <- as.data.frame(lapply(pfz_sim_env[pf_var], function(x) as.numeric(scale(x))))


#-----------------------------------------------------------------------------------------------------------

#set.seed(1)

#-----------------------------------------------------------------------------------------------------------

# BRAY-CURTIS SIMILARITY


#CHECKING COLLINEARITY
all_bcf <- glmmTMB(similarity ~ pH+CICE+Acidity+Ca+Mg+K+P+SA_percent+Zn+Cu+Fe+Mn+percent_N+percent_C+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf + (1 |site_id),
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = pfz_sim_env)
check_collinearity(all_bcf)

all_bcf <- glmmTMB(similarity ~ pH+Ca+K+P+SA_percent+Cu+Mn+percent_N+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf + (1 |site_id),
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = pfz_sim_env)
check_collinearity(all_bcf)
#removed: CICE, %C, Acidity, Zn, Mg, Fe

#--------------------------------

###   step AIC magic: Bray-Curtis similarity

all_bcf <- glmmTMB(similarity ~ pH+Ca+K+P+SA_percent+Cu+Mn+percent_N+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf + (1 |site_id),
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = na.omit(pfz_sim_env))
stepAIC(all_bcf, direction="backward")

#----------------------------------

best_bcf <- glmmTMB(similarity ~   P+Mn+CN_ratio+canhei+lealit+grahei+distf  + (1 |site_id), 
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = pfz_sim_env)

best_bcf2 <- glmmTMB(similarity ~  P+Mn+CN_ratio+canhei+lealit+grahei+distf, #BEST MODEL FIT
                    family = ordbeta(),
                    start=list(psi = c(-1, 1)),
                    data = pfz_sim_env)
#SI I: Table S18
summary(best_bcf2)

AIC(best_bcf, best_bcf2)


#-----------------------------------------------------------------------------------------------------------

# JACCARD SIMILARITY


#CHECKING COLLINEARITY

all_jf <- glmmTMB(jf_similarity ~ pH+CICE+Acidity+Ca+Mg+K+P+SA_percent+Zn+Cu+Fe+Mn+percent_N+percent_C+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf + (1 |site_id),
                 family = ordbeta(),
                 start=list(psi = c(-1, 1)),
                 data = pfz_sim_env)
check_collinearity(all_jf)

all_jf <- glmmTMB(jf_similarity ~ pH+Ca+K+P+SA_percent+Cu+Mn+percent_N+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf + (1 |site_id),
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = pfz_sim_env)
check_collinearity(all_jf)
#CICE, %C, Acidity, Zn, Fe, Mg


#--------------------------------

###   step AIC magic: Jaccard similarity

all_jf <- glmmTMB(jf_similarity ~ pH+Ca+K+P+SA_percent+Cu+Mn+percent_N+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf + (1 |site_id),
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = na.omit(pfz_sim_env))
stepAIC(all_jf, direction="backward")


#--------------------------------

#best stepAIC model
best_jf <- glmmTMB(jf_similarity ~ Ca+P+Mn+CN_ratio+canhei+lealit+grahei+distf+ (1 |site_id), #USE THIS ONE
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = pfz_sim_env)


best_jf2 <- glmmTMB(jf_similarity ~ Ca+P+Mn+CN_ratio+canhei+lealit+grahei+distf, #better fit but convergence problems
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = pfz_sim_env)
#SI I: Table S18
summary(best_jf2)


AIC(best_jf, best_jf2)

#-----------------------------------------------------------------------------------------------------------


# ABSOLUTE RICHNESS


all_richf <- glmmTMB(richness ~ pH+Acidity+Ca+Mg+K+P+CICE+SA_percent+Zn+Cu+Fe+Mn+percent_C+percent_N+CN_ratio+cancov+canhei+bargro+lealit+grass+veg+leadep+grahei+distf + (1 |site_id),
                    family = genpois(link="log"),
                    data =  ztrait_env[ztrait_env$treatment!="P",])
check_collinearity(all_richf) 

all_richf <- glmmTMB(richness ~ pH+Ca+K+P+SA_percent+Mn+percent_N+CN_ratio+cancov+bargro+lealit+grass+veg+leadep+grahei+distf + (1 |site_id),
                    family = genpois(link="log"),
                    data = ztrait_env[ztrait_env$treatment!="P",])
check_collinearity(all_rich) 
#CICE, %C, Acidity, Zn, Fe, Cu, Mg, canhei

#--------------------------------

###   step AIC magic: Absolute Richness

all_richf <- glmmTMB(richness ~ pH+Ca+K+P+SA_percent+Mn+percent_N+CN_ratio+cancov+bargro+lealit+grass+veg+leadep+grahei+distf + (1 |site_id),
                    family = genpois(link="log"),
                    data = na.omit(ztrait_env[ztrait_env$treatment!="P",]))
stepAIC(all_richf, direction="backward")

#--------------------------------

#best stepAIC model
all_richf <- glmmTMB(richness ~ Ca+SA_percent+percent_N+CN_ratio+cancov + (1 | site_id),
                    family = genpois(link="log"),
                    data = ztrait_env[ztrait_env$treatment!="P",])

#weird warning????
all_richf2 <- glmmTMB(richness ~ Ca+SA_percent+percent_N+CN_ratio+cancov, #BEST FIT MODEL!!!
                     family = genpois(link="log"),
                     data = ztrait_env[ztrait_env$treatment!="P",],
                     control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
#SI I: Table S17
summary(all_richf2)

AIC(all_richf, all_richf2)


#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------


#################### CREATING NMDS PLOTS

### BRAY-CURTIS

### COMM DATA vector generation

env_signt <- c("CN_ratio", "canhei", "grass", "leadep", "grahei")

#make community and environmental data nmds object

comm_env_ordi <- zcomm_env[!!rowSums(zcomm_env[names(zcomm_env) %in% DB_species]),]

comm_mds <- metaMDS(comm_env_ordi[names(comm_env_ordi) %in% DB_species], distance = "bray", autotransform = FALSE, trymax=100)
comm_scrs <- as.data.frame(scores(comm_mds, display = "site"))
comm_scrs <- cbind(comm_scrs, site_id = comm_env_ordi$site_id)
comm_scrs <- cbind(comm_scrs, year = comm_env_ordi$year)
comm_scrs <- cbind(comm_scrs, treatment = comm_env_ordi$treatment)
comm_scrs <- cbind(comm_scrs, treatment_year = comm_env_ordi$treatment_year)

comm_envfit <- envfit(comm_mds, comm_env_ordi[names(comm_env_ordi) %in% env_var], permutations = 999, na.rm=TRUE) # this fits environmental vectors to env variables in comm_env_ordi
#unclear if i dropped necessary columns for function df
env_scorest <- as.data.frame(scores(comm_envfit, display = "vectors")) #extracts relevant scores from envifit
env_scorest <- cbind(env_scorest, env.variables = rownames(env_scorest)) #and then gives them their names
comm_env_scores <- env_scorest[env_scorest$env.variables %in% env_signt,] #subsetting by significant GLM predictors

#maskes plotting primary forest control easier for species data
comm_scrs$year <- factor(comm_scrs$year, levels = c("2017", "2019", "2021", "2100"))
comm_scrs[comm_scrs$treatment=="P",]$year <- "2100"



### FUNCTIONAL TRAITS vector generation

env_signf <- c("P", "Mn", "CN_ratio", "canhei", "lealit")

#make community and environmental data nmds object

trait_env_ordi <- ztrait_env[!!rowSums(ztrait_env[names(ztrait_env) %in% functions]),]

trait_mds <- metaMDS(trait_env_ordi[names(trait_env_ordi) %in% functions], distance = "bray", trymax=100, autotransform = FALSE)
trait_scrs <- as.data.frame(scores(trait_mds, display = "site"))
trait_scrs <- cbind(trait_scrs, site_id = trait_env_ordi$site_id)
trait_scrs <- cbind(trait_scrs, year = trait_env_ordi$year)
trait_scrs <- cbind(trait_scrs, treatment = trait_env_ordi$treatment)
trait_scrs <- cbind(trait_scrs, treatment_year = trait_env_ordi$treatment_year)

trait_envfit <- envfit(trait_mds, trait_env_ordi[names(trait_env_ordi) %in% env_var], permutations = 999, na.rm=TRUE) # this fits environmental vectors to env variables in trait_env_ordi

env_scoresf <- as.data.frame(scores(trait_envfit, display = "vectors")) #extracts relevant scores from envifit
env_scoresf <- cbind(env_scoresf, env.variables = rownames(env_scoresf)) #and then gives them their names
trait_env_scores <- env_scoresf[env_scoresf$env.variables %in% env_signf,] #subsetting by significant GLM predictors

trait_scrs$NMDS1 <-  trait_scrs$NMDS1*-1 #flip plot x axis, so that it views similar to species data
trait_scrs$year <- factor(trait_scrs$year, levels = c("2017", "2019", "2021", "2100"))
trait_scrs[trait_scrs$treatment=="P",]$year <- "2100"

#--------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------

### JACCARD

### COMM DATA vector generation

jenv_signt <- c("CN_ratio", "canhei", "grass", "leadep", "grahei")

jzcomm_env <- zcomm_env
jptz_sim_env <- ptz_sim_env

#make community and environmental data nmds object

jcomm_env_ordi <- jzcomm_env[!!rowSums(jzcomm_env[names(jzcomm_env) %in% DB_species]),]

jcomm_mds <- metaMDS(jcomm_env_ordi[names(jcomm_env_ordi) %in% DB_species], distance = "jaccard", autotransform = FALSE, trymax=100)
jcomm_scrs <- as.data.frame(scores(jcomm_mds, display = "site"))
jcomm_scrs <- cbind(jcomm_scrs, site_id = jcomm_env_ordi$site_id)
jcomm_scrs <- cbind(jcomm_scrs, year = jcomm_env_ordi$year)
jcomm_scrs <- cbind(jcomm_scrs, treatment = jcomm_env_ordi$treatment)
jcomm_scrs <- cbind(jcomm_scrs, treatment_year = jcomm_env_ordi$treatment_year)

jcomm_envfit <- envfit(jcomm_mds, jcomm_env_ordi[names(jcomm_env_ordi) %in% env_var], permutations = 999, na.rm=TRUE) # this fits environmental vectors to env variables in comm_env_ordi
#unclear if i dropped necessary columns for function df
jenv_scorest <- as.data.frame(scores(jcomm_envfit, display = "vectors")) #extracts relevant scores from envifit
jenv_scorest <- cbind(jenv_scorest, env.variables = rownames(jenv_scorest)) #and then gives them their names
jcomm_env_scores <- jenv_scorest[jenv_scorest$env.variables %in% jenv_signt,] #subsetting by significant GLM predictors

#maskes plotting primary forest control easier for species data
jcomm_scrs$NMDS2 <-  jcomm_scrs$NMDS2*-1 #flip plot x axis, so that it views similar to species data
jcomm_scrs$year <- factor(jcomm_scrs$year, levels = c("2017", "2019", "2021", "2100"))
jcomm_scrs[jcomm_scrs$treatment=="P",]$year <- "2100"


### FUNCTIONAL TRAITS vector generation

jenv_signf <- c("P", "CN_ratio", "canhei", "lealit")


jztrait_env <- ztrait_env
jpfz_sim_env <- pfz_sim_env

jtrait_env_ordi <- jztrait_env[!!rowSums(jztrait_env[names(jztrait_env) %in% functions]),]

jtrait_mds <- metaMDS(jtrait_env_ordi[names(jtrait_env_ordi) %in% functions], distance = "jaccard", trymax=100, autotransform = FALSE)
jtrait_scrs <- as.data.frame(scores(jtrait_mds, display = "site"))
jtrait_scrs <- cbind(jtrait_scrs, site_id = jtrait_env_ordi$site_id)
jtrait_scrs <- cbind(jtrait_scrs, year = jtrait_env_ordi$year)
jtrait_scrs <- cbind(jtrait_scrs, treatment = jtrait_env_ordi$treatment)
jtrait_scrs <- cbind(jtrait_scrs, treatment_year = jtrait_env_ordi$treatment_year)

jtrait_envfit <- envfit(jtrait_mds, jtrait_env_ordi[names(jtrait_env_ordi) %in% env_var], permutations = 999, na.rm=TRUE) # this fits environmental vectors to env variables in trait_env_ordi

jenv_scoresf <- as.data.frame(scores(jtrait_envfit, display = "vectors")) #extracts relevant scores from envifit
jenv_scoresf <- cbind(jenv_scoresf, env.variables = rownames(jenv_scoresf)) #and then gives them their names
jtrait_env_scores <- jenv_scoresf[jenv_scoresf$env.variables %in% jenv_signf,] #subsetting by significant GLM predictors

jtrait_scrs$NMDS1 <-  jtrait_scrs$NMDS1*-1 #flip plot x axis, so that it views similar to species data
jtrait_scrs$year <- factor(jtrait_scrs$year, levels = c("2017", "2019", "2021", "2100"))
jtrait_scrs[jtrait_scrs$treatment=="P",]$year <- "2100"

#-----------------------------------------------------------------------------------------------------------------------

############ PLOTTING NMDS FIGURES

#MAIN FIGURE 4

library("ggrepel")

#making full panel figure

A <- ggplot()+ ggtitle("A)")+#sets up the plot
  geom_point(data=comm_scrs, aes(x=NMDS1, y=NMDS2, colour=year, shape=year), size = 10)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  #coord_fixed()+
  #labs(colour = "treatment", shape = "year")+ # add legend labels for Management and Landuse
  theme_bw(50)+
  theme(legend.position = "none") +
  scale_shape_manual(breaks = c("2017","2019","2021","2100"), values=c("circle", "circle","circle","diamond")) +
  scale_color_manual(breaks = c("2017","2019","2021","2100"), values=c("darkorange", "yellow3","springgreen2","forestgreen")) +
  geom_segment(data = comm_env_scores, aes(x = 0, xend=3*NMDS1, y=0, yend=3*NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", size=5, linewidth=2) #+ #add vector arrows of significant species
#ggrepel::geom_text_repel(data = comm_env_scores, aes(x=3*NMDS1, y=3*NMDS2, label = env.variables), cex = 15, direction = "both", segment.size = 1) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap


B <- ggplot()+ ggtitle("B)")+ #sets up the plot
  geom_point(data=trait_scrs, aes(x=NMDS1, y=NMDS2, colour=year, shape=year), size = 10)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  #coord_fixed()+ 
  #labs(colour = "treatment", shape = "year")+ # add legend labels for Management and Landuse
  theme_bw(50)+
  theme(legend.position = "none") +
  scale_shape_manual(breaks = c("2017","2019","2021","2100"), values=c("circle", "circle","circle","diamond")) +
  scale_color_manual(breaks = c("2017","2019","2021","2100"), values=c("darkorange", "yellow3","springgreen2","forestgreen")) +
  geom_segment(data = trait_env_scores, aes(x = 0, xend=-3*NMDS1, y=0, yend=3*NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", size=5, linewidth=2) #+ #add vector arrows of significant species #make sure NMDS1 axis is flipped like in 'trait_scrs_plot'
#ggrepel::geom_text_repel(data = trait_env_scores, aes(x=-3*NMDS1, y=3*NMDS2, label = env.variables), cex = 15, direction = "both", segment.size = 1) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap  #make sure NMDS1 axis is flipped like in 'trait_scrs_plot'

#jaccard
C <- ggplot()+ ggtitle("C)")+#sets up the plot
  geom_point(data=jcomm_scrs, aes(x=NMDS1, y=NMDS2, colour=year, shape=year), size = 10)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  #coord_fixed()+
  #labs(colour = "treatment", shape = "year")+ # add legend labels for Management and Landuse
  theme_bw(50)+
  theme(legend.position = "none") +
  scale_shape_manual(breaks = c("2017","2019","2021","2100"), values=c("circle", "circle","circle","diamond")) +
  scale_color_manual(breaks = c("2017","2019","2021","2100"), values=c("darkorange", "yellow3","springgreen2","forestgreen")) +
  geom_segment(data = jcomm_env_scores, aes(x = 0, xend=3*NMDS1, y=0, yend=-3*NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", size=5, linewidth=2) #+ #add vector arrows of significant species
#ggrepel::geom_text_repel(data = comm_env_scores, aes(x=3*NMDS1, y=3*NMDS2, label = env.variables), cex = 15, direction = "both", segment.size = 1) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

#jaccard
D <- ggplot()+ ggtitle("D)")+ #sets up the plot
  geom_point(data=jtrait_scrs, aes(x=NMDS1, y=NMDS2, colour=year, shape=year), size = 10)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  #coord_fixed()+ 
  #labs(colour = "treatment", shape = "year")+ # add legend labels for Management and Landuse
  theme_bw(50)+
  theme(legend.position = "none") +
  scale_shape_manual(breaks = c("2017","2019","2021","2100"), values=c("circle", "circle","circle","diamond")) +
  scale_color_manual(breaks = c("2017","2019","2021","2100"), values=c("darkorange", "yellow3","springgreen2","forestgreen")) +
  geom_segment(data = jtrait_env_scores, aes(x = 0, xend=-3*NMDS1, y=0, yend=3*NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", size=5, linewidth=2) #+ #add vector arrows of significant species #make sure NMDS1 axis is flipped like in 'trait_scrs_plot'


library("gridExtra")
library("grid")
library("cowplot")

plots <- list(A, B, C, D)

plot_layout <- rbind(c(1,2,
                       3,4))

#grid_layout <- grid.arrange(grobs = plots, layout_matrix = plot_layout) # does not align panels

fig4 <- plot_grid(A, B, C, D, align = "v", axis="l", nrow = 2, rel_heights = c(1/2, 1/2), rel_widths = c(1/2, 1/2))

ggplot2::ggsave("./Fig4.png", fig4, height = 30, width = 30)


#CHECK FIG. 4 labels:

#view panel A labels
ggplot()+ ggtitle("A)")+#sets up the plot
  geom_point(data=comm_scrs, aes(x=NMDS1, y=NMDS2, colour=year, shape=year), size = 2)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  #coord_fixed()+
  #labs(colour = "treatment", shape = "year")+ # add legend labels for Management and Landuse
  #theme_bw(50)+
  theme(legend.position = "none") +
  scale_shape_manual(breaks = c("2017","2019","2021","2100"), values=c("circle", "circle","circle","diamond")) +
  scale_color_manual(breaks = c("2017","2019","2021","2100"), values=c("darkorange", "yellow3","springgreen2","forestgreen")) +
  geom_segment(data = comm_env_scores, aes(x = 0, xend=3*NMDS1, y=0, yend=3*NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", size=1, linewidth=1) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = comm_env_scores, aes(x=3*NMDS1, y=3*NMDS2, label = env.variables), cex = 9, direction = "both", segment.size = 1) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

#view panel B labels
  ggplot()+ ggtitle("B)")+ #sets up the plot
  geom_point(data=trait_scrs, aes(x=NMDS1, y=NMDS2, colour=year, shape=year), size = 2)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  #coord_fixed()+ 
  #labs(colour = "treatment", shape = "year")+ # add legend labels for Management and Landuse
  #theme_bw(50)+
  theme(legend.position = "none") +
  scale_shape_manual(breaks = c("2017","2019","2021","2100"), values=c("circle", "circle","circle","diamond")) +
  scale_color_manual(breaks = c("2017","2019","2021","2100"), values=c("darkorange", "yellow3","springgreen2","forestgreen")) +
  geom_segment(data = trait_env_scores, aes(x = 0, xend=-3*NMDS1, y=0, yend=3*NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", size=1, linewidth=1) + #add vector arrows of significant species #make sure NMDS1 axis is flipped like in 'trait_scrs_plot'
  ggrepel::geom_text_repel(data = trait_env_scores, aes(x=-3*NMDS1, y=3*NMDS2, label = env.variables), cex = 9, direction = "both", segment.size = 1) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap  #make sure NMDS1 axis is flipped like in 'trait_scrs_plot'

 ggplot()+ ggtitle("C)")+#sets up the plot
    geom_point(data=jcomm_scrs, aes(x=NMDS1, y=NMDS2, colour=year, shape=year), size = 2)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
    #coord_fixed()+
    #labs(colour = "treatment", shape = "year")+ # add legend labels for Management and Landuse
    #theme_bw(50)+
    theme(legend.position = "none") +
    scale_shape_manual(breaks = c("2017","2019","2021","2100"), values=c("circle", "circle","circle","diamond")) +
    scale_color_manual(breaks = c("2017","2019","2021","2100"), values=c("darkorange", "yellow3","springgreen2","forestgreen")) +
    geom_segment(data = jcomm_env_scores, aes(x = 0, xend=3*NMDS1, y=0, yend=-3*NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", size=1, linewidth=1) + #add vector arrows of significant species
    ggrepel::geom_text_repel(data = jcomm_env_scores, aes(x=3*NMDS1, y=-3*NMDS2, label = env.variables), cex = 9, direction = "both", segment.size = 1) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  
  #jaccard
  ggplot()+ ggtitle("D)")+ #sets up the plot
    geom_point(data=jtrait_scrs, aes(x=NMDS1, y=NMDS2, colour=year, shape=year), size = 2)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
    #coord_fixed()+ 
    #labs(colour = "treatment", shape = "year")+ # add legend labels for Management and Landuse
    #theme_bw(50)+
    theme(legend.position = "none") +
    scale_shape_manual(breaks = c("2017","2019","2021","2100"), values=c("circle", "circle","circle","diamond")) +
    scale_color_manual(breaks = c("2017","2019","2021","2100"), values=c("darkorange", "yellow3","springgreen2","forestgreen")) +
    geom_segment(data = jtrait_env_scores, aes(x = 0, xend=-3*NMDS1, y=0, yend=3*NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", size=1, linewidth=1) + #add vector arrows of significant species #make sure NMDS1 axis is flipped like in 'trait_scrs_plot'
    ggrepel::geom_text_repel(data = jtrait_env_scores, aes(x=-3*NMDS1, y=3*NMDS2, label = env.variables), cex = 9, direction = "both", segment.size = 1) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  
  #--------------------------------------------------------------------------------------------------------------------------------------
  
  
  #-----------------------------------------------------------------------------------------------------------------------
  
  ############ PLOTTING NMDS SI FIGURE
  
  #making function that smooths geom_path
  smooth_it <- function(x, y, n = 1000, method = "natural") 
  {
    t <- seq_along(x)
    new_t <- seq(min(t), max(t), length.out = n)
    new_x <- spline(t, x, xout = new_t, method = method)$y
    new_y <- spline(t, y, xout = new_t, method = method)$y
    data.frame(t = new_t, x = new_x, y = new_y)
  }
  
  #----------------------------------------
  
  #Species Data
  
  #generates centroids by year for species data
  comm_centroid <- comm_scrs %>% group_by(year) %>% dplyr::summarise( #calculate group centroid
    mean_NMDS1 = mean(NMDS1),    
    mean_NMDS2 = mean(NMDS2),)
  
  #smooths species data centroids
  comm_smooth <- smooth_it(comm_centroid$mean_NMDS1, comm_centroid$mean_NMDS2)
  
  #Functional Trait Data
  #generates centroids by year for trait data
  trait_centroid <- trait_scrs %>% group_by(year) %>% dplyr::summarise( #calculate group centroid
    mean_NMDS1 = mean(NMDS1),    
    mean_NMDS2 = mean(NMDS2),)
  
  #smooths trait data centroids
  trait_smooth <- smooth_it(trait_centroid$mean_NMDS1, trait_centroid$mean_NMDS2) 
  
  ###JACCARD
  #generates centroids by year for species data
  jcomm_centroid <- jcomm_scrs %>% group_by(year) %>% dplyr::summarise( #calculate group centroid
    mean_NMDS1 = mean(NMDS1),    
    mean_NMDS2 = mean(NMDS2),)
  
  #smooths species data centroids
  jcomm_smooth <- smooth_it(jcomm_centroid$mean_NMDS1, jcomm_centroid$mean_NMDS2)
  
  #Functional Trait Data
  #generates centroids by year for trait data
  jtrait_centroid <- jtrait_scrs %>% group_by(year) %>% dplyr::summarise( #calculate group centroid
    mean_NMDS1 = mean(NMDS1),    
    mean_NMDS2 = mean(NMDS2),)
  
  #smooths trait data centroids
  jtrait_smooth <- smooth_it(jtrait_centroid$mean_NMDS1, jtrait_centroid$mean_NMDS2) 
  
  
  #SI I: Figure S2
  
  
  library("ggrepel")
  
  #making full panel figure
  
  A1 <- ggplot()+ #sets up the plot
    geom_point(data=comm_scrs, aes(x=NMDS1, y=NMDS2, colour=year, shape=year), size = 10)+ggtitle("C)")+ #adds site points to plot, shape determined by Landuse, colour determined by Management
    #coord_fixed()+
    theme_bw(50)+  #theme(panel.background = linewidth(fill = NA, colour = "black", size = 1, linetype = "solid"))+
    theme(legend.position = "none") +
    #labs(colour = "year", shape = "treatment")+ # add legend labels for Management and Landuse
    #theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) +# add legend at right of plot
    geom_point(data=comm_centroid, aes(x=mean_NMDS1, y=mean_NMDS2, fill=year), colour="black", pch=23, size=30, stroke=2) +
    geom_path(data = comm_smooth, aes(x=x, y=y), linetype=1, size=3, arrow = arrow(angle=15, type="open", length = unit(0.75, "inches"))) +
    #geom_path(data = comm_centroid[!(comm_centroid$year %in% c("2017", "2019")),], aes(x=mean_NMDS1, y=mean_NMDS2), size=3, colour="grey64", lty=6, arrow = arrow(angle=15, type="open", length = unit(1, "inches"))) +
    stat_ellipse(data=comm_scrs, aes(x=NMDS1, y=NMDS2, colour = comm_scrs$year, group=comm_scrs$year), linewidth=3, type = "norm") +
    scale_shape_manual(breaks = c("2017","2019","2021","2100"), values=c("circle", "circle","circle","diamond")) +
    scale_color_manual(breaks = c("2017","2019","2021","2100"), values=c("darkorange", "yellow3","springgreen2","forestgreen")) +
    scale_fill_manual(breaks = c("2017","2019","2021","2100"), values=c("darkorange", "yellow3","springgreen2","forestgreen"))
  
  B1 <- ggplot()+ #sets up the plot
    geom_point(data=trait_scrs, aes(x=NMDS1, y=NMDS2, colour=year, shape=year), size = 10)+ggtitle("D)")+ #adds site points to plot, shape determined by Landuse, colour determined by Management
    #coord_fixed()+
    theme_bw(50)+  #theme(panel.background = linewidth(fill = NA, colour = "black", size = 1, linetype = "solid"))+
    theme(legend.position = "none") +
    #labs(colour = "year", shape = "treatment")+ # add legend labels for Management and Landuse
    #theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) +# add legend at right of plot
    geom_point(data=trait_centroid, aes(x=mean_NMDS1, y=mean_NMDS2, fill=year), colour="black", pch=23, size=30, stroke=2) +
    geom_path(data = trait_smooth, aes(x=x, y=y), linetype=1, size=3, arrow = arrow(angle=15, type="open", length = unit(0.75, "inches"))) +
    #geom_path(data = comm_centroid[!(comm_centroid$year %in% c("2017", "2019")),], aes(x=mean_NMDS1, y=mean_NMDS2), size=3, colour="grey64", lty=6, arrow = arrow(angle=15, type="open", length = unit(1, "inches"))) +
    stat_ellipse(data=trait_scrs, aes(x=NMDS1, y=NMDS2, colour = trait_scrs$year, group=trait_scrs$year), linewidth=3, type = "norm") +
    scale_shape_manual(breaks = c("2017","2019","2021","2100"), values=c("circle", "circle","circle","diamond")) +
    scale_color_manual(breaks = c("2017","2019","2021","2100"), values=c("darkorange", "yellow3","springgreen2","forestgreen")) +
    scale_fill_manual(breaks = c("2017","2019","2021","2100"), values=c("darkorange", "yellow3","springgreen2","forestgreen"))
  
  C1 <- ggplot()+ #sets up the plot
    geom_point(data=jcomm_scrs, aes(x=NMDS1, y=NMDS2, colour=year, shape=year), size = 10)+ggtitle("C)")+ #adds site points to plot, shape determined by Landuse, colour determined by Management
    #coord_fixed()+
    theme_bw(50)+  #theme(panel.background = linewidth(fill = NA, colour = "black", size = 1, linetype = "solid"))+
    theme(legend.position = "none") +
    #labs(colour = "year", shape = "treatment")+ # add legend labels for Management and Landuse
    #theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) +# add legend at right of plot
    geom_point(data=jcomm_centroid, aes(x=mean_NMDS1, y=mean_NMDS2, fill=year), colour="black", pch=23, size=30, stroke=2) +
    geom_path(data = jcomm_smooth, aes(x=x, y=y), linetype=1, size=3, arrow = arrow(angle=15, type="open", length = unit(0.75, "inches"))) +
    #geom_path(data = comm_centroid[!(comm_centroid$year %in% c("2017", "2019")),], aes(x=mean_NMDS1, y=mean_NMDS2), size=3, colour="grey64", lty=6, arrow = arrow(angle=15, type="open", length = unit(1, "inches"))) +
    stat_ellipse(data=jcomm_scrs, aes(x=NMDS1, y=NMDS2, colour = jcomm_scrs$year, group=jcomm_scrs$year), linewidth=3, type = "norm") +
    scale_shape_manual(breaks = c("2017","2019","2021","2100"), values=c("circle", "circle","circle","diamond")) +
    scale_color_manual(breaks = c("2017","2019","2021","2100"), values=c("darkorange", "yellow3","springgreen2","forestgreen")) +
    scale_fill_manual(breaks = c("2017","2019","2021","2100"), values=c("darkorange", "yellow3","springgreen2","forestgreen"))
  
  
  D1 <- ggplot()+ #sets up the plot
    geom_point(data=jtrait_scrs, aes(x=NMDS1, y=NMDS2, colour=year, shape=year), size = 10)+ggtitle("D)")+ #adds site points to plot, shape determined by Landuse, colour determined by Management
    #coord_fixed()+
    theme_bw(50)+  #theme(panel.background = linewidth(fill = NA, colour = "black", size = 1, linetype = "solid"))+
    theme(legend.position = "none") +
    #labs(colour = "year", shape = "treatment")+ # add legend labels for Management and Landuse
    #theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) +# add legend at right of plot
    geom_point(data=jtrait_centroid, aes(x=mean_NMDS1, y=mean_NMDS2, fill=year), colour="black", pch=23, size=30, stroke=2) +
    geom_path(data = jtrait_smooth, aes(x=x, y=y), linetype=1, size=3, arrow = arrow(angle=15, type="open", length = unit(0.75, "inches"))) +
    #geom_path(data = comm_centroid[!(comm_centroid$year %in% c("2017", "2019")),], aes(x=mean_NMDS1, y=mean_NMDS2), size=3, colour="grey64", lty=6, arrow = arrow(angle=15, type="open", length = unit(1, "inches"))) +
    stat_ellipse(data=jtrait_scrs, aes(x=NMDS1, y=NMDS2, colour = jtrait_scrs$year, group=jtrait_scrs$year), linewidth=3, type = "norm") +
    scale_shape_manual(breaks = c("2017","2019","2021","2100"), values=c("circle", "circle","circle","diamond")) +
    scale_color_manual(breaks = c("2017","2019","2021","2100"), values=c("darkorange", "yellow3","springgreen2","forestgreen")) +
    scale_fill_manual(breaks = c("2017","2019","2021","2100"), values=c("darkorange", "yellow3","springgreen2","forestgreen"))
  
  
  library("gridExtra")
  library("grid")
  library("cowplot")
  
  plots <- list(A1, B1, C1, D1)
  
  plot_layout <- rbind(c(1,2,
                         3,4))
  

  grid_layout <- plot_grid(A1, B1, C1, D1, align = "v", axis="l", nrow = 2, rel_heights = c(1/2, 1/2), rel_widths = c(1/2, 1/2))
  
  ggsave("./SI1-FigS2.png", grid_layout, height = 30, width = 30)
  