

# MODEL BUILDING

library("ggplot2")
library("dplyr")
library("vegan")
library(nlme)
library(lme4)
library(car)
library(blmeco)
library(emmeans)
library("glmmTMB")
#library(MixMod) #not available for this R version, 'anovaTAB' provides p-values


set.seed(1)

#data from all years
comm_data <- read.delim("DB_all.txt")

colnames(comm_data)
#comm_data <- comm_data[,-1]

#creating this vector to avoid subsetting by column numbers when using species matrix
DB_species <- colnames(comm_data[-(1:5)])

#creating a new comm_data column just by plot treatment
#extracting treatments from plot_id column
comm_data$treatment <- substr(comm_data$site_id, 1, 1)
unique(comm_data$site_id)
unique(comm_data$treatment)
rownames(comm_data) <- paste(comm_data$site_id, comm_data$year, sep = "_")

#creating treatment_year column
comm_data$treatment_year <- as.factor(paste(comm_data$treatment, comm_data$year, sep="_"))

#converting species columsn from integers to numeric #WHAT A PAIN IN THE ASS
comm_data[names(comm_data) %in% as.factor(DB_species)] <- lapply(comm_data[names(comm_data) %in% as.factor(DB_species)], FUN = function(y){as.numeric(y)})

comm_data$year <- factor(comm_data$year, levels = c("2017", "2019", "2021"))
comm_data$treatment <- factor(comm_data$treatment, levels = c("N", "L", "M", "H"))
comm_data$site_id <- as.factor(comm_data$site_id)

#-----------------------------------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------------------------------

### RICHNESS, no primary forest control

#glmmTMB
tmbr1 <- glmmTMB(richness ~ 1,
                 family = genpois(link="log"),
                 data = comm_data[comm_data$treatment!="P",])

tmbr1.5 <- glmmTMB(richness ~ 1 + (1 |site_id),
                   family = genpois(link="log"),
                   data = comm_data[comm_data$treatment!="P",])

tmbr2 <- glmmTMB(richness ~ treatment,
                 family = genpois(link="log"),
                 data = comm_data[comm_data$treatment!="P",])

tmbr2.5 <- glmmTMB(richness ~ treatment + (1 |site_id),
                 family = genpois(link="log"),
                 data = comm_data[comm_data$treatment!="P",])

tmbr3 <- glmmTMB(richness ~ year,
                 family = genpois(link="log"),
                 data = comm_data[comm_data$treatment!="P",])
#Table 1
emmeans::emmeans(tmbr3, pairwise ~ year, type = "response", data=comm_data[comm_data$treatment!="P",]) #post-hoc
#SI I: Table S5
summary(tmbr3)

tmbr3.5 <- glmmTMB(richness ~ year + (1 |site_id),
                   family = genpois(link="log"),
                   data = comm_data[comm_data$treatment!="P",])
emmeans::emmeans(tmbr3.5, pairwise ~ year, type = "response", data=comm_data[comm_data$treatment!="P",]) #post-hoc
#SI I: Table S5
summary(tmbr3.5)

tmbr4 <- glmmTMB(richness ~ treatment + year,
                 family = genpois(link="log"),
                 data = comm_data[comm_data$treatment!="P",])

tmbr4.5 <- glmmTMB(richness ~ treatment + year + (1 |site_id),
                 family = genpois(link="log"),
                 data = comm_data[comm_data$treatment!="P",])

tmbr5 <- glmmTMB(richness ~ treatment*year,
                 family = genpois(link="log"),
                 data = comm_data[comm_data$treatment!="P",])

tmbr5.5 <- glmmTMB(richness ~ treatment*year + (1 |site_id),
                 family = genpois(link="log"),
                 data = comm_data[comm_data$treatment!="P",])


AIC(tmbr1, tmbr1.5, tmbr2,tmbr2.5, tmbr3, tmbr3.5, tmbr4, tmbr4.5, tmbr5, tmbr5.5)

#-------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------

### ABUNDANCE, no primary forest control

#glmmTMB
tmba1 <- glmmTMB(abundance ~ 1,
                 family = genpois(link="log"),
                 data = comm_data[comm_data$treatment!="P",])

tmba1.5 <- glmmTMB(abundance ~ 1 + (1 |site_id),
                   family = genpois(link="log"),
                   data = comm_data[comm_data$treatment!="P",])
#SI I: Table S4
summary(tmba1.5)
#emmeans::emmeans(tmba1.5, pairwise ~ 1, type = "response", data=comm_data[comm_data$treatment!="P",])

tmba2 <- glmmTMB(abundance ~ treatment,
                 family = genpois(link="log"),
                 data = comm_data[comm_data$treatment!="P",])

tmba2.5 <- glmmTMB(abundance ~ treatment + (1 |site_id),
                   family = genpois(link="log"),
                   data = comm_data[comm_data$treatment!="P",])

tmba3 <- glmmTMB(abundance ~ year,
                 family = genpois(link="log"),
                 data = comm_data[comm_data$treatment!="P",])

tmba3.5 <- glmmTMB(abundance ~ year + (1 |site_id),
                   family = genpois(link="log"),
                   data = comm_data[comm_data$treatment!="P",])

tmba4 <- glmmTMB(abundance ~ treatment + year,
                 family = genpois(link="log"),
                 data = comm_data[comm_data$treatment!="P",])

tmba4.5 <- glmmTMB(abundance ~ treatment + year + (1 |site_id),
                   family = genpois(link="log"),
                   data = comm_data[comm_data$treatment!="P",])

tmba5 <- glmmTMB(abundance ~ treatment*year,
                 family = genpois(link="log"),
                 data = comm_data[comm_data$treatment!="P",])

tmba5.5 <- glmmTMB(abundance ~ treatment*year + (1 |site_id),
                   family = genpois(link="log"),
                   data = comm_data[comm_data$treatment!="P",])

AIC(tmba1, tmba1.5, tmba2,tmba2.5, tmba3, tmba3.5, tmba4, tmba4.5, tmba5, tmba5.5)

#-------------------------------------------------------------------------------------------------------------------------------
