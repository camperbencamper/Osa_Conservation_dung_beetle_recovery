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
library("tidyr")

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
trait_df$site_id <- as.factor(trait_df$site_id)
#---------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------------------


#glmmTMB

tmbr1 <- glmmTMB(richness ~ 1,
                 family = genpois(link="log"),
                 data = trait_df[trait_df$treatment != "P",],
                 control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))) #optimizer argument to avoid warning

tmbr1.5 <- glmmTMB(richness ~ 1 + (1 |site_id),
                   family = genpois(link="log"),
                   data = trait_df[trait_df$treatment != "P",],
                   control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))) #optimizer argument to avoid warning)

tmbr2 <- glmmTMB(richness ~ treatment,
                 family = genpois(link="log"),
                 data = trait_df[trait_df$treatment != "P",],
                 control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))) #optimizer argument to avoid warning)

tmbr2.5 <- glmmTMB(richness ~ treatment + (1 |site_id),
                   family = genpois(link="log"),
                   data = trait_df[trait_df$treatment != "P",],
                   control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))) #optimizer argument to avoid warning)

tmbr3 <- glmmTMB(richness ~ year,
                 family = genpois(link="log"),
                 data = trait_df[trait_df$treatment != "P",],
                 control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))) #optimizer argument to avoid warning)
#Table 1
emmeans::emmeans(tmbr3, pairwise ~ year, type = "response", data=trait_df[trait_df$treatment != "P",]) #post-hoc
#SI I: Table S6
summary(tmbr3)

tmbr3.5 <- glmmTMB(richness ~ year + (1 |site_id),
                   family = genpois(link="log"),
                   data = trait_df[trait_df$treatment != "P",],
                   control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))) #optimizer argument to avoid warning)
emmeans::emmeans(tmbr3.5, pairwise ~ year, type = "response", data=trait_df[trait_df$treatment != "P",]) #post-hoc
#SI I: Table S6
summary(tmbr3.5)

tmbr4 <- glmmTMB(richness ~ treatment + year,
                 family = genpois(link="log"),
                 data = trait_df[trait_df$treatment != "P",],
                 control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))) #optimizer argument to avoid warning)

tmbr4.5 <- glmmTMB(richness ~ treatment + year + (1 |site_id),
                   family = genpois(link="log"),
                   data = trait_df[trait_df$treatment != "P",],
                   control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))) #optimizer argument to avoid warning)

tmbr5 <- glmmTMB(richness ~ treatment*year,
                 family = genpois(link="log"),
                 data = trait_df[trait_df$treatment != "P",],
                 control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))) #optimizer argument to avoid warning)

tmbr5.5 <- glmmTMB(richness ~ treatment*year + (1 |site_id),
                   family = genpois(link="log"),
                   data = trait_df[trait_df$treatment != "P",],
                   control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))) #optimizer argument to avoid warning)


AIC(tmbr1, tmbr1.5, tmbr2,tmbr2.5, tmbr3, tmbr3.5, tmbr4, tmbr4.5, tmbr5, tmbr5.5)

#-----------------------------------------------------------------------------------------------------------------------------------------------

