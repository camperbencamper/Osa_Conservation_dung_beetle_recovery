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

trait_map <- read.table("mafun.txt", sep="\t", header=TRUE)
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

### reading in data
size_matrix <- read.csv("DB_size_matrix_9-28-22.csv", header=TRUE)
row.names(size_matrix) <- size_matrix[,1]
size_matrix[,1] <- NULL
size_matrix$year <- factor(size_matrix$year, levels=c("2017","2019","2021"))
size_matrix$treatment <- factor(size_matrix$treatment, levels=c("N", "L", "M", "H", "P"))
size_matrix$site_id <- as.factor(size_matrix$site_id)
size_matrix$treatment_year <- paste(size_matrix$treatment, size_matrix$year, sep="_")


diet_matrix <- read.csv("DB_diet_matrix_9-28-22.csv", header=TRUE)
row.names(diet_matrix) <- diet_matrix[,1]
diet_matrix[,1] <- NULL
diet_matrix$year <- factor(diet_matrix$year, levels=c("2017","2019","2021"))
diet_matrix$treatment <- factor(diet_matrix$treatment, levels=c("N", "L", "M", "H", "P"))
diet_matrix$site_id <- as.factor(diet_matrix$site_id)
diet_matrix$treatment_year <- paste(diet_matrix$treatment, diet_matrix$year, sep="_")


burial_matrix <- read.csv("DB_burial_matrix_9-28-22.csv", header=TRUE)
row.names(burial_matrix) <- burial_matrix[,1]
burial_matrix[,1] <- NULL
burial_matrix$year <- factor(burial_matrix$year, levels=c("2017","2019","2021"))
burial_matrix$treatment <- factor(burial_matrix$treatment, levels=c("N", "L", "M", "H", "P"))
burial_matrix$site_id <- as.factor(burial_matrix$site_id)
burial_matrix$treatment_year <- paste(burial_matrix$treatment, burial_matrix$year, sep="_")

f_traits <- c("tunneler", "roller", "unknown", "dweller", "small", "medium", "large", "generalist", "coprophagous", "necrophagous") 
metadata <- c("treatment", "year", "treatment_year", "site_id", "site_year")

#------------------------------------------------------------------------------------------------------

#PLOTS 

#------------------------------------------------------------------------------------------------------------

#MAKING SI FIGURE

#SI I: Figure S1

a <- ggplot() + theme_bw(45) + ylab("Small Size Abundance") + xlab("Year") + theme(legend.position = "none") + ggtitle("A)") +
  geom_boxplot(data=size_matrix[size_matrix$treatment!="P",], size=1.75, outlier.size=4, aes(x=year, y=small, colour=treatment)) +
  stat_summary(data=size_matrix[size_matrix$treatment=="P",], mapping = aes(x = year, y = small, colour=treatment), size=4, pch=18, fun = median) +
  stat_summary(data=size_matrix[size_matrix$treatment!="P",], fun = median, geom = 'line', linewidth=1.75, aes(x = year, y = small, colour = treatment, group=treatment), position = position_dodge(width = 0.9)) +
  #geom_boxplot(data=size_matrix[size_matrix$treatment=="P",], size=1.75, outlier.size=4, aes(x=year, y=small, colour=treatment)) +
  scale_color_manual(breaks = c("N","L","M","H","P"), values=c("red2", "dodgerblue","chartreuse3","darkorchid3","darkred"))

b <- ggplot() + theme_bw(45) + ylab("Medium Size Abundance") + xlab("Year") + theme(legend.position = "none") + ggtitle("B)") +
  geom_boxplot(data=size_matrix[size_matrix$treatment!="P",], size=1.75, outlier.size=4, aes(x=year, y=medium, colour=treatment)) +
  stat_summary(data=size_matrix[size_matrix$treatment=="P",], mapping = aes(x = year, y = medium, colour=treatment), size=4, pch=18, fun = median) +
  stat_summary(data=size_matrix[size_matrix$treatment!="P",], fun = median, geom = 'line', linewidth=1.75, aes(x = year, y = medium, colour = treatment, group=treatment), position = position_dodge(width = 0.9)) +
  #geom_boxplot(data=size_matrix[size_matrix$treatment=="P",], size=1.75, outlier.size=4, aes(x=year, y=medium, colour=treatment)) +
  scale_color_manual(breaks = c("N","L","M","H","P"), values=c("red2", "dodgerblue","chartreuse3","darkorchid3","darkred"))

c <- ggplot() + theme_bw(45) + ylab("Large Size Abundance") + xlab("Year") + theme(legend.position = "none") + ggtitle("C)") +
  geom_boxplot(data=size_matrix[size_matrix$treatment!="P",], size=1.75, outlier.size=4, aes(x=year, y=large, colour=treatment)) +
  stat_summary(data=size_matrix[size_matrix$treatment=="P",], mapping = aes(x = year, y = large, colour=treatment), size=4, pch=18, fun = median) +
  stat_summary(data=size_matrix[size_matrix$treatment!="P",], fun = median, geom = 'line', linewidth=1.75, aes(x = year, y = large, colour = treatment, group=treatment), position = position_dodge(width = 0.9)) +
  #geom_boxplot(data=size_matrix[size_matrix$treatment=="P",], size=1.75, outlier.size=4, aes(x=year, y=large, colour=treatment)) +
  scale_color_manual(breaks = c("N","L","M","H","P"), values=c("red2", "dodgerblue","chartreuse3","darkorchid3","darkred"))

d <- ggplot() + theme_bw(45) + ylab("Generalist Diet Abundance") + xlab("Year") + theme(legend.position = "none") + ggtitle("D)") +
  geom_boxplot(data=diet_matrix[diet_matrix$treatment!="P",], size=1.75, outlier.size=4, aes(x=year, y=generalist, colour=treatment)) +
  stat_summary(data=diet_matrix[diet_matrix$treatment=="P",], mapping = aes(x = year, y = generalist, colour=treatment), size=4, pch=18, fun = median) +
  stat_summary(data=diet_matrix[diet_matrix$treatment!="P",], fun = median, geom = 'line', linewidth=1.75, aes(x = year, y = generalist, colour = treatment, group=treatment), position = position_dodge(width = 0.9)) +
  #geom_boxplot(data=diet_matrix[diet_matrix$treatment=="P",], size=1.75, outlier.size=4, aes(x=year, y=generalist, colour=treatment)) +
  scale_color_manual(breaks = c("N","L","M","H","P"), values=c("red2", "dodgerblue","chartreuse3","darkorchid3","darkred"))

e <- ggplot() + theme_bw(45) + ylab("Coprophagous Diet Abundance") + xlab("Year") + theme(legend.position = "none") + ggtitle("E)") +
  geom_boxplot(data=diet_matrix[diet_matrix$treatment!="P",], size=1.75, outlier.size=4, aes(x=year, y=coprophagous, colour=treatment)) +
  stat_summary(data=diet_matrix[diet_matrix$treatment=="P",], mapping = aes(x = year, y = coprophagous, colour=treatment), size=4, pch=18, fun = median) +
  stat_summary(data=diet_matrix[diet_matrix$treatment!="P",], fun = median, geom = 'line', linewidth=1.75, aes(x = year, y = coprophagous, colour = treatment, group=treatment), position = position_dodge(width = 0.9)) +
  #geom_boxplot(data=diet_matrix[diet_matrix$treatment=="P",], size=1.75, outlier.size=4, aes(x=year, y=coprophagous, colour=treatment)) +
  scale_color_manual(breaks = c("N","L","M","H","P"), values=c("red2", "dodgerblue","chartreuse3","darkorchid3","darkred"))

f <- ggplot() + theme_bw(45) + ylab("Necrophagous Diet Abundance") + xlab("Year") + theme(legend.position = "none") + ggtitle("F)") +
  geom_boxplot(data=diet_matrix[diet_matrix$treatment!="P",], size=1.75, outlier.size=4, aes(x=year, y=necrophagous, colour=treatment)) +
  stat_summary(data=diet_matrix[diet_matrix$treatment=="P",], mapping = aes(x = year, y = necrophagous, colour=treatment), size=4, pch=18, fun = median) +
  stat_summary(data=diet_matrix[diet_matrix$treatment!="P",], fun = median, geom = 'line', linewidth=1.75, aes(x = year, y = necrophagous, colour = treatment, group=treatment), position = position_dodge(width = 0.9)) +
  #geom_boxplot(data=diet_matrix[diet_matrix$treatment=="P",], size=1.75, outlier.size=4, aes(x=year, y=necrophagous, colour=treatment)) +
  scale_color_manual(breaks = c("N","L","M","H","P"), values=c("red2", "dodgerblue","chartreuse3","darkorchid3","darkred"))

g <- ggplot() + theme_bw(45) + ylab("Tunneler Abundance") + xlab("Year") + theme(legend.position = "none") + ggtitle("G)") +
  geom_boxplot(data=burial_matrix[burial_matrix$treatment!="P",], size=1.75, outlier.size=4, aes(x=year, y=tunneler, colour=treatment)) +
  stat_summary(data=burial_matrix[burial_matrix$treatment=="P",], mapping = aes(x = year, y = tunneler, colour=treatment), size=4, pch=18, fun = median) +
  stat_summary(data=burial_matrix[burial_matrix$treatment!="P",], fun = median, geom = 'line', linewidth=1.75, aes(x = year, y = tunneler, colour = treatment, group=treatment), position = position_dodge(width = 0.9)) +
  #geom_boxplot(data=burial_matrix[burial_matrix$treatment=="P",], size=1.75, outlier.size=4, aes(x=year, y=tunneler, colour=treatment)) +
  scale_color_manual(breaks = c("N","L","M","H","P"), values=c("red2", "dodgerblue","chartreuse3","darkorchid3","darkred"))

h <- ggplot() + theme_bw(45) + ylab("Roller Abundance") + xlab("Year") + theme(legend.position = "none") + ggtitle("H)") +
  geom_boxplot(data=burial_matrix[burial_matrix$treatment!="P",], size=1.75, outlier.size=4, aes(x=year, y=roller, colour=treatment)) +
  stat_summary(data=burial_matrix[burial_matrix$treatment=="P",], mapping = aes(x = year, y = roller, colour=treatment), size=4, pch=18, fun = median) +
  stat_summary(data=burial_matrix[burial_matrix$treatment!="P",], fun = median, geom = 'line', linewidth=1.75, aes(x = year, y = roller, colour = treatment, group=treatment), position = position_dodge(width = 0.9)) +
  #geom_boxplot(data=burial_matrix[burial_matrix$treatment=="P",], size=1.75, outlier.size=4, aes(x=year, y=roller, colour=treatment)) +
  scale_color_manual(breaks = c("N","L","M","H","P"), values=c("red2", "dodgerblue","chartreuse3","darkorchid3","darkred"))

i <- ggplot() + theme_bw(45) + ylab("Dweller Abundance") + xlab("Year") + theme(legend.position = "none") + ggtitle("I)") +
  geom_boxplot(data=burial_matrix[burial_matrix$treatment!="P",], size=1.75, outlier.size=4, aes(x=year, y=dweller, colour=treatment)) +
  stat_summary(data=burial_matrix[burial_matrix$treatment=="P",], mapping = aes(x = year, y = dweller, colour=treatment), size=4, pch=18, fun = median) +
  stat_summary(data=burial_matrix[burial_matrix$treatment!="P",], fun = median, geom = 'line', linewidth=1.75, aes(x = year, y = dweller, colour = treatment, group=treatment), position = position_dodge(width = 0.9)) +
  #geom_boxplot(data=burial_matrix[burial_matrix$treatment=="P",], size=1.75, outlier.size=4, aes(x=year, y=dweller, colour=treatment)) +
  scale_color_manual(breaks = c("N","L","M","H","P"), values=c("red2", "dodgerblue","chartreuse3","darkorchid3","darkred"))

library("gridExtra")
library("grid")
library("cowplot")

plots <- list(a, b, c, d, e, f, g, h, i)

plot_layout <- rbind(c(1,2,3),
                     c(4,5,6),
                     c(7,8,9))

#grid_layout <- grid.arrange(grobs = plots, layout_matrix = plot_layout) # does not align panels

grid_layout <- plot_grid(a, b, c, d, e, f, g, h, i, align = "v", nrow = 3, rel_heights = c(1/3, 1/3, 1/3), rel_widths = c(1/3, 1/3, 1/3))

ggsave("./SI1_FigS1.png", grid_layout, height = 30, width = 30)


#------------------------------------------------------------------------------------------------------------


#////////////////////

#####################       MODELS

#////////////////////


#-------------------------------------------------------------------------------------------------------------------------


# SIZE

### SMALL ABUNDANCE, no primary forest control

# GLMER

gTA_dung<-glmer(small~treatment + (1 |site_id) , family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gTA_dung)
gYA_dung<-glmer(small~year + (1 |site_id) , family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gYA_dung)
gTYA_dung<-glmer(small~treatment + year + (1 |site_id) , family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gTYA_dung)
gTYIA_dung<-glmer(small~treatment*year + (1 |site_id) , family="poisson", data=size_matrix[size_matrix$treatment != "P",])
#WARNING
dispersion_glmer(gTYIA_dung)


AIC(gTA_dung, gYA_dung, gTYA_dung, gTYIA_dung)
#any of the interaction models best fit

#-------------------------------------------------------------------------------------------------------------------------------

# GLM
aglm1<-glm(small~treatment, family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glm(glm1)
aglm2<-glm(small~year, family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glm(glm2)
aglm3<-glm(small~treatment + year, family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glm(glm3)
aglm4<-glm(small~treatment*year, family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glm(glm4) #error produced from interaction term

AIC(aglm1, aglm2, aglm3, aglm4)
#gYR_dung best fit
plot(glmer2)

#-------------------------------------------------------------------------------------------------------------------------------


# LMER

lTA_dung<-lmer(small~treatment + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(lTA_dung)
lYA_dung<-lmer(small~year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(lYA_dung)
lTYA_dung<-lmer(small~treatment + year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(lTYA_dung)
lTYIA_dung<-lmer(small~treatment*year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(lTYIA_dung) #error produced from interaction term

AIC(lTA_dung, lYA_dung, lTYA_dung, lTYIA_dung)
#any of the interaction models best fit

#---------------------------------------------------------------------------------------------------------------------------------

# GLMER.nb

gnbTA_dung<-glmer.nb(small~treatment + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gnbTA_dung)
gnbYA_dung<-glmer.nb(small~year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gnbYA_dung)
summary(gnbYA_dung)
glmernb_posthoc <- emmeans::emmeans(gnbYA_dung, pairwise ~ year, type = "response") #post-hoc
#pairs(j_posthoc2, adjust = "bonf") #https://stats.stackexchange.com/questions/9425/bonferroni-or-tukey-when-does-the-number-of-comparisons-become-large
plot(glmernb_posthoc, comparisons = TRUE)
histogram(residuals(gnbYA_dung))
plot(residuals(gnbYA_dung))
testZeroInflation(gnbYA_dung)
testDispersion(gnbYA_dung)
testOverdispersion(gnbYA_dung)
plot(simulateResiduals(gnbYA_dung, n=500))
library("statmod")
res <- qresiduals(gnbYA_dung)
plot(year, res)

gnbTYA_dung<-glmer.nb(small~treatment + year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gnbTYA_dung)
gnbTYIA_dung<-glmer.nb(small~treatment*year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gnbTYIA_dung)
#warning

AIC(gnbTA_dung, gnbYA_dung, gnbTYA_dung, gnbTYIA_dung)
#gnbYA_dung
plot(gnbYA_dung)

#---------------------------------------------------------------------------------------------------------------------

AIC(gTYIA_dung, gTTYIA_dung, gYTYIA_dung, gTYTYITA_dung, lTYIA_dung, lTTYIA_dung, lYTYIA_dung, lTYTYITA_dung, gnbYA_dung)
#gnbYA_dung best model for small
plot(gnbYA_dung)

#all small models
AIC(gTA_dung, gYA_dung, gTYA_dung, gTYIA_dung,lTA_dung, lYA_dung, lTYA_dung, lTYIA_dung, gnbTA_dung, gnbYA_dung, gnbTYA_dung, gnbTYIA_dung, aglm1, aglm2, aglm3, aglm4)


#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------


# SIZE

### MEDIUM ABUNDANCE, no primary forest control

# GLMER

gTA_dung<-glmer(medium~treatment + (1 |site_id) , family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gTA_dung)
gYA_dung<-glmer(medium~year + (1 |site_id) , family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gYA_dung)
gTYA_dung<-glmer(medium~treatment + year + (1 |site_id) , family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gTYA_dung)
gTYIA_dung<-glmer(medium~treatment*year + (1 |site_id) , family="poisson", data=size_matrix[size_matrix$treatment != "P",])
#WARNING
dispersion_glmer(gTYIA_dung)


AIC(gTA_dung, gYA_dung, gTYA_dung, gTYIA_dung)
#any of the interaction models best fit

#-------------------------------------------------------------------------------------------------------------------------------

# GLM
aglm1<-glm(medium~treatment, family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glm(glm1)
aglm2<-glm(medium~year, family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glm(glm2)
aglm3<-glm(medium~treatment + year, family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glm(glm3)
aglm4<-glm(medium~treatment*year, family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glm(glm4) #error produced from interaction term

AIC(aglm1, aglm2, aglm3, aglm4)
#gYR_dung best fit
plot(glmer2)

#-------------------------------------------------------------------------------------------------------------------------------


# LMER

lTA_dung<-lmer(medium~treatment + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(lTA_dung)
lYA_dung<-lmer(medium~year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(lYA_dung)
lTYA_dung<-lmer(medium~treatment + year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(lTYA_dung)
lTYIA_dung<-lmer(medium~treatment*year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(lTYIA_dung) #error produced from interaction term

AIC(lTA_dung, lYA_dung, lTYA_dung, lTYIA_dung)
#any of the interaction models best fit

#---------------------------------------------------------------------------------------------------------------------------------

# GLMER.nb

gnbTA_dung<-glmer.nb(medium~treatment + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gnbTA_dung)
gnbYA_dung<-glmer.nb(medium~year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gnbYA_dung)
summary(gnbYA_dung)
glmernb_posthoc <- emmeans::emmeans(gnbYA_dung, pairwise ~ year, type = "response") #post-hoc
#pairs(j_posthoc2, adjust = "bonf") #https://stats.stackexchange.com/questions/9425/bonferroni-or-tukey-when-does-the-number-of-comparisons-become-large
plot(glmernb_posthoc, comparisons = TRUE)
histogram(residuals(gnbYA_dung))
plot(residuals(gnbYA_dung))
testZeroInflation(gnbYA_dung)
testDispersion(gnbYA_dung)
testOverdispersion(gnbYA_dung)
plot(simulateResiduals(gnbYA_dung, n=500))
library("statmod")
res <- qresiduals(gnbYA_dung)
plot(year, res)

gnbTYA_dung<-glmer.nb(medium~treatment + year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gnbTYA_dung)
gnbTYIA_dung<-glmer.nb(medium~treatment*year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gnbTYIA_dung)
#warning

AIC(gnbTA_dung, gnbYA_dung, gnbTYA_dung, gnbTYIA_dung)
#gnbYA_dung
plot(gnbYA_dung)

#---------------------------------------------------------------------------------------------------------------------

AIC(gTYIA_dung, gTTYIA_dung, gYTYIA_dung, gTYTYITA_dung, lTYIA_dung, lTTYIA_dung, lYTYIA_dung, lTYTYITA_dung, gnbYA_dung)
#gnbYA_dung best model for medium
plot(gnbYA_dung)

#all medium models
AIC(gTA_dung, gYA_dung, gTYA_dung, gTYIA_dung,lTA_dung, lYA_dung, lTYA_dung, lTYIA_dung, gnbTA_dung, gnbYA_dung, gnbTYA_dung, gnbTYIA_dung, aglm1, aglm2, aglm3, aglm4)


#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------


# SIZE

### LARGE ABUNDANCE, no primary forest control

# GLMER

gTA_dung<-glmer(large~treatment + (1 |site_id) , family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gTA_dung)
gYA_dung<-glmer(large~year + (1 |site_id) , family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gYA_dung)
gTYA_dung<-glmer(large~treatment + year + (1 |site_id) , family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gTYA_dung)
gTYIA_dung<-glmer(large~treatment*year + (1 |site_id) , family="poisson", data=size_matrix[size_matrix$treatment != "P",])
#WARNING
dispersion_glmer(gTYIA_dung)


AIC(gTA_dung, gYA_dung, gTYA_dung, gTYIA_dung)
#any of the interaction models best fit

#-------------------------------------------------------------------------------------------------------------------------------

# GLM
aglm1<-glm(large~treatment, family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glm(glm1)
aglm2<-glm(large~year, family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glm(glm2)
aglm3<-glm(large~treatment + year, family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glm(glm3)
aglm4<-glm(large~treatment*year, family="poisson", data=size_matrix[size_matrix$treatment != "P",])
dispersion_glm(glm4) #error produced from interaction term

AIC(aglm1, aglm2, aglm3, aglm4)
#gYR_dung best fit
plot(glmer2)

#-------------------------------------------------------------------------------------------------------------------------------


# LMER

lTA_dung<-lmer(large~treatment + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(lTA_dung)
lYA_dung<-lmer(large~year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(lYA_dung)
lTYA_dung<-lmer(large~treatment + year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(lTYA_dung)
lTYIA_dung<-lmer(large~treatment*year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(lTYIA_dung) #error produced from interaction term

AIC(lTA_dung, lYA_dung, lTYA_dung, lTYIA_dung)
#any of the interaction models best fit

#---------------------------------------------------------------------------------------------------------------------------------

# GLMER.nb

gnbTA_dung<-glmer.nb(large~treatment + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gnbTA_dung)
gnbYA_dung<-glmer.nb(large~year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gnbYA_dung)
summary(gnbYA_dung)
glmernb_posthoc <- emmeans::emmeans(gnbYA_dung, pairwise ~ year, type = "response") #post-hoc
#pairs(j_posthoc2, adjust = "bonf") #https://stats.stackexchange.com/questions/9425/bonferroni-or-tukey-when-does-the-number-of-comparisons-become-large
plot(glmernb_posthoc, comparisons = TRUE)
histogram(residuals(gnbYA_dung))
plot(residuals(gnbYA_dung))
testZeroInflation(gnbYA_dung)
testDispersion(gnbYA_dung)
testOverdispersion(gnbYA_dung)
plot(simulateResiduals(gnbYA_dung, n=500))
library("statmod")
res <- qresiduals(gnbYA_dung)
plot(year, res)

gnbTYA_dung<-glmer.nb(large~treatment + year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gnbTYA_dung)
gnbTYIA_dung<-glmer.nb(large~treatment*year + (1 |site_id), data=size_matrix[size_matrix$treatment != "P",])
dispersion_glmer(gnbTYIA_dung)
#warning

AIC(gnbTA_dung, gnbYA_dung, gnbTYA_dung, gnbTYIA_dung)
#gnbYA_dung
plot(gnbYA_dung)

#---------------------------------------------------------------------------------------------------------------------

AIC(gTYIA_dung, gTTYIA_dung, gYTYIA_dung, gTYTYITA_dung, lTYIA_dung, lTTYIA_dung, lYTYIA_dung, lTYTYITA_dung, gnbYA_dung)
#gnbYA_dung best model for large
plot(gnbYA_dung)

#all large models
AIC(gTA_dung, gYA_dung, gTYA_dung, gTYIA_dung,lTA_dung, lYA_dung, lTYA_dung, lTYIA_dung, gnbTA_dung, gnbYA_dung, gnbTYA_dung, gnbTYIA_dung, aglm1, aglm2, aglm3, aglm4)


#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------


# DIET

### GENERALIST ABUNDANCE, no primary forest control

# GLMER

gTA_dung<-glmer(generalist~treatment + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gTA_dung)
gYA_dung<-glmer(generalist~year + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gYA_dung)
gTYA_dung<-glmer(generalist~treatment + year + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gTYA_dung)
gTYIA_dung<-glmer(generalist~treatment*year + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
#WARNING
dispersion_glmer(gTYIA_dung)


AIC(gTA_dung, gYA_dung, gTYA_dung, gTYIA_dung)
#any of the interaction models best fit

#-------------------------------------------------------------------------------------------------------------------------------

# GLM
aglm1<-glm(generalist~treatment, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm1)
aglm2<-glm(generalist~year, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm2)
aglm3<-glm(generalist~treatment + year, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm3)
aglm4<-glm(generalist~treatment*year, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm4) #error produced from interaction term

AIC(aglm1, aglm2, aglm3, aglm4)
#gYR_dung best fit
plot(glmer2)

#-------------------------------------------------------------------------------------------------------------------------------


# LMER

lTA_dung<-lmer(generalist~treatment + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lTA_dung)
lYA_dung<-lmer(generalist~year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lYA_dung)
lTYA_dung<-lmer(generalist~treatment + year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lTYA_dung)
lTYIA_dung<-lmer(generalist~treatment*year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lTYIA_dung) #error produced from interaction term

AIC(lTA_dung, lYA_dung, lTYA_dung, lTYIA_dung)
#any of the interaction models best fit

#---------------------------------------------------------------------------------------------------------------------------------

# GLMER.nb

gnbTA_dung<-glmer.nb(generalist~treatment + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbTA_dung)
gnbYA_dung<-glmer.nb(generalist~year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbYA_dung)
summary(gnbYA_dung)
glmernb_posthoc <- emmeans::emmeans(gnbYA_dung, pairwise ~ year, type = "response") #post-hoc
#pairs(j_posthoc2, adjust = "bonf") #https://stats.stackexchange.com/questions/9425/bonferroni-or-tukey-when-does-the-number-of-comparisons-become-generalist
plot(glmernb_posthoc, comparisons = TRUE)
histogram(residuals(gnbYA_dung))
plot(residuals(gnbYA_dung))
testZeroInflation(gnbYA_dung)
testDispersion(gnbYA_dung)
testOverdispersion(gnbYA_dung)
plot(simulateResiduals(gnbYA_dung, n=500))
library("statmod")
res <- qresiduals(gnbYA_dung)
plot(year, res)

gnbTYA_dung<-glmer.nb(generalist~treatment + year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbTYA_dung)
gnbTYIA_dung<-glmer.nb(generalist~treatment*year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbTYIA_dung)
#warning

AIC(gnbTA_dung, gnbYA_dung, gnbTYA_dung, gnbTYIA_dung)
#gnbYA_dung
plot(gnbYA_dung)

#---------------------------------------------------------------------------------------------------------------------

AIC(gTYIA_dung, gTTYIA_dung, gYTYIA_dung, gTYTYITA_dung, lTYIA_dung, lTTYIA_dung, lYTYIA_dung, lTYTYITA_dung, gnbYA_dung)
#gnbYA_dung best model for generalist
plot(gnbYA_dung)

#all generalist models
AIC(gTA_dung, gYA_dung, gTYA_dung, gTYIA_dung,lTA_dung, lYA_dung, lTYA_dung, lTYIA_dung, gnbTA_dung, gnbYA_dung, gnbTYA_dung, gnbTYIA_dung, aglm1, aglm2, aglm3, aglm4)



#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------


# DIET

### COPROPHAGOUS ABUNDANCE, no primary forest control

# GLMER

gTA_dung<-glmer(coprophagous~treatment + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gTA_dung)
gYA_dung<-glmer(coprophagous~year + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gYA_dung)
gTYA_dung<-glmer(coprophagous~treatment + year + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gTYA_dung)
gTYIA_dung<-glmer(coprophagous~treatment*year + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
#WARNING
dispersion_glmer(gTYIA_dung)


AIC(gTA_dung, gYA_dung, gTYA_dung, gTYIA_dung)
#any of the interaction models best fit

#-------------------------------------------------------------------------------------------------------------------------------

# GLM
aglm1<-glm(coprophagous~treatment, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm1)
aglm2<-glm(coprophagous~year, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm2)
aglm3<-glm(coprophagous~treatment + year, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm3)
aglm4<-glm(coprophagous~treatment*year, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm4) #error produced from interaction term

AIC(aglm1, aglm2, aglm3, aglm4)
#gYR_dung best fit
plot(glmer2)

#-------------------------------------------------------------------------------------------------------------------------------


# LMER

lTA_dung<-lmer(coprophagous~treatment + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lTA_dung)
lYA_dung<-lmer(coprophagous~year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lYA_dung)
lTYA_dung<-lmer(coprophagous~treatment + year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lTYA_dung)
lTYIA_dung<-lmer(coprophagous~treatment*year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lTYIA_dung) #error produced from interaction term

AIC(lTA_dung, lYA_dung, lTYA_dung, lTYIA_dung)
#any of the interaction models best fit

#---------------------------------------------------------------------------------------------------------------------------------

# GLMER.nb

gnbTA_dung<-glmer.nb(coprophagous~treatment + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbTA_dung)
gnbYA_dung<-glmer.nb(coprophagous~year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbYA_dung)
summary(gnbYA_dung)
glmernb_posthoc <- emmeans::emmeans(gnbYA_dung, pairwise ~ year, type = "response") #post-hoc
#pairs(j_posthoc2, adjust = "bonf") #https://stats.stackexchange.com/questions/9425/bonferroni-or-tukey-when-does-the-number-of-comparisons-become-coprophagous
plot(glmernb_posthoc, comparisons = TRUE)
histogram(residuals(gnbYA_dung))
plot(residuals(gnbYA_dung))
testZeroInflation(gnbYA_dung)
testDispersion(gnbYA_dung)
testOverdispersion(gnbYA_dung)
plot(simulateResiduals(gnbYA_dung, n=500))
library("statmod")
res <- qresiduals(gnbYA_dung)
plot(year, res)

gnbTYA_dung<-glmer.nb(coprophagous~treatment + year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbTYA_dung)
gnbTYIA_dung<-glmer.nb(coprophagous~treatment*year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbTYIA_dung)
#warning

AIC(gnbTA_dung, gnbYA_dung, gnbTYA_dung, gnbTYIA_dung)
#gnbYA_dung
plot(gnbYA_dung)

#---------------------------------------------------------------------------------------------------------------------

AIC(gTYIA_dung, gTTYIA_dung, gYTYIA_dung, gTYTYITA_dung, lTYIA_dung, lTTYIA_dung, lYTYIA_dung, lTYTYITA_dung, gnbYA_dung)
#gnbYA_dung best model for coprophagous
plot(gnbYA_dung)

#all coprophagous models
AIC(gTA_dung, gYA_dung, gTYA_dung, gTYIA_dung,lTA_dung, lYA_dung, lTYA_dung, lTYIA_dung, gnbTA_dung, gnbYA_dung, gnbTYA_dung, gnbTYIA_dung, aglm1, aglm2, aglm3, aglm4)


#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------


# DIET

### NECROPHAGOUS ABUNDANCE, no primary forest control

# GLMER

gTA_dung<-glmer(necrophagous~treatment + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gTA_dung)
gYA_dung<-glmer(necrophagous~year + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gYA_dung)
gTYA_dung<-glmer(necrophagous~treatment + year + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gTYA_dung)
gTYIA_dung<-glmer(necrophagous~treatment*year + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
#WARNING
dispersion_glmer(gTYIA_dung)


AIC(gTA_dung, gYA_dung, gTYA_dung, gTYIA_dung)
#any of the interaction models best fit

#-------------------------------------------------------------------------------------------------------------------------------

# GLM
aglm1<-glm(necrophagous~treatment, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm1)
aglm2<-glm(necrophagous~year, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm2)
aglm3<-glm(necrophagous~treatment + year, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm3)
aglm4<-glm(necrophagous~treatment*year, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm4) #error produced from interaction term

AIC(aglm1, aglm2, aglm3, aglm4)
#gYR_dung best fit
plot(glmer2)

#-------------------------------------------------------------------------------------------------------------------------------


# LMER

lTA_dung<-lmer(necrophagous~treatment + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lTA_dung)
lYA_dung<-lmer(necrophagous~year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lYA_dung)
lTYA_dung<-lmer(necrophagous~treatment + year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lTYA_dung)
lTYIA_dung<-lmer(necrophagous~treatment*year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lTYIA_dung) #error produced from interaction term

AIC(lTA_dung, lYA_dung, lTYA_dung, lTYIA_dung)
#any of the interaction models best fit

#---------------------------------------------------------------------------------------------------------------------------------

# GLMER.nb

gnbTA_dung<-glmer.nb(necrophagous~treatment + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbTA_dung)
gnbYA_dung<-glmer.nb(necrophagous~year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbYA_dung)
summary(gnbYA_dung)
glmernb_posthoc <- emmeans::emmeans(gnbYA_dung, pairwise ~ year, type = "response") #post-hoc
#pairs(j_posthoc2, adjust = "bonf") #https://stats.stackexchange.com/questions/9425/bonferroni-or-tukey-when-does-the-number-of-comparisons-become-necrophagous
plot(glmernb_posthoc, comparisons = TRUE)
histogram(residuals(gnbYA_dung))
plot(residuals(gnbYA_dung))
testZeroInflation(gnbYA_dung)
testDispersion(gnbYA_dung)
testOverdispersion(gnbYA_dung)
plot(simulateResiduals(gnbYA_dung, n=500))
library("statmod")
res <- qresiduals(gnbYA_dung)
plot(year, res)

gnbTYA_dung<-glmer.nb(necrophagous~treatment + year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbTYA_dung)
gnbTYIA_dung<-glmer.nb(necrophagous~treatment*year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbTYIA_dung)
#warning

AIC(gnbTA_dung, gnbYA_dung, gnbTYA_dung, gnbTYIA_dung)
#gnbYA_dung
plot(gnbYA_dung)

#---------------------------------------------------------------------------------------------------------------------

AIC(gTYIA_dung, gTTYIA_dung, gYTYIA_dung, gTYTYITA_dung, lTYIA_dung, lTTYIA_dung, lYTYIA_dung, lTYTYITA_dung, gnbYA_dung)
#gnbYA_dung best model for necrophagous
plot(gnbYA_dung)

#all necrophagous models
AIC(gTA_dung, gYA_dung, gTYA_dung, gTYIA_dung,lTA_dung, lYA_dung, lTYA_dung, lTYIA_dung, gnbTA_dung, gnbYA_dung, gnbTYA_dung, gnbTYIA_dung, aglm1, aglm2, aglm3, aglm4)


#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------

# BURIAL

### TUNNELER ABUNDANCE, no primary forest control

# GLMER

gTA_dung<-glmer(tunneler~treatment + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gTA_dung)
gYA_dung<-glmer(tunneler~year + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gYA_dung)
gTYA_dung<-glmer(tunneler~treatment + year + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gTYA_dung)
gTYIA_dung<-glmer(tunneler~treatment*year + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
#WARNING
dispersion_glmer(gTYIA_dung)


AIC(gTA_dung, gYA_dung, gTYA_dung, gTYIA_dung)
#any of the interaction models best fit

#-------------------------------------------------------------------------------------------------------------------------------

# GLM
aglm1<-glm(tunneler~treatment, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm1)
aglm2<-glm(tunneler~year, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm2)
aglm3<-glm(tunneler~treatment + year, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm3)
aglm4<-glm(tunneler~treatment*year, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm4) #error produced from interaction term

AIC(aglm1, aglm2, aglm3, aglm4)
#gYR_dung best fit
plot(glmer2)

#-------------------------------------------------------------------------------------------------------------------------------


# LMER

lTA_dung<-lmer(tunneler~treatment + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lTA_dung)
lYA_dung<-lmer(tunneler~year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lYA_dung)
lTYA_dung<-lmer(tunneler~treatment + year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lTYA_dung)
lTYIA_dung<-lmer(tunneler~treatment*year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lTYIA_dung) #error produced from interaction term

AIC(lTA_dung, lYA_dung, lTYA_dung, lTYIA_dung)
#any of the interaction models best fit

#---------------------------------------------------------------------------------------------------------------------------------

# GLMER.nb

gnbTA_dung<-glmer.nb(tunneler~treatment + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbTA_dung)
gnbYA_dung<-glmer.nb(tunneler~year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbYA_dung)
summary(gnbYA_dung)
glmernb_posthoc <- emmeans::emmeans(gnbYA_dung, pairwise ~ year, type = "response") #post-hoc
#pairs(j_posthoc2, adjust = "bonf") #https://stats.stackexchange.com/questions/9425/bonferroni-or-tukey-when-does-the-number-of-comparisons-become-tunneler
plot(glmernb_posthoc, comparisons = TRUE)
histogram(residuals(gnbYA_dung))
plot(residuals(gnbYA_dung))
testZeroInflation(gnbYA_dung)
testDispersion(gnbYA_dung)
testOverdispersion(gnbYA_dung)
plot(simulateResiduals(gnbYA_dung, n=500))
library("statmod")
res <- qresiduals(gnbYA_dung)
plot(year, res)

gnbTYA_dung<-glmer.nb(tunneler~treatment + year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbTYA_dung)
gnbTYIA_dung<-glmer.nb(tunneler~treatment*year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbTYIA_dung)
#warning

AIC(gnbTA_dung, gnbYA_dung, gnbTYA_dung, gnbTYIA_dung)
#gnbYA_dung
plot(gnbYA_dung)

#---------------------------------------------------------------------------------------------------------------------

AIC(gTYIA_dung, gTTYIA_dung, gYTYIA_dung, gTYTYITA_dung, lTYIA_dung, lTTYIA_dung, lYTYIA_dung, lTYTYITA_dung, gnbYA_dung)
#gnbYA_dung best model for tunneler
plot(gnbYA_dung)

#all tunneler models
AIC(gTA_dung, gYA_dung, gTYA_dung, gTYIA_dung,lTA_dung, lYA_dung, lTYA_dung, lTYIA_dung, gnbTA_dung, gnbYA_dung, gnbTYA_dung, gnbTYIA_dung, aglm1, aglm2, aglm3, aglm4)


#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------


# BURIAL

### ROLLER ABUNDANCE, no primary forest control

# GLMER

gTA_dung<-glmer(roller~treatment + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gTA_dung)
gYA_dung<-glmer(roller~year + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gYA_dung)
gTYA_dung<-glmer(roller~treatment + year + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gTYA_dung)
gTYIA_dung<-glmer(roller~treatment*year + (1 |site_id) , family="poisson", data=comm_data[comm_data$treatment != "P",])
#WARNING
dispersion_glmer(gTYIA_dung)


AIC(gTA_dung, gYA_dung, gTYA_dung, gTYIA_dung)
#any of the interaction models best fit

#-------------------------------------------------------------------------------------------------------------------------------

# GLM
aglm1<-glm(roller~treatment, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm1)
aglm2<-glm(roller~year, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm2)
aglm3<-glm(roller~treatment + year, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm3)
aglm4<-glm(roller~treatment*year, family="poisson", data=comm_data[comm_data$treatment != "P",])
dispersion_glm(glm4) #error produced from interaction term

AIC(aglm1, aglm2, aglm3, aglm4)
#gYR_dung best fit
plot(glmer2)

#-------------------------------------------------------------------------------------------------------------------------------


# LMER

lTA_dung<-lmer(roller~treatment + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lTA_dung)
lYA_dung<-lmer(roller~year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lYA_dung)
lTYA_dung<-lmer(roller~treatment + year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lTYA_dung)
lTYIA_dung<-lmer(roller~treatment*year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(lTYIA_dung) #error produced from interaction term

AIC(lTA_dung, lYA_dung, lTYA_dung, lTYIA_dung)
#any of the interaction models best fit

#---------------------------------------------------------------------------------------------------------------------------------

# GLMER.nb

gnbTA_dung<-glmer.nb(roller~treatment + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbTA_dung)
gnbYA_dung<-glmer.nb(roller~year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbYA_dung)
summary(gnbYA_dung)
glmernb_posthoc <- emmeans::emmeans(gnbYA_dung, pairwise ~ year, type = "response") #post-hoc
#pairs(j_posthoc2, adjust = "bonf") #https://stats.stackexchange.com/questions/9425/bonferroni-or-tukey-when-does-the-number-of-comparisons-become-roller
plot(glmernb_posthoc, comparisons = TRUE)
histogram(residuals(gnbYA_dung))
plot(residuals(gnbYA_dung))
testZeroInflation(gnbYA_dung)
testDispersion(gnbYA_dung)
testOverdispersion(gnbYA_dung)
plot(simulateResiduals(gnbYA_dung, n=500))
library("statmod")
res <- qresiduals(gnbYA_dung)
plot(year, res)

gnbTYA_dung<-glmer.nb(roller~treatment + year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbTYA_dung)
gnbTYIA_dung<-glmer.nb(roller~treatment*year + (1 |site_id), data=comm_data[comm_data$treatment != "P",])
dispersion_glmer(gnbTYIA_dung)
#warning

AIC(gnbTA_dung, gnbYA_dung, gnbTYA_dung, gnbTYIA_dung)
#gnbYA_dung
plot(gnbYA_dung)

#---------------------------------------------------------------------------------------------------------------------

AIC(gTYIA_dung, gTTYIA_dung, gYTYIA_dung, gTYTYITA_dung, lTYIA_dung, lTTYIA_dung, lYTYIA_dung, lTYTYITA_dung, gnbYA_dung)
#gnbYA_dung best model for roller
plot(gnbYA_dung)

#all roller models
AIC(gTA_dung, gYA_dung, gTYA_dung, gTYIA_dung,lTA_dung, lYA_dung, lTYA_dung, lTYIA_dung, gnbTA_dung, gnbYA_dung, gnbTYA_dung, gnbTYIA_dung, aglm1, aglm2, aglm3, aglm4)

