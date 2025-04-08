
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

### SAMPLE COVERAGE

library("iNEXT")

comm_data2 <- data.frame(t(comm_data[names(comm_data) %in% DB_species]))
comm_data2 <- comm_data2[colSums(comm_data2)>0]
#rich2 <- colSums(comm_data2)


samp_cov <- estimateD(comm_data2, q=0, datatype="abundance", base="coverage")
samp_cov$relative_width <- (samp_cov$qD.UCL - samp_cov$qD.LCL) / samp_cov$qD
median(samp_cov$relative_width)
histogram(samp_cov$relative_width)

median(samp_cov$SC)
histogram(samp_cov$SC)
sum(samp_cov$SC < 0.78)#only 4 samples below 78%

SClow <- samp_cov[samp_cov$SC < 0.80,]
histogram(SClow$SC)

problem_sites <- which(colSums(comm_data2 > 0) == 1)
problem_sites #only have one species present; N=30
length(problem_sites)


#------------------------------------

#adding adjusted richness data to main dataset

comm_data$rich_mod <- samp_cov$qD[match(rownames(comm_data), samp_cov$Assemblage)]
comm_data$rich_mod[is.na(comm_data$rich_mod)] <- 0
comm_data$SC <- samp_cov$SC[match(rownames(comm_data), samp_cov$Assemblage)]
#comm_data$rich2 <- rowSums(comm_data[names(comm_data) %in% DB_species] > 0) #checking to make sure the 'richness' column is correct
comm_data[names(comm_data) %in% c("richness", "rich_mod", "SC")]

mod <- lm(rich_mod ~ richness, comm_data)
summary(mod)
ggplot(comm_data, aes(richness, rich_mod)) +
  geom_point() + geom_smooth(method="lm", se=FALSE) +
  annotate("text", x=quantile(comm_data$richness, .75), y=quantile(comm_data$rich_mod, .25),
           label=paste0("R²=", round(summary(mod)$r.squared, 2)))

#----

comm_data$treatment <- as.character(comm_data$treatment)
comm_data$treatment[is.na(comm_data$treatment)] <- "P"
comm_data$treatment <- as.factor(comm_data$treatment)
comm_data$treatment <- factor(comm_data$treatment, levels = c("N", "L", "M", "H", "P"))

#RICHNESS

#SI II: Figure S1
p1 <- ggplot() + theme_bw(40) + ylab("Species Richness Adjusted") + xlab("Year") + ggtitle("A)") + theme(legend.position = "none") + #ylim(0,8) +
  geom_boxplot(data=comm_data[comm_data$treatment!="P",], size=1.75, outlier.size=4, aes(x=year, y=rich_mod, colour=treatment)) +
  geom_boxplot(data=comm_data[comm_data$treatment=="P",], size=1.75, outlier.size=4, aes(x=year, y=rich_mod, colour=treatment)) +
  scale_color_manual(breaks = c("N","L","M","H","P"), values=c("red2", "dodgerblue","chartreuse3","darkorchid3","darkred"))

#------------------------------------------------------

### RICHNESS, no primary forest control

#glmmTMB
tmbr1 <- glmmTMB(rich_mod ~ 1,
                 family = gaussian(link="identity"),
                 data = comm_data[comm_data$treatment!="P",])

tmbr1.5 <- glmmTMB(rich_mod ~ 1 + (1 |site_id),
                   family = gaussian(link="identity"),
                   data = comm_data[comm_data$treatment!="P",])

tmbr2 <- glmmTMB(rich_mod ~ treatment,
                 family = gaussian(link="identity"),
                 data = comm_data[comm_data$treatment!="P",])

tmbr2.5 <- glmmTMB(rich_mod ~ treatment + (1 |site_id),
                 family = gaussian(link="identity"),
                 data = comm_data[comm_data$treatment!="P",])

tmbr3 <- glmmTMB(rich_mod ~ year,
                 family = gaussian(link="identity"),
                 data = comm_data[comm_data$treatment!="P",])
emmeans::emmeans(tmbr3, pairwise ~ year, type = "response", data=comm_data[comm_data$treatment!="P",]) #post-hoc
summary(tmbr3)
library("DHARMa")
res_tmbr3 <- simulateResiduals(tmbr3)
plot(res_tmbr3)

AIC(tmbr3, tmbr3_idk)

tmbr3.5 <- glmmTMB(rich_mod ~ year + (1 |site_id),
                   family = gaussian(link="identity"),
                   data = comm_data[comm_data$treatment!="P",])
#SI II: Table S1
emmeans::emmeans(tmbr3.5, pairwise ~ year, type = "response", data=comm_data[comm_data$treatment!="P",]) #post-hoc
#SI II: Table S2
summary(tmbr3.5)

tmbr4 <- glmmTMB(rich_mod ~ treatment + year,
                 family = gaussian(link="identity"),
                 data = comm_data[comm_data$treatment!="P",])

tmbr4.5 <- glmmTMB(rich_mod ~ treatment + year + (1 |site_id),
                 family = gaussian(link="identity"),
                 data = comm_data[comm_data$treatment!="P",])

tmbr5 <- glmmTMB(rich_mod ~ treatment*year,
                 family = gaussian(link="identity"),
                 data = comm_data[comm_data$treatment!="P",])

tmbr5.5 <- glmmTMB(rich_mod ~ treatment*year + (1 |site_id),
                 family = gaussian(link="identity"),
                 data = comm_data[comm_data$treatment!="P",])


AIC(tmbr1, tmbr1.5, tmbr2,tmbr2.5, tmbr3, tmbr3.5, tmbr4, tmbr4.5, tmbr5, tmbr5.5)

#-------------------------------------------------------------------------------------------------------------------------------

