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
trait_map$functional_group <- paste(trait_map$size_group, trait_map$diet_breadth, trait_map$burial_strategy, sep="_")

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


### SAMPLE COVERAGE

library("iNEXT")

trait_data2 <- data.frame(t(trait_df[names(trait_df) %in% functions]))
trait_data2 <- trait_data2[colSums(trait_data2)>0]

samp_cov <- estimateD(trait_data2, q=0, datatype="abundance", base="coverage")
median(samp_cov$SC)


#---

trait_df$rich_mod <- samp_cov$qD[match(rownames(trait_df), samp_cov$Assemblage)]
trait_df$rich_mod[is.na(trait_df$rich_mod)] <- 0
trait_df$SC <- samp_cov$SC[match(rownames(trait_df), samp_cov$Assemblage)]
#comm_data$rich2 <- rowSums(trait_df[names(trait_df) %in% DB_species] > 0) #checking to make sure the 'richness' column is correct
trait_df[names(trait_df) %in% c("richness", "rich_mod", "SC")]

mod <- lm(rich_mod ~ richness, trait_df)
summary(mod)
ggplot(trait_df, aes(richness, rich_mod)) +
  geom_point() + geom_smooth(method="lm", se=FALSE) +
  annotate("text", x=quantile(trait_df$richness, .75), y=quantile(comm_data$rich_mod, .25),
           label=paste0("R²=", round(summary(mod)$r.squared, 2)))



library("glmmTMB")

#glmmTMB

tmbr1 <- glmmTMB(rich_mod ~ 1,
                 family = gaussian(link="identity"),
                 data = trait_df[trait_df$treatment != "P",])

tmbr1.5 <- glmmTMB(rich_mod ~ 1 + (1 |site_id),
                   family = gaussian(link="identity"),
                   data = trait_df[trait_df$treatment != "P",])

tmbr2 <- glmmTMB(rich_mod ~ treatment,
                 family = gaussian(link="identity"),
                 data = trait_df[trait_df$treatment != "P",])

tmbr2.5 <- glmmTMB(rich_mod ~ treatment + (1 |site_id),
                   family = gaussian(link="identity"),
                   data = trait_df[trait_df$treatment != "P",])

tmbr3 <- glmmTMB(rich_mod ~ year,
                 family = gaussian(link="identity"),
                 data = trait_df[trait_df$treatment != "P",])
emmeans::emmeans(tmbr3, pairwise ~ year, type = "response", data=trait_df[trait_df$treatment != "P",]) #post-hoc
summary(tmbr3)

tmbr3_idk <- glmmTMB(rich_mod ~ year,
                     family = Gamma(link="log"),
                     data = trait_df[trait_df$treatment != "P",])
AIC(tmbr3, tmbr3_idk)

tmbr3.5 <- glmmTMB(rich_mod ~ year + (1 |site_id),
                   family = gaussian(link="identity"),
                   data = trait_df[trait_df$treatment != "P",])
#SI II: Table S1
emmeans::emmeans(tmbr3.5, pairwise ~ year, type = "response", data=trait_df[trait_df$treatment != "P",]) #post-hoc
#SI II: Table S3
summary(tmbr3.5)

tmbr4 <- glmmTMB(rich_mod ~ treatment + year,
                 family = gaussian(link="identity"),
                 data = trait_df[trait_df$treatment != "P",])

tmbr4.5 <- glmmTMB(rich_mod ~ treatment + year + (1 |site_id),
                   family = gaussian(link="identity"),
                   data = trait_df[trait_df$treatment != "P",])

tmbr5 <- glmmTMB(rich_mod ~ treatment*year,
                 family = gaussian(link="identity"),
                 data = trait_df[trait_df$treatment != "P",])

tmbr5.5 <- glmmTMB(rich_mod ~ treatment*year + (1 |site_id),
                   family = gaussian(link="identity"),
                   data = trait_df[trait_df$treatment != "P",])


AIC(tmbr1, tmbr1.5, tmbr2,tmbr2.5, tmbr3, tmbr3.5, tmbr4, tmbr4.5, tmbr5, tmbr5.5)

#-----------------------------------------------------------------------------------------------------------------------------------------------

#SI II: Figure S1

#p1 from "taxonomy script"DB_sampling_coverage_taxonomy_richness.R"

p2 <- ggplot() + theme_bw(40) + ylab("Functional Richness Adjusted") + xlab("Year") + ggtitle("B)") + theme(legend.position = "none") + #ylim(0,8) +
  geom_boxplot(data=trait_df[trait_df$treatment!="P",], size=1.75, outlier.size=4, aes(x=year, y=rich_mod, colour=treatment)) +
  geom_boxplot(data=trait_df[trait_df$treatment=="P",], size=1.75, outlier.size=4, aes(x=year, y=rich_mod, colour=treatment)) +
  scale_color_manual(breaks = c("N","L","M","H","P"), values=c("red2", "dodgerblue","chartreuse3","darkorchid3","darkred"))

library("gridExtra")
library("grid")
library("cowplot")

plots <- list(p1, p2)

plot_layout <- rbind(c(1,2))


grid_layout <- plot_grid(p1, p2, align = "v", nrow = 2, rel_heights = c(1/4, 1/4), rel_widths = c(1/4, 1/4))

ggsave("./SI2-FigureS1_sample-coverage.png", grid_layout, height = 20, width = 30)


#-----------------------------------------------------------------------------------------------------------------------------------------------
