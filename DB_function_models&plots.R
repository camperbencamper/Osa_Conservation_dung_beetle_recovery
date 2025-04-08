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
library("glmmTMB")
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


#---------------------------------------------------------------------------------------------------------

f_traits <- c("tunneler", "roller", "unknown", "dweller", "small", "medium", "large", "generalist", "coprophagous", "necrophagous") 
metadata <- c("treatment", "year", "treatment_year", "site_id", "site_year")

#------------------------------------------------------------------------------------------------------------------------------------------------

#restoration analysis; all functional traits

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


#------------------------------------------------------------------------------------------------------------------

bc_bet <- glmmTMB(similarity ~ 1,
               family = ordbeta(),
               start=list(psi = c(-1, 1)),
               data = p_sim)


bc_bet1 <- glmmTMB(similarity ~ 1 + (1|site_id),
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = p_sim)


bc_bet2 <- glmmTMB(similarity ~ year, #SECOND BEST MODEL; -345.78
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = p_sim)
summary(bc_bet2)
emmeans::emmeans(bc_bet2, pairwise ~ year, type = "response", data=p_sim) #post-hoc
histogram(residuals(bc_bet2))
plot(residuals(bc_bet2))
testZeroInflation(bc_bet2)
testDispersion(bc_bet2)
testOverdispersion(bc_bet2)
plot(simulateResiduals(bc_bet2, n=500))


bc_bet3 <- glmmTMB(similarity ~ year + (1|site_id), #BEST MODEL; -346.7
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = p_sim)
summary(bc_bet3)
bc_posthoc3 <- emmeans::emmeans(bc_bet3, pairwise ~ year, type = "response", data=p_sim) #post-hoc
plot(bc_posthoc3, comparisons = TRUE)
histogram(residuals(bc_bet3))
plot(residuals(bc_bet3))
testZeroInflation(bc_bet3)
testDispersion(bc_bet3)
testOverdispersion(bc_bet3)
plot(simulateResiduals(bc_bet3, n=500))


bc_bet4 <- glmmTMB(similarity ~ treatment,
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = p_sim)


bc_bet5 <- glmmTMB(similarity ~ treatment + (1|site_id),
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = p_sim)


bc_bet6 <- glmmTMB(similarity ~ treatment + year,
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = p_sim)


bc_bet7 <- glmmTMB(similarity ~ treatment + year + (1|site_id),
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = p_sim)


bc_bet8 <- glmmTMB(similarity ~ treatment*year,
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = p_sim)


bc_bet9 <- glmmTMB(similarity ~ treatment*year + (1|site_id),
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = p_sim)


AIC(bc_bet, bc_bet1, bc_bet2, bc_bet3, bc_bet4, bc_bet5, bc_bet6, bc_bet7, bc_bet8, bc_bet9)

#---------------------------------------------------------------

### JACCARD

j_bet <- glmmTMB(similarity ~ 1,
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = j_sim)


j_bet1 <- glmmTMB(similarity ~ 1 + (1|site_id),
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = j_sim)


j_bet2 <- glmmTMB(similarity ~ year, # SECOND BEST MODEL; -464.58
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = j_sim)
summary(j_bet2)
emmeans::emmeans(j_bet2, pairwise ~ year, type = "response", data=j_sim) #post-hoc
histogram(residuals(j_bet2))
plot(residuals(j_bet2))
testZeroInflation(j_bet2)
testDispersion(j_bet2)
testOverdispersion(j_bet2)
plot(simulateResiduals(j_bet2, n=500))


j_bet3 <- glmmTMB(similarity ~ year + (1|site_id), # BEST MODEL; -465.24
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = j_sim)
summary(j_bet3)
j_posthoc3 <- emmeans::emmeans(j_bet3, pairwise ~ year, type = "response", data=j_sim) #post-hoc
plot(j_posthoc3, comparisons = TRUE)
histogram(residuals(j_bet3))
plot(residuals(j_bet3))
testZeroInflation(j_bet3)
testDispersion(j_bet3)
testOverdispersion(j_bet3)
plot(simulateResiduals(j_bet3, n=500))
emmip(j_bet3, ~ year, CIs = TRUE, type = "response") +
  geom_point(aes(x = year, y = similarity), data = p_sim, pch = 2, color = "blue")


j_bet4 <- glmmTMB(similarity ~ treatment,
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = j_sim)


j_bet5 <- glmmTMB(similarity ~ treatment + (1|site_id),
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = j_sim)


j_bet6 <- glmmTMB(similarity ~ treatment + year,
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = j_sim)


j_bet7 <- glmmTMB(similarity ~ treatment + year + (1|site_id),
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = j_sim)


j_bet8 <- glmmTMB(similarity ~ treatment*year,
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = j_sim)


j_bet9 <- glmmTMB(similarity ~ treatment*year + (1|site_id),
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = j_sim)


AIC(j_bet, j_bet1, j_bet2, j_bet3, j_bet4, j_bet5, j_bet6, j_bet7, j_bet8, j_bet9)


#------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------


#Figure 3

#panels p1 (A) and p2 (B) pulled from "DB_taxonomy_models&plots.R"

#bray-curtis
p3 <- ggplot(data=p_sim, aes(x=year, y=similarity, colour=treatment)) + ylab("Bray-Curtis Similarity") + xlab("Year") + ggtitle("A)") +
  geom_boxplot(position = position_dodge(width = 0.9), size=1.75, outlier.size=4) + theme_bw(40) + theme(legend.position = "none") + ylim(0, 0.4) +
  stat_summary(fun = median, geom = 'line', linewidth=1.75, aes(group = treatment),position = position_dodge(width = 0.9)) +
  scale_color_manual(breaks = c("N","L","M","H"), values=c("red2", "dodgerblue","chartreuse3","darkorchid3"))

#jaccard
p4 <- ggplot(data=j_sim, aes(x=year, y=similarity, colour=treatment)) + ylab("Jaccard Similarity") + xlab("Year") + ggtitle("B)") +
  geom_boxplot(position = position_dodge(width = 0.9), size=1.75, outlier.size=4) + theme_bw(40) + theme(legend.position = "none") + ylim(0, 0.4) +
  stat_summary(fun = median, geom = 'line', linewidth=1.75, aes(group = treatment),position = position_dodge(width = 0.9))+
  scale_color_manual(breaks = c("N","L","M","H"), values=c("red2", "dodgerblue","chartreuse3","darkorchid3"))


library("gridExtra")
library("grid")
library("cowplot")

plots <- list(p1, p2, p3, p4)

plot_layout <- rbind(c(1,2),
                     c(3,4))

grid_layout <- plot_grid(p1, p2, p3, p4, align = "v", nrow = 2, rel_heights = c(1/4, 1/4, 1/4, 1/4), rel_widths = c(1/4, 1/4, 1/4, 1/4))

ggsave("./Figure3.png", grid_layout, height = 20, width = 30)


#----------------------------------------------------------------------------------------------------------------------------------

#richness and abundance plotting

#Figure 2: panel C

#panels p1 (A) and p2 (B) pulled from "DB_taxonomy_models&plots.R"

#richness
p3 <- ggplot() + theme_bw(40) + ylab("Functional Richness") + xlab("Year") + ggtitle("C)") + theme(legend.position = "none") +
  geom_boxplot(data=trait_df[trait_df$treatment!="P",], size=1.75, outlier.size=4, aes(x=as.factor(year), y=richness, colour=treatment)) +
  #geom_boxplot(data=trait_df[trait_df$treatment=="P",], size=1.75, outlier.size=4, aes(x=as.factor(year), y=richness, colour=treatment)) +
  scale_color_manual(breaks = c("N","L","M","H","P"), values=c("red2", "dodgerblue","chartreuse3","darkorchid3","darkred"))


library("gridExtra")
library("grid")
library("cowplot")

plots <- list(p1, p2, p3)

plot_layout <- rbind(c(1,2),
                     c(3,NA))


grid_layout <- plot_grid(p1, p2, p3, align = "v", nrow = 2, rel_heights = c(1/4, 1/4, 1/4), rel_widths = c(1/4, 1/4, 1/4))

ggsave("./Figure2.png", grid_layout, height = 20, width = 30)

