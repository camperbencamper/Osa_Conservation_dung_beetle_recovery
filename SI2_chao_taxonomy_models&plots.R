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

#chao-jaccard
p_dissim <- vegdist(comm_data_ordi[,names(comm_data_ordi) %in% DB_species], method="chao", binary=FALSE)
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
#-------------------------------------------------------------------------------------------------------------------------

bc_bet <- glmmTMB(similarity ~ 1,
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = p_sim)


bc_bet1 <- glmmTMB(similarity ~ 1 + (1|site_id),
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = p_sim)


bc_bet2 <- glmmTMB(similarity ~ year,
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = p_sim)
summary(bc_bet2)
emmeans::emmeans(bc_bet2, pairwise ~ year, type = "response", data=p_sim) #post-hoc
summary(b_posthoc2)
plot(b_posthoc2, comparisons = TRUE)
histogram(residuals(bc_bet2))
plot(residuals(bc_bet2))
testZeroInflation(bc_bet2)
testDispersion(bc_bet2)
testOverdispersion(bc_bet2)
plot(simulateResiduals(bc_bet2, n=500))

ggplot(posthoc, aes(x=year, y= emmean))


bc_bet3 <- glmmTMB(similarity ~ year + (1|site_id),
                   family = ordbeta(),
                   start=list(psi = c(-1, 1)),
                   data = p_sim)


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

### creating jaccard similarity matrix
j_dissim <- vegdist(comm_data_ordi[,names(comm_data_ordi) %in% DB_species], method="chao", binary=TRUE)
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
j_sim$year <- factor(j_sim$year, levels=c("2017","2019","2021"))


### JACCARD

j_bet <- glmmTMB(similarity ~ 1,
                 family = ordbeta(),
                 start=list(psi = c(-1, 1)),
                 data = j_sim)


j_bet1 <- glmmTMB(similarity ~ 1 + (1|site_id),
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = j_sim)


j_bet2 <- glmmTMB(similarity ~ year,
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = j_sim)
summary(j_bet2)
j_posthoc2 <- emmeans::emmeans(j_bet2, pairwise ~ year, type = "response", data=j_sim) #post-hoc
plot(j_posthoc2, comparisons = TRUE)
histogram(residuals(j_bet2))
plot(residuals(j_bet2))
testZeroInflation(j_bet2)
testDispersion(j_bet2)
testOverdispersion(j_bet2)
plot(simulateResiduals(j_bet2, n=500))


j_bet3 <- glmmTMB(similarity ~ year + (1|site_id), 
                  family = ordbeta(),
                  start=list(psi = c(-1, 1)),
                  data = j_sim)


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

#----------------------------------------------------------------------------------------------------------------

### PLOTTING

#SI II: Figure 2, panels A and B

#pull in functional model script: 'DB_chao_functional_models&plots.R"

#abundance
p1 <- ggplot(data=p_sim, aes(x=year, y=similarity, colour=treatment)) + ylab("Chao-Jaccard Similarity (Abundance)") + xlab("Year") + ggtitle("A)") +
  geom_boxplot(position = position_dodge(width = 0.9), size=1.75, outlier.size=4) + theme_bw(40) + theme(legend.position = "none") + ylim(0, 0.9) +
  stat_summary(fun = median, geom = 'line', linewidth=1.75, aes(group = treatment),position = position_dodge(width = 0.9)) +
  scale_color_manual(breaks = c("N","L","M","H"), values=c("red2", "dodgerblue","chartreuse3","darkorchid3"))

#incidence
p2 <- ggplot(data=j_sim, aes(x=year, y=similarity, colour=treatment)) + ylab("Chao-Jaccard Similarity (Incidence)") + xlab("Year") + ggtitle("B)") +
  geom_boxplot(position = position_dodge(width = 0.9), size=1.75, outlier.size=4) + theme_bw(40) + theme(legend.position = "none") + ylim(0, 0.9) +
  stat_summary(fun = median, geom = 'line', linewidth=1.75, aes(group = treatment),position = position_dodge(width = 0.9))+
  scale_color_manual(breaks = c("N","L","M","H"), values=c("red2", "dodgerblue","chartreuse3","darkorchid3"))
