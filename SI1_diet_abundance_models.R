
#read data in from functional_abundance_plots.R

set.seed(1)

# DIET

### GENERALIST ABUNDANCE, no primary forest control

# glmmTMB

tmbg <- glmmTMB(generalist ~ 1 + (1 |site_id),
                family = genpois(link="log"),
                data = diet_matrix[diet_matrix$treatment!="P",])

tmbg0.5 <- glmmTMB(generalist ~ 1,
                   family = genpois(link="log"),
                   data = diet_matrix[diet_matrix$treatment!="P",])

tmbg1 <- glmmTMB(generalist ~ treatment + (1 |site_id),
                 family = genpois(link="log"),
                 data = diet_matrix[diet_matrix$treatment!="P",])

tmbg1.5 <- glmmTMB(generalist ~ treatment,
                   family = genpois(link="log"),
                   data = diet_matrix[diet_matrix$treatment!="P",])

tmbg2 <- glmmTMB(generalist ~ year + (1 |site_id),      
                 family = genpois(link="log"),
                 data = diet_matrix[diet_matrix$treatment!="P",])
summary(tmbg2)
emmeans(tmbg2, pairwise ~ year, type = "response")

tmbg2.5 <- glmmTMB(generalist ~ year,                #best model, 163.6738              
                   family = genpois(link="log"),
                   data = diet_matrix[diet_matrix$treatment!="P",])
#SI I: Table S9
summary(tmbg2.5)
emmeans(tmbg2.5, pairwise ~ year, type = "response")

tmbg3 <- glmmTMB(generalist ~ treatment + year + (1 |site_id),
                 family = genpois(link="log"),
                 data = diet_matrix[diet_matrix$treatment!="P",])

tmbg3.5 <- glmmTMB(generalist ~ treatment + year,
                   family = genpois(link="log"),
                   data = diet_matrix[diet_matrix$treatment!="P",])

tmbg4 <- glmmTMB(generalist ~ treatment*year + (1 |site_id),
                 family = genpois(link="log"),
                 data = diet_matrix[diet_matrix$treatment!="P",])

tmbg4.5 <- glmmTMB(generalist ~ treatment*year,
                   family = genpois(link="log"),
                   data = diet_matrix[diet_matrix$treatment!="P",])


AIC(tmbg, tmbg0.5, tmbg1, tmbg1.5, tmbg2, tmbg2.5, tmbg3, tmbg3.5, tmbg4, tmbg4.5)


#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------


# DIET

### COPROPHAGOUS ABUNDANCE, no primary forest control

# glmmTMB

tmbc <- glmmTMB(coprophagous ~ 1 + (1 |site_id),                 #EQUAL BEST MODEL, 749.5514
                family = genpois(link="log"),
                data = diet_matrix[diet_matrix$treatment!="P",])
summary(tmbc)

tmbc0.5 <- glmmTMB(coprophagous ~ 1,                            #BEST MODEL, 749.4125
                   family = genpois(link="log"),
                   data = diet_matrix[diet_matrix$treatment!="P",])
#SI I: Table S9
summary(tmbc0.5)

tmbc1 <- glmmTMB(coprophagous ~ treatment + (1 |site_id),
                 family = genpois(link="log"),
                 data = diet_matrix[diet_matrix$treatment!="P",])

tmbc1.5 <- glmmTMB(coprophagous ~ treatment,
                   family = genpois(link="log"),
                   data = diet_matrix[diet_matrix$treatment!="P",])

tmbc2 <- glmmTMB(coprophagous ~ year + (1 |site_id),      
                 family = genpois(link="log"),
                 data = diet_matrix[diet_matrix$treatment!="P",])
summary(tmbc2)
emmeans(tmbc2, pairwise ~ year, type = "response")

tmbc2.5 <- glmmTMB(coprophagous ~ year,               
                   family = genpois(link="log"),
                   data = diet_matrix[diet_matrix$treatment!="P",])
summary(tmbc2.5)
emmeans(tmbc2.5, pairwise ~ year, type = "response")

tmbc3 <- glmmTMB(coprophagous ~ treatment + year + (1 |site_id),
                 family = genpois(link="log"),
                 data = diet_matrix[diet_matrix$treatment!="P",])

tmbc3.5 <- glmmTMB(coprophagous ~ treatment + year,
                   family = genpois(link="log"),
                   data = diet_matrix[diet_matrix$treatment!="P",])

tmbc4 <- glmmTMB(coprophagous ~ treatment*year + (1 |site_id),
                 family = genpois(link="log"),
                 data = diet_matrix[diet_matrix$treatment!="P",])

tmbc4.5 <- glmmTMB(coprophagous ~ treatment*year,
                   family = genpois(link="log"),
                   data = diet_matrix[diet_matrix$treatment!="P",])


AIC(tmbc, tmbc0.5, tmbc1, tmbc1.5, tmbc2, tmbc2.5, tmbc3, tmbc3.5, tmbc4, tmbc4.5)
#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------


# DIET

### NECROPHAGOUS ABUNDANCE, no primary forest control

# glmmTMB

#low sample size; not used

tmbn <- glmmTMB(necrophagous ~ 1 + (1 |site_id),                
                family = genpois(link="log"),
                data = diet_matrix[diet_matrix$treatment!="P",])
summary(tmbn)

tmbn0.5 <- glmmTMB(necrophagous ~ 1,                            
                   family = genpois(link="log"),
                   data = diet_matrix[diet_matrix$treatment!="P",])
summary(tmbn0.5)

tmbn1 <- glmmTMB(necrophagous ~ treatment + (1 |site_id),
                 family = genpois(link="log"),
                 data = diet_matrix[diet_matrix$treatment!="P",])

tmbn1.5 <- glmmTMB(necrophagous ~ treatment,
                   family = genpois(link="log"),
                   data = diet_matrix[diet_matrix$treatment!="P",])

tmbn2 <- glmmTMB(necrophagous ~ year + (1 |site_id),      
                 family = genpois(link="log"),
                 data = diet_matrix[diet_matrix$treatment!="P",])
summary(tmbn2)
emmeans(tmbn2, pairwise ~ year, type = "response")

tmbn2.5 <- glmmTMB(necrophagous ~ year,               
                   family = genpois(link="log"),
                   data = diet_matrix[diet_matrix$treatment!="P",])
summary(tmbn2.5)
emmeans(tmbn2.5, pairwise ~ year, type = "response")


tmbn3 <- glmmTMB(necrophagous ~ treatment + year + (1 |site_id),
                 family = genpois(link="log"),
                 data = diet_matrix[diet_matrix$treatment!="P",])

tmbn3.5 <- glmmTMB(necrophagous ~ treatment + year,
                   family = genpois(link="log"),
                   data = diet_matrix[diet_matrix$treatment!="P",])

tmbn4 <- glmmTMB(necrophagous ~ treatment*year + (1 |site_id),
                 family = genpois(link="log"),
                 data = diet_matrix[diet_matrix$treatment!="P",])

tmbn4.5 <- glmmTMB(necrophagous ~ treatment*year,
                   family = genpois(link="log"),
                   data = diet_matrix[diet_matrix$treatment!="P",])


AIC(tmbn, tmbn0.5, tmbn1, tmbn1.5, tmbn2, tmbn2.5, tmbn3, tmbn3.5, tmbn4, tmbn4.5)
#-------------------------------------------------------------------------------------------------------------------------
