
#read data in from functional_abundance_plots.R

set.seed(1)

# SIZE

### SMALL ABUNDANCE, no primary forest control

# glmmTMB
tmbs <- glmmTMB(small ~ 1 + (1 |site_id),
                family = genpois(link="log"),
                data = size_matrix[size_matrix$treatment!="P",])

tmb0.5s <- glmmTMB(small ~ 1,
                family = genpois(link="log"),
                data = size_matrix[size_matrix$treatment!="P",])

tmb1s <- glmmTMB(small ~ treatment + (1 |site_id),
        family = genpois(link="log"),
        data = size_matrix[size_matrix$treatment!="P",])

tmb1.5s <- glmmTMB(small ~ treatment,
                 family = genpois(link="log"),
                 data = size_matrix[size_matrix$treatment!="P",])

tmb2s <- glmmTMB(small ~ year + (1 |site_id),       #equally best model, 341.1020
        family = genpois(link="log"),
        data = size_matrix[size_matrix$treatment!="P",])
#SI I: Table S7
summary(tmb2s)
#SI I: Table S8
emmeans(tmb2s, pairwise ~ year, type = "response")

tmb2.5s <- glmmTMB(small ~ year,                              #BesT MODEL, 339.2723
                 family = genpois(link="log"),
                 data = size_matrix[size_matrix$treatment!="P",])
#SI I: Table S7
summary(tmb2.5s)
#SI I: Table S8
emmeans(tmb2.5s, pairwise ~ year, type = "response")

tmb3s <- glmmTMB(small ~ treatment + year + (1 |site_id),
        family = genpois(link="log"),
        data = size_matrix[size_matrix$treatment!="P",])

tmb3.5s <- glmmTMB(small ~ treatment + year,
                    family = genpois(link="log"),
                    data = size_matrix[size_matrix$treatment!="P",])
  
tmb4s <- glmmTMB(small ~ treatment*year + (1 |site_id),
        family = genpois(link="log"),
        data = size_matrix[size_matrix$treatment!="P",])

tmb4.5s <- glmmTMB(small ~ treatment*year,
                 family = genpois(link="log"),
                 data = size_matrix[size_matrix$treatment!="P",])


AIC(tmbs, tmb0.5s, tmb1s, tmb1.5s, tmb2s, tmb2.5s, tmb3s, tmb3.5s, tmb4s, tmb4.5s)

#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------


# SIZE

### MEDIUM ABUNDANCE, no primary forest control

# glmmTMB

tmbm <- glmmTMB(medium ~ 1 + (1 |site_id),
                family = genpois(link="log"),
                data = size_matrix[size_matrix$treatment!="P",])

tmb0.5m <- glmmTMB(medium ~ 1,
                   family = genpois(link="log"),
                   data = size_matrix[size_matrix$treatment!="P",])

tmb1m <- glmmTMB(medium ~ treatment + (1 |site_id),
                 family = genpois(link="log"),
                 data = size_matrix[size_matrix$treatment!="P",])

tmb1.5m <- glmmTMB(medium ~ treatment,
                   family = genpois(link="log"),
                   data = size_matrix[size_matrix$treatment!="P",])

tmb2m <- glmmTMB(medium ~ year + (1 |site_id),                  #equally best model, 675.6887
                 family = genpois(link="log"),
                 data = size_matrix[size_matrix$treatment!="P",])
#SI I: Table S7
summary(tmb2m)
#SI I: Table S8
emmeans(tmb2m, pairwise ~ year, type = "response")

tmb2.5m <- glmmTMB(medium ~ year,                              #BEST MODEL, 674.0532
                   family = genpois(link="log"),
                   data = size_matrix[size_matrix$treatment!="P",])
#SI I: Table S7
summary(tmb2.5m)
#SI I: Table S8
emmeans(tmb2.5m, pairwise ~ year, type = "response")

tmb3m <- glmmTMB(medium ~ treatment + year + (1 |site_id),
                 family = genpois(link="log"),
                 data = size_matrix[size_matrix$treatment!="P",])

tmb3.5m <- glmmTMB(medium ~ treatment + year,
                   family = genpois(link="log"),
                   data = size_matrix[size_matrix$treatment!="P",])

tmb4m <- glmmTMB(medium ~ treatment*year + (1 |site_id),
                 family = genpois(link="log"),
                 data = size_matrix[size_matrix$treatment!="P",])

tmb4.5m <- glmmTMB(medium ~ treatment*year,
                   family = genpois(link="log"),
                   data = size_matrix[size_matrix$treatment!="P",])


AIC(tmbm, tmb0.5m, tmb1m, tmb1.5m, tmb2m, tmb2.5m, tmb3m, tmb3.5m, tmb4m, tmb4.5m)



#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------


# SIZE

### LARGE ABUNDANCE, no primary forest control

# glmmTMB

tmbl <- glmmTMB(large ~ 1 + (1 |site_id),
                family = genpois(link="log"),
                data = size_matrix[size_matrix$treatment!="P",])

tmb0.5l <- glmmTMB(large ~ 1,
                   family = genpois(link="log"),
                   data = size_matrix[size_matrix$treatment!="P",])

tmb1l <- glmmTMB(large ~ treatment + (1 |site_id),
                 family = genpois(link="log"),
                 data = size_matrix[size_matrix$treatment!="P",])

tmb1.5l <- glmmTMB(large ~ treatment,
                   family = genpois(link="log"),
                   data = size_matrix[size_matrix$treatment!="P",])

tmb2l <- glmmTMB(large ~ year + (1 |site_id),                  
                 family = genpois(link="log"),
                 data = size_matrix[size_matrix$treatment!="P",])

tmb2.5l <- glmmTMB(large ~ year,                              #BEST MODEL, 166.9272
                   family = genpois(link="log"),
                   data = size_matrix[size_matrix$treatment!="P",])
#SI I: Table S7
summary(tmb2.5l)
emmeans(tmb2.5l, pairwise ~ year, type = "response")

tmb3l <- glmmTMB(large ~ treatment + year + (1 |site_id),
                 family = genpois(link="log"),
                 data = size_matrix[size_matrix$treatment!="P",])

tmb3.5l <- glmmTMB(large ~ treatment + year,
                   family = genpois(link="log"),
                   data = size_matrix[size_matrix$treatment!="P",])

tmb4l <- glmmTMB(large ~ treatment*year + (1 |site_id),
                 family = genpois(link="log"),
                 data = size_matrix[size_matrix$treatment!="P",])

tmb4.5l <- glmmTMB(large ~ treatment*year,
                   family = genpois(link="log"),
                   data = size_matrix[size_matrix$treatment!="P",])


AIC(tmbl, tmb0.5l, tmb1l, tmb1.5l, tmb2l, tmb2.5l, tmb3l, tmb3.5l, tmb4l, tmb4.5l)

#-------------------------------------------------------------------------------------------------------------------------
