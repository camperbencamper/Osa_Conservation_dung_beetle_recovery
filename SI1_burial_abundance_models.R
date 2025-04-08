
#read data in from functional_abundance_plots.R

set.seed(1)
# burial

### tunneler ABUNDANCE, no primary forest control

# glmmTMB

tmbt <- glmmTMB(tunneler ~ 1 + (1 |site_id),
                family = genpois(link="log"),
                data = burial_matrix[burial_matrix$treatment!="P",])

tmbt0.5 <- glmmTMB(tunneler ~ 1,
                   family = genpois(link="log"),
                   data = burial_matrix[burial_matrix$treatment!="P",])

tmbt1 <- glmmTMB(tunneler ~ treatment + (1 |site_id),
                 family = genpois(link="log"),
                 data = burial_matrix[burial_matrix$treatment!="P",])

tmbt1.5 <- glmmTMB(tunneler ~ treatment,
                   family = genpois(link="log"),
                   data = burial_matrix[burial_matrix$treatment!="P",])

tmbt2 <- glmmTMB(tunneler ~ year + (1 |site_id),                  #equal best model, 637.7714
                 family = genpois(link="log"),
                 data = burial_matrix[burial_matrix$treatment!="P",])
#SI I: Table S10
summary(tmbt2)
#SI I: Table S12
emmeans(tmbt2, pairwise ~ year, type = "response")

tmbt2.5 <- glmmTMB(tunneler ~ year,                #best model, 635.8415
                   family = genpois(link="log"),
                   data = burial_matrix[burial_matrix$treatment!="P",])
#SI I: Table S10
summary(tmbt2.5)
#SI I: Table S12
emmeans(tmbt2.5, pairwise ~ year, type = "response")

tmbt3 <- glmmTMB(tunneler ~ treatment + year + (1 |site_id),
                 family = genpois(link="log"),
                 data = burial_matrix[burial_matrix$treatment!="P",])

tmbt3.5 <- glmmTMB(tunneler ~ treatment + year,
                   family = genpois(link="log"),
                   data = burial_matrix[burial_matrix$treatment!="P",])

tmbt4 <- glmmTMB(tunneler ~ treatment*year + (1 |site_id),
                 family = genpois(link="log"),
                 data = burial_matrix[burial_matrix$treatment!="P",])

tmbt4.5 <- glmmTMB(tunneler ~ treatment*year,
                   family = genpois(link="log"),
                   data = burial_matrix[burial_matrix$treatment!="P",])


AIC(tmbt, tmbt0.5, tmbt1, tmbt1.5, tmbt2, tmbt2.5, tmbt3, tmbt3.5, tmbt4, tmbt4.5)


#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------


# burial

### roller ABUNDANCE, no primary forest control

# glmmTMB


tmbr <- glmmTMB(roller ~ 1 + (1 |site_id),
                family = genpois(link="log"),
                data = burial_matrix[burial_matrix$treatment!="P",])

tmbr0.5 <- glmmTMB(roller ~ 1,
                   family = genpois(link="log"),
                   data = burial_matrix[burial_matrix$treatment!="P",])

tmbr1 <- glmmTMB(roller ~ treatment + (1 |site_id),
                 family = genpois(link="log"),
                 data = burial_matrix[burial_matrix$treatment!="P",])

tmbr1.5 <- glmmTMB(roller ~ treatment,
                   family = genpois(link="log"),
                   data = burial_matrix[burial_matrix$treatment!="P",])

tmbr2 <- glmmTMB(roller ~ year + (1 |site_id),                
                 family = genpois(link="log"),
                 data = burial_matrix[burial_matrix$treatment!="P",])
summary(tmbr2)
emmeans(tmbr2, pairwise ~ year, type = "response")

tmbr2.5 <- glmmTMB(roller ~ year,                #best model, 564.2734
                   family = genpois(link="log"),
                   data = burial_matrix[burial_matrix$treatment!="P",])
#SI I: Table S10
summary(tmbr2.5)
#SI I: Table S12
emmeans(tmbr2.5, pairwise ~ year, type = "response")

tmbr3 <- glmmTMB(roller ~ treatment + year + (1 |site_id),
                 family = genpois(link="log"),
                 data = burial_matrix[burial_matrix$treatment!="P",])

tmbr3.5 <- glmmTMB(roller ~ treatment + year,               #equal best model, 566.2658
                   family = genpois(link="log"),
                   data = burial_matrix[burial_matrix$treatment!="P",])
#SI I: Table S11
summary(tmbr3.5)

tmbr4 <- glmmTMB(roller ~ treatment*year + (1 |site_id),
                 family = genpois(link="log"),
                 data = burial_matrix[burial_matrix$treatment!="P",])

tmbr4.5 <- glmmTMB(roller ~ treatment*year,
                   family = genpois(link="log"),
                   data = burial_matrix[burial_matrix$treatment!="P",])


AIC(tmbr, tmbr0.5, tmbr1, tmbr1.5, tmbr2, tmbr2.5, tmbr3, tmbr3.5, tmbr4, tmbr4.5)

#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------


# burial

### dweller ABUNDANCE, no primary forest control

# glmmTMB

#small sample size; not used


tmbd <- glmmTMB(dweller ~ 1 + (1 |site_id),
                family = genpois(link="log"),
                data = burial_matrix[burial_matrix$treatment!="P",])

tmbd0.5 <- glmmTMB(dweller ~ 1,
                   family = genpois(link="log"),
                   data = burial_matrix[burial_matrix$treatment!="P",])

tmbd1 <- glmmTMB(dweller ~ treatment + (1 |site_id),
                 family = genpois(link="log"),
                 data = burial_matrix[burial_matrix$treatment!="P",])

tmbd1.5 <- glmmTMB(dweller ~ treatment,
                   family = genpois(link="log"),
                   data = burial_matrix[burial_matrix$treatment!="P",])

tmbd2 <- glmmTMB(dweller ~ year + (1 |site_id),                
                 family = genpois(link="log"),
                 data = burial_matrix[burial_matrix$treatment!="P",])
summary(tmbd2)
emmeans(tmbd2, pairwise ~ year, type = "response")

tmbd2.5 <- glmmTMB(dweller ~ year,               
                   family = genpois(link="log"),
                   data = burial_matrix[burial_matrix$treatment!="P",])
summary(tmbd2.5)
emmeans(tmbd2.5, pairwise ~ year, type = "response")

tmbd3 <- glmmTMB(dweller ~ treatment + year + (1 |site_id),
                 family = genpois(link="log"),
                 data = burial_matrix[burial_matrix$treatment!="P",])

tmbd3.5 <- glmmTMB(dweller ~ treatment + year,               
                   family = genpois(link="log"),
                   data = burial_matrix[burial_matrix$treatment!="P",])
summary(tmbd3.5)

tmbd4 <- glmmTMB(dweller ~ treatment*year + (1 |site_id),
                 family = genpois(link="log"),
                 data = burial_matrix[burial_matrix$treatment!="P",])

tmbd4.5 <- glmmTMB(dweller ~ treatment*year,
                   family = genpois(link="log"),
                   data = burial_matrix[burial_matrix$treatment!="P",])


AIC(tmbd, tmbd0.5, tmbd1, tmbd1.5, tmbd2, tmbd2.5, tmbd3, tmbd3.5, tmbd4, tmbd4.5)
#-------------------------------------------------------------------------------------------------------------------------
