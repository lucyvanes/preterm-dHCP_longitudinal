# 27/01/2023
# Analysis script for Vanes et al. (2022)
# Longitudinal neonatal brain development and socio-demographic correlates of infant outcomes following preterm birth
# https://www.medrxiv.org/content/10.1101/2022.08.11.22278469v1

# Overview
#==============
# 1) PCA on 18-month neurodevelopmental outcomes 
# 2) Analysis of socio-demographic factors and behaviour
# 3) Analysis and plots for structural and functional brain-behaviour analysis


#===============================================================================
#    
#    1)          PCA on 18-month neurodevelopmental outcomes 
# 
#===============================================================================

# read and prep data
#=====================
setwd("C:/Users/k1456336/Dropbox/Github/preterm-dHCP_longitudinal/data")
dat <- read.csv("dHCP_behav.csv", header=T)
dat$id <- factor(dat$id)

cbcl_vars <- c( "cbcl_emot_react_tscore","cbcl_anx_dep_tscore","cbcl_somatic_tscore", "cbcl_withdrawn_tscore",
                "cbcl_sleep_tscore", "cbcl_attention_tscore","cbcl_aggessive_tscore")
bayley_vars <- c("cog_comp","language_comp","motor_comp")
psych_vars <- c(cbcl_vars, bayley_vars)

control_vars <- c("age_at_assess")

all_vars <- c("id",psych_vars, control_vars)

dat_prep <- dat[,all_vars]
dat_prep <- dat_prep[complete.cases(dat_prep),]


# Regress out age (from cbcl)
#============================
for (v in cbcl_vars){
  lm1 <- lm(paste(v, "~ age_at_assess"), dat_prep)
  dat_prep[v] <- lm1$resid
}
dat_prep$age <- NULL


# Run PCA
#================
pcadat <- dat_prep
pcadat$id <- NULL
pcadat_scaled <- as.data.frame(scale(pcadat, center=T, scale=T))
dim(pcadat_scaled)

# permutation testing
#======================
# for each permutation, shuffle rows in each column, re-compute PCA
# http://bioops.info/2015/01/permutation-pca/

# the fuction to assess the significance of the principal components.
sign.pc<-function(x,R=5000,s=10, cor=T,...){
  pc.out<-princomp(x,cor=cor,...)
  pve=(pc.out$sdev^2/sum(pc.out$sdev^2))[1:s]
  pve.perm<-matrix(NA,ncol=s,nrow=R)
  for(i in 1:R){
    x.perm<-apply(x,2,sample)
    pc.perm.out<-princomp(x.perm,cor=cor,...)
    pve.perm[i,]=(pc.perm.out$sdev^2/sum(pc.perm.out$sdev^2))[1:s]
  }
  pval<-apply(t(pve.perm)>pve,1,sum)/R
  return(list(pve=pve,pval=pval))
}

pca_sign <- sign.pc(pcadat_scaled,cor=T)
pca_sign

# first 2 PCs are significant (p < 0.05)

# Inspect PCA loadings
#==========================

out <- princomp(pcadat_scaled)
summary(out)

PCs_orig <- as.data.frame(out$scores)
PC_nos <- which(pca_sign$pval < 0.05)
PC_nos

PCs_orig_sig <- PCs_orig[,PC_nos]
names(PCs_orig_sig) <- paste("orig_PC", PC_nos, sep="")

loadings <- out$loadings
loadings <- as.data.frame.matrix(loadings)
loadings

# # http://strata.uga.edu/8370/lecturenotes/principalComponents.html
threshold <- sqrt(1/ncol(pcadat))  # cutoff for 'important' loadings
threshold 

loadings[abs(loadings) < threshold] <- NA
loadings


# Add PCs to dataset
#==================================
dat_prep$PC1 <- out$scores[,1]
dat_prep$PC2 <- out$scores[,2]


dat$id <- factor(dat$id)
dat_prep$id <- factor(dat_prep$id)

dat$PC1 <- NA
dat$PC2 <- NA

for (i in levels(dat_prep$id)){
  dat$PC1[dat$id==i] <- dat_prep$PC1[dat_prep$id==i]
  dat$PC2[dat$id==i] <- dat_prep$PC2[dat_prep$id==i]
}


# Plot heatmaps
#==================

library(corrplot)
library(RColorBrewer)
library(dplyr)

# prepare correlations
pcs_indx <- which(!is.na(match(names(dat_prep), c("PC1","PC2"))))
vars_indx <- which(!is.na(match(names(dat_prep), c(psych_vars))))

n_p <- length(pcs_indx)
n_v <- length(vars_indx)

cor_data <- dat_prep[,c(pcs_indx,vars_indx)]
M <- cor(cor_data)
M <- M[(n_p+1):(n_v+n_p), 1:n_p]

dimnames(M)[[1]] <- c(
  # "Surgency","Negative affect","Effortful control",
  # "QCHAT total",
  "Emotional reactivity","Anxiety/depression","Somatic","Withdrawn","Sleep","Attention","Aggressive",
  # "eye Social","eye orient","eye memory","eye cont","eye nonsoc",
  "Cognitive","Language","Motor")


res1 <- cor.mtest(cor_data, conf.level = .95)
res1$p <- res1$p[(n_p+1):(n_v+n_p), 1:n_p]
res1$lowCI <- res1$lowCI[(n_p+1):(n_v+n_p), 1:n_p]
res1$uppCI <- res1$uppCI[(n_p+1):(n_v+n_p), 1:n_p]

loadings <- as.matrix(loadings[,1:2])

col1=rev(colorRampPalette(brewer.pal(n=11, name="RdBu"))(100))
col2=colorRampPalette(brewer.pal(n=9, name="YlOrRd"))(100)


corrplot(M, p.mat = res1$p,method="color", col=col1,
         tl.col="black",tl.cex=1, # tl.offset=0.5,
         cl.align="l",cl.ratio = 1, 
         insig = "blank",
         addgrid.col="grey")

#===============================================================================
#    
#    2)          Analysis of socio-demographic factors and behaviour
# 
#===============================================================================

setwd("C:/Users/k1456336/Dropbox/Github/preterm-dHCP_longitudinal/data")
dat <- read.csv("dHCP_behav.csv", header=T)
dat$id <- factor(dat$id)

lm1 <- lm(PC1 ~ GA + age_at_assess + gender + imd + sps + PC2, dat)
lm2 <- lm(PC2 ~ GA + age_at_assess + gender + imd + sps + PC1, dat)

summary(lm1)
summary(lm2)

library(jtools)
library(ggplot2)

# Figure 5A (Psychopathology)
#=================================
effect_plot(lm1, pred=imd, plot.points=T, interval=T, x.label="Index of Multiple Deprivation", y.label="PC1 (partial residuals)",
            point.size=2, point.alpha=0.5, partial.residuals=T) +
  theme(axis.text=element_text(size=16, face="bold"), axis.title=element_text(size=20,face="bold"))



effect_plot(lm1, pred=sps, plot.points=T, interval=T, x.label="Stimulating Parenting Scale", y.label="PC1 (partial residuals)",
            point.size=2, point.alpha=0.5, partial.residuals=T) +
  theme(axis.text=element_text(size=16, face="bold"), axis.title=element_text(size=20,face="bold"))


# Figure 5B (Psychomotor functioning)
#====================================
effect_plot(lm2, pred=GA, plot.points=T, interval=T, x.label="Gestational age at birth", y.label="PC2 (partial residuals)",
            point.size=2, point.alpha=0.5, partial.residuals=T) +
  theme(axis.text=element_text(size=16, face="bold"), axis.title=element_text(size=20,face="bold"))


effect_plot(lm2, pred=imd, plot.points=T, interval=T, x.label="Index of Multiple Deprivation", y.label="PC2 (partial residuals)",
            point.size=2, point.alpha=0.5, partial.residuals=T) +
  theme(axis.text=element_text(size=16, face="bold"), axis.title=element_text(size=20,face="bold"))



#====================================================================================
#    
#    3)  Analysis and plots for structural and functional brain-behaviour analysis
# 
#====================================================================================

# Preparation and analysis scripts for FSL-SwE can be found under the "SwE_analyses" folders.
# Following SwE analyses, mean values from significant clusters are extracted using fslmeants and inserted into the design data for visualisation


# Median split function for visualisations
#===========================================
asfactor_cat2 <- function(var){
  var_median <- median(var, na.rm=T)
  data <- NULL
  data$var <- var
  data$cat_var <- as.numeric()
  data$cat_var[data$var <= var_median] <- "LOW"
  data$cat_var[data$var > var_median] <- "HIGH"
  
  data$cat_var <- factor(data$cat_var, levels=c("LOW","HIGH"))
  return(data$cat_var)
}

#===========================
# STRUCTURAL ANALYSIS PLOTS
#============================
setwd("C:/Users/k1456336/Dropbox/Github/preterm-dHCP_longitudinal/data")

dat <- read.csv("struct_longitudinal_data_bothPCs_IMD_SPS.csv", header=T)
dat$id <- factor(dat$id)

dat$PC1 <- dat$PC1 * -1
dat$PC1_cont <- dat$PC1
dat$Psychopathology <- asfactor_cat2(dat$PC1)
dat$PC2_cont <- dat$PC2

library(nlme)
library(ggplot2)
library(ggthemes)

# TSTAT3: positive effect of  PC2
# ==================================

# CLUSTER 1 (right cerebellum)
lme1 <- lme(tstat3_cluster1 ~ imd_c + sps_c + sex_c + day_bias_c + ga_c + pma_crosssectional + pma_longitudinal * PC1_cont+ pma_longitudinal * PC2_cont, dat, ~1|id,method="ML")
summary(lme1)

dat$tstat3_cluster1_pred <-predict(lme1, level=0, newdata=dat)

ggplot(dat, aes(x=PC2_cont, y=tstat3_cluster1, group=id)) + geom_point(size=2,col="grey",colour=NULL) + geom_line(col="grey",colour=NULL) +
  geom_smooth(aes(x=PC2_cont, y=tstat3_cluster1_pred, group=NULL), method='lm', size=2) +
  theme_stata() + scale_color_stata() + 
  theme(plot.title = element_text(size = 20, face = "bold"), text = element_text(size = 20),legend.position="right") + 
  xlab("Psychomotor functioning (PC2)") + ylab("Log-Jacobian") + ggtitle("Right Cerebellum") +
  guides(color = guide_legend(override.aes = list(size = 5)))


# TSTAT5: PC1 x PMA longitudinal
# ==============================

# CLUSTER  1 (R Thalamus)
lme1 <- lme(tstat5_cluster1 ~ imd_c + sps_c + sex_c + day_bias_c + ga_c + pma_crosssectional + pma_longitudinal * PC1_cont+ pma_longitudinal * PC2_cont, dat, ~1|id,method="ML")
summary(lme1)

dat$tstat5_cluster1_pred <-predict(lme1, level=0, newdata=dat)

ggplot(dat, aes(x=pma_longitudinal, y=tstat5_cluster1, group=id, colour=Psychopathology)) + geom_point(size=2,col="grey",colour=NULL) + geom_line(col="grey",colour=NULL) +
  geom_smooth(aes(x=pma_longitudinal, y=tstat5_cluster1_pred, group=NULL, colour=Psychopathology), method='lm', size=2) +
  labs(colour="PC1") +
  theme(plot.title = element_text(size = 20, face = "bold"), text = element_text(size = 20),legend.position="right") + 
  xlab("") + ylab("Log-Jacobian") + ggtitle("R Thalamus") +
  guides(color = guide_legend(override.aes = list(size = 5))) 



# CLUSTER  2 (bilateral STN)
lme2 <- lme(tstat5_cluster2 ~ imd_c + sps_c + sex_c + day_bias_c + ga_c + pma_crosssectional + pma_longitudinal * PC1_cont+ pma_longitudinal * PC2_cont, dat, ~1|id,method="ML")
summary(lme2)

dat$tstat5_cluster2_pred <-predict(lme2, level=0, newdata=dat)


ggplot(dat, aes(x=pma_longitudinal, y=tstat5_cluster2, group=id, colour=Psychopathology)) + geom_point(size=2,col="grey",colour=NULL) + geom_line(col="grey",colour=NULL) +
  geom_smooth(aes(x=pma_longitudinal, y=tstat5_cluster2_pred, group=NULL, colour=Psychopathology), method='lm', size=2) +
  labs(colour="PC1") +
  theme(plot.title = element_text(size = 20, face = "bold"), text = element_text(size = 20),legend.position="right") + 
  xlab("Time (centered, weeks)") + ylab("Log-Jacobian") + ggtitle("Bilateral STN") +
  guides(color = guide_legend(override.aes = list(size = 5))) 


#============================
# FUNCTIONAL ANALYSIS PLOTS
#============================
setwd("C:/Users/k1456336/Dropbox/Github/preterm-dHCP_longitudinal/data")

dat <- read.csv("func_longitudinal_data_bothPCs_IMD_SPS.csv", header=T)
dat$id <- factor(dat$id)


dat$PC1_cont <- dat$PC1
dat$PC2_cont <- dat$PC2
dat$Psychomotor <- asfactor_cat2(dat$PC2)

library(nlme)
library(ggplot2)
library(ggthemes)

# TSTAT 8: PC2 x PMA longitudinal
#=======================================

# CLUSTER 5 (left sensorimotor)
lme1 <- lme(tstat8_cluster5 ~ imd_c + sps_c + sex_c + motion + ga_c + pma_crosssectional + pma_longitudinal * PC1_cont+ pma_longitudinal * PC2_cont, dat, ~1|id,method="ML")
summary(lme1)

dat$tstat8_cluster5_pred <-predict(lme1, level=0, newdata=dat)

ggplot(dat, aes(x=pma_longitudinal, y=tstat8_cluster5, group=id, colour=Psychomotor)) + geom_point(size=2,col="grey",colour=NULL) + geom_line(col="grey",colour=NULL) +
  geom_smooth(aes(x=pma_longitudinal, y=tstat8_cluster5_pred, group=NULL, colour=Psychomotor), method='lm',size=2) +
  labs(colour="PC2") +
  theme(plot.title = element_text(size = 20, face = "bold"), text = element_text(size = 20),legend.position="bottom") + 
  xlab("Longitudinal PMA") + ylab("Degree Centrality") + ggtitle("Left sensorimotor cortex")  +
  guides(color = guide_legend(override.aes = list(size = 5)))







