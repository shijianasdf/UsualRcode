#!/usr/bin/env Rscript

library(data.table)
library(car)

##########################################################################################

# Load data: genomic features w/ N >=30 events
wgd.dat <- fread("C:/Users/lenovo/Desktop/GD-master/data/FeatureMatrix.txt")
dim(wgd.dat) # 9181 x 270
head(wgd.dat,5)
# WGD                                Cancer_Type_Merged ABL1_lof AKT1_hotspot AKT2_amp
# 1:   0 Renal Cell Carcinoma | Renal Clear Cell Carcinoma        0            0        0
# 2:   0                             Ovarian Cancer | Rare        0            0        0
# 3:   0                       Bone Cancer | Ewing Sarcoma        0            0        0
# APC_hotspot APC_lof AR_lof AR_amp ARID1A_hotspot ARID1A_lof ARID1B_lof ARID2_lof
# 1:           0       0      0      0              0          0          0         0
# 2:           0       0      0      0              0          0          0         0
# 3:           0       0      0      0              0          0          0         0
# ARID5B_lof ASXL1_lof ASXL1_amp ASXL2_lof ATM_hotspot ATM_lof ATR_lof ATRX_lof AURKA_amp
# 1:          0         0         0         0           0       0       0        0         0
# 2:          0         0         0         0           0       0       0        0         0
# 3:          0         0         0         0           0       0       0        0         0
# AXIN1_lof AXIN2_amp B2M_lof BAP1_lof BCL2L1_amp BCOR_lof BRAF_hotspot BRCA1_lof
# 1:         0         0       0        0          0        0            0         0
# 2:         0         0       0        0          0        0            0         0
# 3:         0         0       0        0          0        0            0         0
# BRCA2_lof BRIP1_amp BTK_lof CARD11_lof CARD11_amp CASP8_lof CBFB_lof CBL_lof CCND1_amp
# 1:         0         0       0          0          0         0        0       0         0
# 2:         0         0       0          0          0         0        0       0         0
# 3:         0         0       0          0          0         0        0       0         0
# CCND2_amp CCND3_amp CCNE1_amp CD274_lof CD274_amp CD79B_amp CDH1_hotspot CDH1_lof
# 1:         0         0         0         0         0         0            0        0
# 2:         0         0         0         0         0         0            0        0
# 3:         0         0         0         0         0         0            0        0
# Run logistic regression
model.full <-
  glm(WGD ~ ., family = binomial(link = 'logit'),
      data = wgd.dat)
summary(model.full)
# Check variance inflation factors
model.vif = as.data.table(car::vif(model.full))
model.vif[, Test := (`GVIF^(1/(2*Df))`) ^ 2]
any(model.vif$Test > 10) # FALSE

# Extract model variables
model.vars <-
  as.data.table(summary(model.full)$coefficients, 
                keep.rownames = T)

# Compute 95% confidence intervals for odds ratios
model.vars[, OR := exp(Estimate)]
model.vars[, OR_lower := exp(Estimate - 1.96 * `Std. Error`)]
model.vars[, OR_upper := exp(Estimate + 1.96 * `Std. Error`)]
model.vars = model.vars[order(`Pr(>|z|)`)]
