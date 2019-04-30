library("survival")
library("survminer")
res.cox <- coxph(Surv(time, status) ~ age + sex + wt.loss, data =  lung)
res.cox
#To test for the proportional-hazards (PH) assumption, type this:
test.ph <- cox.zph(res.cox)
test.ph
# rho chisq     p
# age     -0.0483 0.378 0.538
# sex      0.1265 2.349 0.125
# wt.loss  0.0126 0.024 0.877
# GLOBAL       NA 2.846 0.416
# From the output above, the test is not statistically significant for each of the covariates, and the global test is also not statistically significant. Therefore, we can assume the proportional hazards.
# It’s possible to do a graphical diagnostic using the function ggcoxzph() [in the survminer package], which produces,
# for each covariate, graphs of the scaled Schoenfeld residuals against the transformed time.
ggcoxzph(test.ph)
ggcoxdiagnostics(res.cox, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())
ggcoxdiagnostics(res.cox, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())
ggcoxfunctional(Surv(time, status) ~ age + log(age) + sqrt(age), data = lung)
