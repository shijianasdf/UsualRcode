################################################################################################
########################################logistic回归############################################
#A logistic regression model allows us to establish a relationship between a binary outcome variable and a group of predictor variables. 
mydata <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
## view the first few rows of the data
head(mydata)
##   admit gre  gpa rank
## 1     0 380 3.61    3
## 2     1 660 3.67    3
## 3     1 800 4.00    1
## 4     1 640 3.19    4
## 5     0 520 2.93    4
## 6     1 760 3.00    2
summary(mydata);
mydata$rank <- factor(mydata$rank);
mylogit <- glm(admit ~ gre + gpa + rank, data = mydata, family = "binomial");
#逐步logistic回归
mylogit <- step(glm(admit ~ gre + gpa + rank, data = mydata, family = "binomial"),direction = "backward");
tmp <- summary(mylogit); #list
tmp;
#Call:
#glm(formula = admit ~ gre + gpa + rank, family = "binomial", 
#    data = mydata)
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.6268  -0.8662  -0.6388   1.1490   2.0790  

#Coefficients:
#             Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -3.989979   1.139951  -3.500 0.000465 ***
#gre          0.002264   0.001094   2.070 0.038465 *  
#gpa          0.804038   0.331819   2.423 0.015388 *  
#rank2       -0.675443   0.316490  -2.134 0.032829 *  
#rank3       -1.340204   0.345306  -3.881 0.000104 ***
#rank4       -1.551464   0.417832  -3.713 0.000205 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for binomial family taken to be 1)

#    Null deviance: 499.98  on 399  degrees of freedom
#Residual deviance: 458.52  on 394  degrees of freedom
#AIC: 470.52

#Number of Fisher Scoring iterations: 4
## odds ratios and 95% CI
exp(cbind.data.frame(OR = coef(mylogit), confint(mylogit)));
##P值
tmp$coefficients[,4];
##组合OR CI P值
cbind.data.frame(exp(cbind.data.frame(OR = coef(mylogit), confint(mylogit))),P=tmp$coefficients[,4]);
#########################################################################################################