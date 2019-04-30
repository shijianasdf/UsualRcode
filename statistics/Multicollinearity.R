library(tidyverse)
library(caret)
# Load the data
data("Boston", package = "MASS")
# Split the data into training and test set
set.seed(123)
training.samples <- Boston$medv %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- Boston[training.samples, ]
test.data <- Boston[-training.samples, ]
# Build the model
model1 <- lm(medv ~., data = train.data)
# Make predictions
predictions <- model1 %>% predict(test.data)
# Model performance
data.frame(
  RMSE = RMSE(predictions, test.data$medv),
  R2 = R2(predictions, test.data$medv)
)
car::vif(model1)
#In our example, the VIF score for the predictor variable tax is very high (VIF = 9.16). This might be problematic.
# Build a model excluding the tax variable
model2 <- lm(medv ~. -tax, data = train.data)
# Make predictions
predictions <- model2 %>% predict(test.data)
# Model performance
data.frame(
  RMSE = RMSE(predictions, test.data$medv),
  R2 = R2(predictions, test.data$medv)
)
# This chapter describes how to detect and deal with multicollinearity in regression models. Multicollinearity problems consist of including, in the model, different variables that have a similar predictive relationship with the outcome. This can be assessed for each predictor by computing the VIF value.
# Any variable with a high VIF value (above 5 or 10) should be removed from the model. This leads to a simpler model without compromising the model accuracy, which is good.
# Note that, in a large data set presenting multiple correlated predictor variables, you can perform principal component regression and partial least square regression strategies. See Chapter @ref(pcr-and-pls-regression).