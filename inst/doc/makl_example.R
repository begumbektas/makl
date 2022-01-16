## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(MAKL)
set.seed(64327) #midas
df <- matrix(rnorm(6000, 0, 1), nrow = 1000)
colnames(df) <- c("F1", "F2", "F3", "F4", "F5", "F6")

## -----------------------------------------------------------------------------
# check colnames(df) for them to be matching with group members
groups <- list()
groups[[1]] <- c("F1", "F5", "F6")
groups[[2]] <- c("F2", "F3", "F4")

## -----------------------------------------------------------------------------
y <- c()
for(i in 1:nrow(df)) {
  if((df[i, 2] + df[i, 3] + df[i, 4]) > 0) {
    y[i] <- +1
  } else {
    y[i] <- -1
  }
}

## -----------------------------------------------------------------------------
makl_model <- makl_train(X = df, y = y, D = 2, sigma_N = 1000, CV = 1, membership = groups, lambda_set = c(0.9, 0.8, 0.7, 0.6))

## -----------------------------------------------------------------------------
makl_model$model$coefficients

## -----------------------------------------------------------------------------
df_test <- matrix(rnorm(600, 0, 1), nrow = 100)
colnames(df_test) <- c("F1", "F2", "F3", "F4", "F5", "F6")
y_test <- c()
for(i in 1:nrow(df_test)) {
  if((df_test[i, 2] + df_test[i, 3] + df_test[i, 4]) > 0) {
    y_test[i] <- +1
  } else {
    y_test[i] <- -1
  }
}
result <-makl_test(X = df_test, y = y_test, makl_model = makl_model)

## -----------------------------------------------------------------------------
result$auroc_kernel_number

