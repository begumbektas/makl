#' Function to Train a Multiple Approximate Kernel Learning (MAKL) Model
#'
#'
#'
#'
#' Trains a MAKL model to be used as an input to makl_test() function.
#'
#'
#'
#'
#' @param X data matrix of size N x d, containing the training instances.
#' @param y response vector of length N, containing only -1 and 1.
#' @param D numeric value related to the number of random features to be used for approximation.
#' @param sigma_N numeric value preferably smaller than N, used to calculate sigma to create random features, it is a subset of the number of rows of X.
#' @param CV integer value between 0 and N. If CV is equal to 0 or 1, no cross validation is performed. If CV is greater than or equal to 2, CV is assigned as fold count in the cross validation.
#' @param lambda_set a continuous number between 0 and 1, used for regularization.
#' @param membership a list of length of number of groups, containing feature memberships to each group.
#'
#' @return a list containing the MAKL model and related parameters to be used in makl_test().
#' @export
#'
#' @examples
makl_train <- function(X, y, D = 100, sigma_N = 1000, CV = 1,
                       lambda_set = c(0.9, 0.8, 0.7, 0.6), membership) {

  N_group <- length(membership)
  feature_groups <- rep(1:N_group, each = 2 * D)
  feature_names <- sort(unique(unlist(sapply(1:N_group, FUN = function(x) {membership[[x]]}))))
  X <- X[, which(colnames(X) %in% feature_names)]

  pdist <- function(X1, X2) {
    if(identical(X1, X2) == TRUE) {
      D <- as.matrix(stats::dist(X1))
    }
    else {
      D <- as.matrix(stats::dist(rbind(X1, X2)))
      D <- D[1:nrow(X1), (nrow(X1) + 1):(nrow(X1) + nrow(X2))]
    }
    return(D)
  }

  if(CV >= 2) {
    fold_count <- CV
    negative_indices <- which(y == -1)
    positive_indices <- which(y == +1)

    negative_allocation <- sample(rep(1:fold_count, ceiling(length(negative_indices) / fold_count)), length(negative_indices))
    positive_allocation <- sample(rep(1:fold_count, ceiling(length(positive_indices) / fold_count)), length(positive_indices))

    auroc_matrix <- matrix(NA, nrow = fold_count, ncol = length(lambda_set), dimnames = list(1:fold_count, sprintf("%g", lambda_set)))

    for(fold in 1:fold_count) {
      train_indices <- c(negative_indices[which(negative_allocation != fold)], positive_indices[which(positive_allocation != fold)])
      test_indices <- c(negative_indices[which(negative_allocation == fold)], positive_indices[which(positive_allocation == fold)])
      X_train <- X[train_indices,]
      X_test <- X[test_indices,]
      valid_features <- as.numeric(which(apply(X_train, 2, sd) != 0))
      X_train <- X_train[, valid_features]
      X_test <- X_test[, valid_features]
      X_train <- scale(X_train)
      X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)

      Z_train <- matrix(0, nrow = nrow(X_train), ncol = 2 * D * N_group)
      Z_test <- matrix(0, nrow = nrow(X_test), ncol = 2 * D * N_group)

      for(m in 1:N_group) {
        feature_indices <- which(colnames(X_train) %in% membership[[m]])
        distance_indices <- sample(1:nrow(X_train), sigma_N)
        D_train <- pdist(X_train[distance_indices, feature_indices], X_train[distance_indices, feature_indices])
        sigma <- mean(D_train)

        gamma <- 1 / (2 * sigma^2)
        sigma_p <- sqrt(2 * gamma)

        W <- matrix(stats::rnorm(n = length(feature_indices) * D, mean = 0, sd = sigma_p), nrow = length(feature_indices), ncol = D)
        b <- stats::runif(D, min = 0, max = 2 * pi)
        Z <- sqrt(1 / D) * cbind(cos(X_train[, feature_indices] %*% W + matrix(b, nrow = nrow(X_train), ncol = D, byrow = TRUE)), sin(X_train[, feature_indices] %*% W + matrix(b, nrow = nrow(X_train), ncol = D, byrow = TRUE)))
        Z_train[,((m - 1) * 2 * D + 1):(m * 2 * D)] <- Z
        Z <- sqrt(1 / D) * cbind(cos(X_test[, feature_indices] %*% W + matrix(b, nrow = nrow(X_test), ncol = D, byrow = TRUE)), sin(X_test[, feature_indices] %*% W + matrix(b, nrow = nrow(X_test), ncol = D, byrow = TRUE)))
        Z_test[,((m - 1) * 2 * D + 1):(m * 2 * D)] <- Z
      }

      y_train <- y[train_indices]
      y_test <- y[test_indices]

      lambdamax <- grplasso::lambdamax(x = cbind(Z_train, 1), y = y_train / 2 + 0.5, index = c(feature_groups, NA), center = TRUE)
      lambda_set_val <- lambdamax * lambda_set
      model <- grplasso::grplasso(x = cbind(Z_train, 1), y = y_train / 2 + 0.5, index = c(feature_groups, NA), lambda = lambda_set_val, center = TRUE)
      y_predicted <- stats::predict(model, cbind(Z_test, 1))

      for(lambda in lambda_set) {
        k <- which(lambda == lambda_set)
        auroc_matrix[fold, k] <- AUC::auc(AUC::roc(predictions = y_predicted[, k], labels = as.factor(y_test / 2 + 0.5)))
      }
    }
    lambda_set <- lambda_set[max.col(t(colMeans(auroc_matrix, na.rm = TRUE)))]
  }

  valid_features <- as.numeric(which(apply(X, 2, sd) != 0))
  X_train <- X[, valid_features]
  X_train <- scale(X_train)
  y_train <- y
  mean <- attr(X_train, "scaled:center")
  sd <- attr(X_train, "scaled:scale")
  Z_train <- matrix(0, nrow = nrow(X_train), ncol = 2 * D * N_group)
  W <- vector("list", N_group)
  b <- vector("list", N_group)

  for(m in 1:N_group) {
    feature_indices <- which(colnames(X_train) %in% membership[[m]])
    distance_indices <- sample(1:nrow(X_train), sigma_N)
    D_train <- pdist(X_train[distance_indices, feature_indices], X_train[distance_indices, feature_indices])
    sigma <- mean(D_train)

    gamma <- 1 / (2 * sigma^2)
    sigma_p <- sqrt(2 * gamma)

    W[[m]] <- matrix(stats::rnorm(n = length(feature_indices) * D, mean = 0, sd = sigma_p), nrow = length(feature_indices), ncol = D)
    b[[m]] <- stats::runif(D, min = 0, max = 2 * pi)
    Z <- sqrt(1 / D) * cbind(cos(X_train[, feature_indices] %*% W[[m]] + matrix(b[[m]], nrow = nrow(X_train), ncol = D, byrow = TRUE)), sin(X_train[, feature_indices] %*% W[[m]] + matrix(b[[m]], nrow = nrow(X_train), ncol = D, byrow = TRUE)))
    Z_train[,((m - 1) * 2 * D + 1):(m * 2 * D)] <- Z
  }

  lambdamax <- grplasso::lambdamax(x = cbind(Z_train, 1), y = y_train / 2 + 0.5, index = c(feature_groups, NA), center = TRUE)
  lambda_set_val <- lambdamax * lambda_set

  model <- grplasso::grplasso(x = cbind(Z_train, 1), y = y_train / 2 + 0.5, index = c(feature_groups, NA), lambda = lambda_set_val, center = TRUE)

  makl_model <- list()
  makl_model$mean <- mean
  makl_model$sd <- sd
  makl_model$W <- W
  makl_model$b <- b
  makl_model$model <- model
  makl_model$valid_features <- valid_features
  makl_model$D <- D
  makl_model$membership <- membership
  makl_model$lambda_set <- lambda_set
  makl_model$lambda_set_val <- lambda_set_val

  return(makl_model)
}
