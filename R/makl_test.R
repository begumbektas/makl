#' Function for Binary Classification using a Multiple Approximate Kernel Learning (MAKL) Model
#'
#'
#'
#'
#' Classifies the test data, using the MAKL model resulted from makl_train().
#'
#'
#'
#'
#' @param X data matrix of size T x d, containing the test instances.
#' @param y response vector of length T, containing only -1 and 1.
#' @param makl_model a list containing the MAKL model returning from makl_train() function.
#'
#' @return a list containing the predictions for test instances and the area under the ROC curve (auroc) values with corresponding number of used kernels for prediction.
#' @export
#'
#' @examples
makl_test <- function(X, y, makl_model) {

  N_group <- length(makl_model$membership)
  feature_names <- sort(unique(unlist(sapply(1:N_group, FUN = function(x) {makl_model$membership[[x]]}))))
  X <- X[, which(colnames(X) %in% feature_names)]
  X_test <- X[, makl_model$valid_features]
  X_test <- (X_test - matrix(makl_model$mean, nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(makl_model$sd, nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)

  Z_test <- matrix(0, nrow = nrow(X_test), ncol = 2 * makl_model$D * N_group)
  for(m in 1:N_group) {
    feature_indices <- which(colnames(X_test) %in% makl_model$membership[[m]])
    Z <- sqrt(1 / makl_model$D) * cbind(cos(X_test[, feature_indices] %*% makl_model$W[[m]] + matrix(makl_model$b[[m]], nrow = nrow(X_test), ncol = makl_model$D, byrow = TRUE)), sin(X_test[, feature_indices] %*% makl_model$W[[m]] + matrix(makl_model$b[[m]], nrow = nrow(X_test), ncol = makl_model$D, byrow = TRUE)))
    Z_test[,((m - 1) * 2 * makl_model$D + 1):(m * 2 * makl_model$D)] <- Z
  }
  y_predicted <- stats::predict(makl_model$model, cbind(Z_test, 1))

  auroc_array <- array(0, dim = length(makl_model$lambda_set), dimnames = list(sprintf("%g", makl_model$lambda_set)))
  for(lambda in makl_model$lambda_set) {
    k <- which(makl_model$lambda_set == lambda)
    auroc_array[k] <- AUC::auc(AUC::roc(predictions = y_predicted[, k], labels = as.factor(y / 2 + 0.5)))
  }

  coefficient_matrix <- array(0, dim = c(2 * makl_model$D, N_group, length(makl_model$lambda_set)), dimnames = list(1:(2 * makl_model$D), 1:N_group, sprintf("%g", makl_model$lambda_set)))
  for(m in 1:N_group) {
    coefficient_matrix[,m,] <- makl_model$model$coefficients[((m - 1) * 2 * makl_model$D + 1):(m * 2 * makl_model$D),]
  }

  pathway_norms <- matrix(NA, nrow = N_group, length(makl_model$lambda_set))
  colnames(pathway_norms) <-  sprintf("%g", makl_model$lambda_set)
  n_selected_kernels <- array(NA, dim = length(makl_model$lambda_set))

  for(lambda in makl_model$lambda_set) {
    k <- which(makl_model$lambda_set == lambda)
    pathway_norms[, k] <- sapply(1:N_group, function(p) {sqrt(sum(coefficient_matrix[,p,k]^2))})
    n_selected_kernels[k] <- length(which(pathway_norms[,k] > 0))
  }

  auroc_vs_kernel_number <- cbind(auroc_array, n_selected_kernels)
  return(list(y_predicted = y_predicted, auroc_kernel_number = auroc_vs_kernel_number))

}
