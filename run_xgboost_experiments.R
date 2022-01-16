source("classification_helper.R")
library(AUC)
library(xgboost)

data_path <- "./data"
result_path <- "./results"
pathway <- "hallmark" #replace with "pid" if you would like to use PID pathways
problem <- "stage" #replace with "survival" or with "sc" if you would like to perform the experiments on these problems.

for(replication in 1:100) {
  if(problem == "stage") {
    if(file.exists(sprintf("%s/xgboost_pathway_%s_%s_measure_AUROC_replication_%d_result.RData", result_path, problem, pathway, replication)) == FALSE) {
      load(sprintf("%s/all_clinical_stage.RData", data_path))
      load(sprintf("%s/all_mrna_stage.RData", data_path))

      common_patients <- intersect(rownames(all_clinical_stage)[which(is.na(all_clinical_stage$pathologic_stage) == FALSE & all_clinical_stage$pathologic_stage != "Stage X" & all_clinical_stage$pathologic_stage != "IS")], rownames(all_mrna_stage))

      X <- log2(all_mrna_stage[common_patients,] + 1)
      y <- rep(NA, length(common_patients))

      y[all_clinical_stage[common_patients, "pathologic_stage"] %in% c("Stage I",  "Stage IA",  "Stage IB",  "Stage IC")] <- +1
      y[all_clinical_stage[common_patients, "pathologic_stage"] %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC",
                                                                       "Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC",
                                                                       "Stage IV",  "Stage IVA",  "Stage IVB",  "Stage IVC")] <- -1
    }
  }

  if(problem == "survival") {
    if(file.exists(sprintf("%s/xgboost_pathway_%s_%s_measure_AUROC_replication_%d_result.RData", result_path, problem, pathway, replication)) == FALSE) {
      load(sprintf("%s/all_clinical_survival.RData", data_path))
      load(sprintf("%s/all_mrna_survival.RData", data_path))

      patients_with_mrna <- rownames(all_mrna_survival)
      patients_with_survival <- rownames(all_clinical_survival)[which(is.na(all_clinical_survival$vital_status) == FALSE & ((all_clinical_survival$vital_status == "Dead" & is.na(all_clinical_survival$days_to_death) == FALSE & all_clinical_survival$days_to_death > 0) | (all_clinical_survival$vital_status == "Alive" & ((is.na(all_clinical_survival$days_to_last_followup) == FALSE & all_clinical_survival$days_to_last_followup >= 365 * 2) | (is.na(all_clinical_survival$days_to_last_known_alive) == FALSE & all_clinical_survival$days_to_last_known_alive >= 365 *2)))))]
      common_patients <- intersect(patients_with_mrna, patients_with_survival)
      X <- log2(all_mrna_survival[common_patients,] + 1)

      Y <- all_clinical_survival[common_patients, c("vital_status", "days_to_death", "days_to_last_followup", "days_to_last_known_alive")]
      y <- array(1, dim = nrow(Y))
      y[which(Y$vital_status == "Dead" & Y$days_to_death < 365 * 2)] <- -1
      y <- as.numeric(y)
    }
  }

  if(problem == "sc") {
    if(file.exists(sprintf("%s/xgboost_pathway_%s_%s_measure_AUROC_replication_%d_result.RData", result_path, problem, pathway, replication)) == FALSE) {
      load(sprintf("%s/melanoma_immunotherapy_data.RData", data_path))
      load(sprintf("%s/all_labels.RData", data_path))

      common_patients <- intersect(rownames(all_labels), rownames(melanoma_immunotherapy_data))
      X <- melanoma_immunotherapy_data[common_patients,]
      y <- all_labels$V2
      y <- as.numeric(y)
    }
  }

  class(X) <- "numeric"
  negative_indices <- which(y == -1)
  positive_indices <- which(y == +1)

  fold_count <- 4
  train_ratio <- 0.8
  max_depth_set <- c(1:3) * 2

  if(pathway != "none") {
    pathways <- read_pathways(pathway)
    gene_names <- sort(unique(unlist(sapply(1:length(pathways), FUN = function(x) {pathways[[x]]$symbols}))))
    X <- X[, which(colnames(X) %in% gene_names)]
  }

  set.seed(1606 * replication)
  train_negative_indices <- sample(negative_indices, ceiling(train_ratio * length(negative_indices)))
  train_positive_indices <- sample(positive_indices, ceiling(train_ratio * length(positive_indices)))

  auroc_matrix <- matrix(NA, nrow = fold_count, ncol = length(max_depth_set), dimnames = list(1:fold_count, sprintf("%g", max_depth_set)))

  negative_allocation <- sample(rep(1:fold_count, ceiling(length(train_negative_indices) / fold_count)), length(train_negative_indices))
  positive_allocation <- sample(rep(1:fold_count, ceiling(length(train_positive_indices) / fold_count)), length(train_positive_indices))

  for(fold in 1:fold_count) {
    train_indices <- c(train_negative_indices[which(negative_allocation != fold)], train_positive_indices[which(positive_allocation != fold)])
    test_indices <- c(train_negative_indices[which(negative_allocation == fold)], train_positive_indices[which(positive_allocation == fold)])
    X_train <- X[train_indices,]
    X_test <- X[test_indices,]
    valid_features <- as.numeric(which(apply(X_train, 2, sd) != 0))
    X_train <- X_train[, valid_features]
    X_test <- X_test[, valid_features]
    X_train <- scale(X_train)
    X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)

    y_train <- y[train_indices]
    y_test <- y[test_indices]

    for(max_depth in max_depth_set){
      print(sprintf("running fold = %d, max_depth = %g", fold, max_depth))
      model <- xgboost(data = X_train, label = y_train / 2 + 0.5, objective = "binary:logistic", eta = .2, nrounds = 1000, max_depth = max_depth, early_stopping_rounds = 15)
      y_predicted <- predict(model, X_test)
      auroc_matrix[fold, sprintf("%g", max_depth)] <- auc(roc(predictions = y_predicted, labels = as.factor(y_test / 2 + 0.5)))
    }
  }

  max_depth_star_AUROC <- max_depth_set[max.col(t(colMeans(auroc_matrix, na.rm = TRUE)), ties.method = "last")]

  train_indices <- c(train_negative_indices, train_positive_indices)
  test_indices <- setdiff(1:length(y), train_indices)
  X_train <- X[train_indices,]
  X_test <- X[test_indices,]
  valid_features <- as.numeric(which(apply(X_train, 2, sd) != 0))
  X_train <- X_train[, valid_features]
  X_test <- X_test[, valid_features]
  X_train <- scale(X_train)
  X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)

  y_train <- y[train_indices]
  y_test <- y[test_indices]

  model <- xgboost(data = X_train, label = y_train / 2 + 0.5, objective = "binary:logistic", eta = .2, nrounds = 1000, max_depth = max_depth_star_AUROC, early_stopping_rounds = 15)
  y_predicted <- predict(model, X_test)

  result <- list()
  result$y_predicted <- y_predicted
  result$AUROC <- auc(roc(predictions = y_predicted, labels = as.factor(y_test / 2 + 0.5)))
  result$max_depth_star_AUROC <- max_depth_star_AUROC

  save("result", file = sprintf("%s/xgboost_pathway_%s_%s_measure_AUROC_replication_%d_result.RData", result_path, problem, pathway, replication))
}
