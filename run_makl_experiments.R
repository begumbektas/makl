library(makl)

data_path <- "./data"
result_path <- "./results"
pathway <- "hallmark" #replace with "pid" if you would like to use PID pathways
problem <- "stage" #replace with "survival" or with "sc" if you would like to perform the experiments on these problems.

read_pathways <- function(name) {
  symbols_lines <- read.table(sprintf("msigdb/%s.gmt", name), header = FALSE, sep = ",", stringsAsFactor = FALSE)
  pathways <- vector("list", nrow(symbols_lines))
  for(line in 1:nrow(symbols_lines)) {
    symbols_entries <- strsplit(symbols_lines[line, 1], "\t")
    names(pathways)[line] <- symbols_entries[[1]][1]
    pathways[[line]]<- sort(symbols_entries[[1]][-2:-1])
  }
  return(pathways)
}

for(replication in 1:100) {
  if(problem == "stage") {
    if(file.exists(sprintf("%s/makl_%s_%s_measure_AUROC_replication_%d_result.RData", result_path, problem, pathway, replication)) == FALSE) {
      load(sprintf("%s/all_clinical_stage.RData", data_path))
      load(sprintf("%s/all_mrna_stage.RData", data_path))

      common_patients <- intersect(rownames(all_clinical_stage_e1)[which(is.na(all_clinical_stage_e1$pathologic_stage) == FALSE & all_clinical_stage_e1$pathologic_stage != "Stage X" & all_clinical_stage_e1$pathologic_stage != "IS")], rownames(all_mrna_stage_e1))

      X <- log2(all_mrna_stage_e1[common_patients,] + 1)
      y <- rep(NA, length(common_patients))

      y[all_clinical_stage_e1[common_patients, "pathologic_stage"] %in% c("Stage I",  "Stage IA",  "Stage IB",  "Stage IC")] <- +1
      y[all_clinical_stage_e1[common_patients, "pathologic_stage"] %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC",
                                                                          "Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC",
                                                                          "Stage IV",  "Stage IVA",  "Stage IVB",  "Stage IVC")] <- -1
    }
  }

  if(problem == "survival") {
    if(file.exists(sprintf("%s/makl_%s_%s_measure_AUROC_replication_%d_result.RData", result_path, problem, pathway, replication)) == FALSE) {
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
    if(file.exists(sprintf("%s/makl_%s_%s_measure_AUROC_replication_%d_result.RData", result_path, problem, pathway, replication)) == FALSE) {
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
  train_ratio <- 0.8
  set.seed(1606 * replication)
  train_negative_indices <- sample(negative_indices, ceiling(train_ratio * length(negative_indices)))
  train_positive_indices <- sample(positive_indices, ceiling(train_ratio * length(positive_indices)))
  train_indices <- c(train_negative_indices, train_positive_indices)
  test_indices <- setdiff(1:length(y), train_indices)

  X_train <- X[train_indices,]
  y_train <- y[train_indices]
  X_test <- X[test_indices,]
  y_test <- y[test_indices]
  pathways <- read_pathways(pathway)

  makl_model <- makl_train(X = X_train, y = y_train, D = 100, sigma_N = 1000, CV = TRUE,
                           lambda_set = c(0.9, 0.8, 0.7, 0.6), membership = pathways)
  result <- makl_test(X = X_test, y = y_test,makl_model = makl_model)
  save("result", file = sprintf("%s/makl_%s_%s_measure_AUROC_replication_%d_result.RData", result_path, problem, pathway, replication))
}
