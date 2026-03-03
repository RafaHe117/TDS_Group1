library(dplyr)
library(glmnet)
library(pROC)
library(caret)
library(ggplot2)
library(Matrix)

# Paths + load
setwd("/rds/general/project/hda_25-26/live/TDS/TDS_Group1")
out_dir <- "modelling_script/outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ukb_full <- readRDS("ukb_G1_imputed_final.rds")

# Cohort
if ("cvd_prevalent" %in% names(ukb_full))
  ukb_full <- ukb_full %>% filter(cvd_prevalent == 0)

# Outcome
outcome <- "cvd_incident"

# Drop
drop_cols <- c(
  "eid",
  "cvd_prevalent",
  "cvd_event",
  "cvd_first_date",
  "dod",
  "age_of_death"
)

preds <- setdiff(names(ukb_full), c(outcome, drop_cols))

# Split
set.seed(123)
train_id <- createDataPartition(ukb_full[[outcome]], p = 0.7, list = FALSE)
ukb_train <- ukb_full[train_id, ]
ukb_test  <- ukb_full[-train_id, ]

y_tr <- as.integer(ukb_train[[outcome]])
y_te <- as.integer(ukb_test[[outcome]])

# Sparse design matrices (fast + memory-safe)
x_tr_df <- ukb_train[, preds, drop = FALSE]
x_te_df <- ukb_test[,  preds, drop = FALSE]

char_cols <- names(x_tr_df)[vapply(x_tr_df, is.character, logical(1))]
for (cc in char_cols) {
  x_tr_df[[cc]] <- factor(x_tr_df[[cc]])
  x_te_df[[cc]] <- factor(x_te_df[[cc]], levels = levels(x_tr_df[[cc]]))
}

X_tr <- sparse.model.matrix(~ . - 1, data = x_tr_df)
X_te <- sparse.model.matrix(~ . - 1, data = x_te_df)

# Align columns
miss <- setdiff(colnames(X_tr), colnames(X_te))
if (length(miss) > 0) {
  X_te <- cbind(
    X_te,
    Matrix(0, nrow(X_te), length(miss), sparse = TRUE,
           dimnames = list(NULL, miss))
  )
}
X_te <- X_te[, colnames(X_tr), drop = FALSE]

# =========================================================
# 5) CV LASSO logistic
# =========================================================
set.seed(123)
cvfit <- cv.glmnet(
  x = X_tr,
  y = y_tr,
  family = "binomial",
  alpha = 1,
  nfolds = 10,
  type.measure = "auc",
  standardize = TRUE
)

nz <- cvfit$nzero

lambda_use <- cvfit$lambda.1se   # or cvfit$lambda.min

# =========================================================
# 6) Test AUC + ROC plot
# =========================================================
p_te <- as.numeric(predict(cvfit, newx = X_te, s = lambda_use, type = "response"))

roc_obj <- pROC::roc(y_te, p_te, quiet = TRUE)
auc_te  <- as.numeric(pROC::auc(roc_obj))

roc_df <- data.frame(
  tpr = roc_obj$sensitivities,
  fpr = 1 - roc_obj$specificities
)

p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_line(size = 1.2, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  labs(
    title = paste0("ROC - LASSO (AUC = ", round(auc_te, 3), ")"),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_minimal(base_size = 14)

print(p_roc)

ggsave(
  filename = file.path(out_dir, "roc_lasso.png"),
  plot = p_roc, width = 7, height = 6, dpi = 300
)

writeLines(
  paste0("Test AUC (", outcome, ") = ", round(auc_te, 4),
         "\nlambda_use = ", format(lambda_use, scientific = TRUE),
         "\nnonzero = ", tail(nz, 1)),
  con = file.path(out_dir, "auc_lasso.txt")
)

# =========================================================
# 7) Selected variables
# =========================================================
coef_mat <- as.matrix(coef(cvfit, s = lambda_use))

coef_df <- data.frame(
  term = rownames(coef_mat),
  beta = as.numeric(coef_mat[, 1]),
  row.names = NULL
) %>%
  filter(term != "(Intercept)", beta != 0) %>%
  mutate(abs_beta = abs(beta)) %>%
  arrange(desc(abs_beta))

print(coef_df)

write.csv(coef_df, file.path(out_dir, "selected_variables_lasso.csv"), row.names = FALSE)
saveRDS(cvfit, file.path(out_dir, "cvfit_lasso.rds"))