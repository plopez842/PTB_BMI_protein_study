# Functions are modified from Caroline Loos systems serology Package
### modifed loadings_bar function to have the option to rank by LV or VIP score


# VIP Plot code ---
VIP_plot<-function(model_sel,X,y,my_colors){ #put in model, X, and y
  df_vip <- as.data.frame(getVipVn(model_sel))
  colnames(df_vip)[1] <- "VIP"
  df_vip$features <- rownames(df_vip)
  df_vip$enriched <- ""
  temp<-subset(t_tests_v3,AptName %in% df_vip$features)
  df_vip$genes<-temp[match(df_vip$features,temp$AptName),]$EntrezGeneSymbol
  
  for (ind_feat in 1:nrow(df_vip)) {
    tmp_mean <- rep(NA, length = nlevels(y))
    for (ind_class in 1:nlevels(y)) {
      tmp_mean[ind_class] <- mean(X[which(y == levels(y)[ind_class]),
                                    which(colnames(X) == df_vip$features[ind_feat])])
    }
    df_vip$enriched[ind_feat] <- levels(y)[which.max(tmp_mean)]
  }
  df_vip$enriched <- factor(df_vip$enriched, levels = levels(y))
  df_vip<-subset(df_vip,VIP >1.0)
  
  #reduce to only 10
  
  
  df_vip <- df_vip[(order(-df_vip$VIP)),]
  df_vip<-df_vip[1:15,]
  df_vip$features <- factor(df_vip$features, levels = (df_vip$features))
  df_vip$features <- factor(df_vip$features, labels=df_vip$genes)
  
  plt_bar <- ggplot(data = df_vip, aes(x = features, y = VIP, fill = enriched)) +
    scale_fill_manual(values = my_colors[[1]]) +
    geom_bar(stat = "identity", color = "black") +
    #coord_flip() +
    xlab("") + ylab("VIP score") +
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 15),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = grid::unit(rep(5,4), "mm"),
          axis.title = element_text(size = 16,face="bold"),
          legend.position = "none")#, angle = 45, hjust = 1))
  
  
  return(plt_bar)
}

# Visualize Bar Loadings ----
visualize_ropls_loadings_bar <- function(model, options = list(),order_by="VIP") {
  
  # ----------------- BEGIN OPTIONS ----------------- #
  if (!("LV_ind" %in% names(options))) {
    options$LV_ind <- c(1)
  }
  if ("y" %in% names(options)) {
    y <- options$y
    if (is.factor(y)) {
      n_groups <- nlevels(y)
    } else {
      n_groups <- NA
    }
  } else {
    n_groups <- NA
  }
  if (!("mark_enrichment" %in% names(options)) | is.na(n_groups)) {
    options$mark_enrichment <- FALSE
  }
  if (options$mark_enrichment & (is.na(n_groups) | !("X" %in% names(options))))  {
    stop("Enrichment only works for classification and when X and y are provided")
  }
  
  # color for the scores and name of the grouping
  if (!("y_name" %in% names(options))) {
    y_name <- "y"
  } else {
    y_name <- options$y_name
  }
  # color for the scores and name of the grouping
  if (!("colors" %in% names(options)) | length(grep(y_name, names(options$colors))) == 0) {
    if (is.factor(y)) {
      tmp <- rep(NA, length = nlevels(y))
      names(tmp) <- levels(y)
      for (ind in 1:nlevels(y)) {
        tmp[ind] <- RColorBrewer::brewer.pal(n = max(3, nlevels(y)), name = 'Dark2')[ind]
      }
      options$colors <- list()
      options$colors[[y_name]] <- tmp
    } else {
      # For regression, a color palette needs to be provided
      options$colors$y <- list(low = "#C7E4F9", high = "#004D7F")
    }
  }
  
  if (ropls::getSummaryDF(model)$pre +
      ropls::getSummaryDF(model)$ort < options$LV_ind) {
    stop("required LV exceed existing LVs")
  }
  
  if (!("df_features" %in% names(options))) {
    options$df_features <- data.frame(name = rownames(model@loadingMN),
                                      label = rownames(model@loadingMN))
  }
  
  # ----------------- END OPTIONS ----------------- #
  
  # check first whether its a orthogonal PLS or a regular PLS
  if (ropls::getSummaryDF(model)$ort > 0) {
    stop("orthogonal PLS-DA not supported yet")
    # if (options$LV_ind[1] == 1) {
    #   df_loadings <- data.frame(LV1 = ropls::getLoadingMN(model),
    #                             LV2 = ropls::getLoadingMN(model, orthoL = TRUE)[,options$LV_ind[2] - 1])
    # } else {
    #   df_loadings <- data.frame(LV1 = ropls::getLoadingMN(model, orthoL = TRUE)[, options$LV_ind[1] - 1],
    #                             LV2 = ropls::getLoadingMN(model, orthoL = TRUE)[, options$LV_ind[2] - 1])
    # }
  } else {
    df_loadings <- data.frame(LV = ropls::getLoadingMN(model)[,options$LV_ind[1]],
                              vip_scores = ropls::getVipVn(model))
    df_loadings$features <- rownames(df_loadings)
    df_loadings$labels <- options$df_features$label[match(rownames(df_loadings), options$df_features$name)]
  }
  
  
  
  # TODO: catch if its an orthogonal
  
  if (options$mark_enrichment & !is.na(n_groups)) {
    df_loadings$mark <- NA
    X <- options$X
    
    for (ind_feat in 1:nrow(df_loadings)) {
      tmp_mean <- rep(NA, length = nlevels(y))
      for (ind_class in 1:nlevels(y)) {
        tmp_mean[ind_class] <- mean(X[which(y == levels(y)[ind_class]),
                                      which(colnames(X) == df_loadings$features[ind_feat])])
      }
      df_loadings$mark[ind_feat] <- levels(y)[which.max(tmp_mean)]
    }
    df_loadings$mark  <- factor(df_loadings$mark, levels = levels(y))
  }
  
  if(order_by =="VIP"){
    df_loadings <- df_loadings[order(df_loadings$vip_scores), ]
  }
  else(
    df_loadings <- df_loadings[order(df_loadings$LV,decreasing=FALSE), ]
  )
 # df_loadings <- df_loadings[order(df_loadings$vip_scores), ]
  df_loadings$features <- factor(df_loadings$features, levels = unique(df_loadings$features))
  
  # plot loadings sorted according to the VIP score and color coding it
  # according to enrichent in classes
  if (options$mark_enrichment) {
    plt_bar <- ggplot2::ggplot(data = df_loadings, ggplot2::aes(x = features, y = LV, fill = mark)) +
      ggplot2::scale_fill_manual(values = options$colors[[y_name]])
  } else {
    plt_bar <- ggplot2::ggplot(data = df_loadings, ggplot2::aes(x = features, y = LV))
  }
  plt_bar <- plt_bar +
    ggplot2::geom_bar(stat = "identity", color = "black") +
    ggplot2::coord_flip() +
    ggplot2::xlab("") +
    ggplot2::ylab(paste("LV", options$LV_ind, " loadings", sep = "")) +
    ggplot2::labs(fill = "enriched in") +
    ggplot2::scale_x_discrete(labels = df_loadings$labels) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(color = "black"))#,
  #axis.text.y = element_text(colour = as.character(feature_annot$useColor[match(dfBar$features[order(dfBar$vipScores)],
  # rownames(feature_annot))])))
  
}


# Score Accuracy Func ----

score_accuracy <- function(y, y_pred) {
  num <- as.numeric(y)
  num_pred <- as.numeric(y_pred)
  correct <- which(num == num_pred)
  accuracy <- length(correct) / length(y)
  return(accuracy)
}

# Training ROPLS model ----

train_ropls <- function(X, y, options = list()) {
  # suppress annoying "error"s from ropls
  #sink(file = tempfile())
  if ("n_LV" %in% names(options)) {
    predI <- options$n_LV
  } else {
    predI <- NA
  }
  try_out <- try(model <- ropls::opls(X, y,
                                      permI = 0, # no permutation and other output to save computation time
                                      predI = predI,
                                      scaleC = 'none',
                                      info.txtC = "none",
                                      fig.pdfC = "none"
  )
  )
  if (is(try_out, "try-error")) {
    # No model could be build, provide a model with one latent variable
    # (even if this is not significant)
    model <- ropls::opls(X, y, predI = 1,
                         permI = 0,
                         scaleC = 'none',
                         info.txtC = "none",
                         fig.pdfC = "none"
    )
  }
  return(model)
}

# repeating regulariztion ----
select_repeat <- function(X, y, selector, options = list()) {
  
  # ----------------- BEGIN OPTIONS ----------------- #
  # How often it should be repeated
  if (!("n_trials" %in% names(options))) {
    options$n_trials <- 100
  }
  if (!("threshold" %in% names(options))) {
    options$threshold <- 0.8
  }
  # returns the whole data frame of how often a feature was selected
  if (!("return_count" %in% names(options))) {
    options$return_count <- FALSE
  }
  if (!("force_select" %in% names(options))) {
    options$force_select <- TRUE
  }
  
  if (!("alpha" %in% names(options))) {
    options$alpha <- 1
  }
  # ----------------- END OPTIONS ----------------- #
  
  # vector counting how often each feature is selected
  feature_count <- rep(0, ncol(X))
  names(feature_count) <- colnames(X)
  
  # run the feature selector trials times and increment the counters
  # for the features that are selected
  for (trial in 1:options$n_trials) {
    features <- selector(X, y, options=list(alpha = options$alpha))
    feature_count[features] <- 1 + feature_count[features]
    #print(feature_count)
  }
  
  # keep those features that were selected in more than threshold
  # percent of the trials
  selected <- feature_count[-which(feature_count <= options$threshold * options$n_trials)]
  
  # if a selection is forced, return the features which are selected the most
  if (length(selected) == 0 & options$force_select) {
    selected <- feature_count[which(feature_count == max(feature_count))]
  }
  
  if (options$return_count) {
    return(list(feature_count = feature_count, sel_features = names(selected)))
  } else {
    return(names(selected))
  }
}

# Select Lasso ----
select_lasso <- function(X, y, options = list()) {
  # decide on type of GLM depending on type of y
  # 2-level factor -> binomial, n-level factor -> multinomial
  # vector -> gaussian, anything else gives an error
  if (is.factor(y)) {
    if (nlevels(y) == 1) {
      stop("y is a factor with only one level")
    } else if (nlevels(y) == 2) {
      fam <- "binomial"
    } else {
      fam <- "multinomial"
    }
  } else if (is.numeric(y)) {
    fam <- "gaussian"
  } else {
    stop("y must be of type factor or numeric vector")
  }
  
  # default alpha to 1 if it is not set
  if (!("alpha" %in% names(options))) {
    options$alpha <- 1
  }
  
  # cv.glmnet needs at least 3 folds, so we need at least three features
  n_samples <- nrow(X)
  if (n_samples < 3) {
    stop("select_lasso() requires more than three samples for internal cross-validation")
  }
  
  # default subfolds if it is not set
  # if there are 5 or more samples, default to 5
  # otherwise default to 3
  # while were at it, ensure that subfolds is in [3, n_samples]
  # and print a warning if it wasn't
  if (!("subfolds" %in% names(options))) {
    if (n_samples >= 5) {
      options$subfolds <- 5
    } else {
      options$subfolds <- 3
    }
  } else {
    if (options$subfolds > n_samples) {
      message("Warning in select_lasso():")
      message("    options$subfolds greater than number of samples")
      message("    setting options$subfolds = number of samples")
      options$subfolds <- n_samples
    }
    if (options$subfolds < 3) {
      message("Warning in select_lasso():")
      message("    options$subfolds was less than 3")
      message("    setting options$subfolds = 3")
      options$subfolds <- 3
    }
  }
  
  # fit an appropriate lasso model with a number of trials corresponding to
  # different values of lambda
  lasso <- glmnet::cv.glmnet(X, y, type.measure = "mse", alpha = options$alpha,
                             family = fam, type.multinomial = "grouped",
                             nfolds = options$subfolds)
  
  # lasso$lambda[k] is the value of lambda in the k-th trial
  # lasso$nzero[k] is the number of non-zero coefficients in the fitted model
  # (= number of features not including the intercept) in the k-th trial
  # things are arranged such that the first entry of nzero is 0
  # and the final entry of nzero equals the number of total features
  # lasso$cvm is the cross-validated score of the k-th trial
  
  # find the model with the smallest error that has at least one non-zero
  # coefficient other than the intercept
  indices <- which(lasso$nzero > 0)
  lambdas <- lasso$lambda[indices]
  scores <- lasso$cvm[indices]
  # if there is more than one index attaining the minimum which.min picks
  # the smallest one - this corresponds to choosing best score with the
  # minimal number of features selected
  best <- which.min(scores)
  lasso_coeffs <- coef(lasso, s = lambdas[best])
  
  if (fam == "multinomial") {
    # if the data has multiple responses, the coefficients are a matrix
    # that is returned as a list of columns. type.multinomial = "grouped"
    # forced features to be selected for all responses or for none, so we
    # can get the selected features also by only considering the first
    # column. we just replace lasso_coeffs by this column and proceed as usual
    lasso_coeffs <- lasso_coeffs[[1]]
  }
  
  # remove the intercept and the entries that are zero
  lasso_coeffs <- lasso_coeffs[-1,]
  lasso_coeffs <- lasso_coeffs[which(lasso_coeffs != 0)]
  
  # return the names of the selected features. previous code turned
  # coefficients into a vector, so use names() rather than rownames()
  return(names(lasso_coeffs))
}
   

# Visualize PLSDA biplot ----


visualize_ropls_scores <- function(model, y, options = list()) {
  
  # ----------------- OPTIONS ----------------- #
  if (!("alpha" %in% names(options))) {
    options$alpha <- 1
  }
  if (!("size" %in% names(options))) {
    options$size <- 2.5
  }
  if (!("stroke" %in% names(options))) {
    options$stroke <- 0.5
  }
  # level at which to draw the ellipse for the data ellipse
  if (!("level" %in% names(options))) {
    options$level <- 0.95
  }
  if ("y" %in% names(options)) {
    y <- options$y
    if (is.factor(y)) {
      n_groups <- nlevels(y)
    } else {
      n_groups <- NA
    }
  } else {
    n_groups <- NA
  }
  if (!("y_name" %in% names(options))) {
    y_name <- "y"
  } else {
    y_name <- options$y_name
  }
  # color for the scores and name of the grouping
  if (!("colors" %in% names(options)) | length(grep(y_name, names(options$colors))) == 0) {
    if (is.factor(y)) {
      tmp <- rep(NA, length = nlevels(y))
      names(tmp) <- levels(y)
      for (ind in 1:nlevels(y)) {
        tmp[ind] <- RColorBrewer::brewer.pal(n = max(3, nlevels(y)), name = 'Dark2')[ind]
      }
      options$colors <- list()
      options$colors[[y_name]] <- tmp
    } else {
      # For regression, a color palette needs to be provided
      options$colors$y <- list(low = "#C7E4F9", high = "#004D7F")
    }
  }
  
  
  # which latent variables to check, defaults to the first two
  # if they are provided, ensure that the model has the required number of LVs
  if (!("LV_ind" %in% names(options))) {
    options$LV_ind <- c(1,2)
  } else if (ropls::getSummaryDF(model)$pre +
             ropls::getSummaryDF(model)$ort < max(options$LV_ind)) {
    stop("required LV exceed existing LVs")
  } else if (!(length(options$LV_ind) == 2)) {
    stop("two LVs required")
  }
  # ----------------- END OPTIONS ----------------- #
  
  # ----------------- GET SCORES ----------------- #
  # check first whether its a orthogonal PLS or a regular PLS
  if (ropls::getSummaryDF(model)$ort > 0) {
    if (options$LV_ind[1] == 1) {
      df_scores <- data.frame(LV1 = ropls::getScoreMN(model),
                              LV2 = ropls::getScoreMN(model, orthoL = TRUE)[, options$LV_ind[2] - 1],
                              y = y)
    } else {
      df_scores <- data.frame(LV1 = ropls::getScoreMN(model, orthoL = TRUE)[, options$LV_ind[1] - 1],
                              LV2 = ropls::getScoreMN(model, orthoL = TRUE)[, options$LV_ind[2] - 1],
                              y = y)
    }
  } else {
    df_scores <- data.frame(LV1 = ropls::getScoreMN(model)[,options$LV_ind[1]],
                            LV2 = ropls::getScoreMN(model)[,options$LV_ind[2]],
                            y = y)
  }
  # --------------------- END GET SCORES ------------------- #
  
  
  # ---------------------- BEGIN PLOT ---------------------- #
  plt_scores <- ggplot2::ggplot(df_scores, ggplot2::aes(LV1, LV2, fill = y)) +
    ggplot2::geom_vline(xintercept = 0, size = 0.3) +
    ggplot2::geom_hline(yintercept = 0, size = 0.3) +
    ggplot2::geom_point(color = "black",
                        size = options$size,
                        alpha = options$alpha,
                        stroke = options$stroke,
                        shape = 21,
                        show.legend = TRUE) +
    ggplot2::labs(x = paste("scores on LV", options$LV_ind[1], " (",
                            toString(round(model@modelDF$R2X[options$LV_ind[1]] * 100)), "%)", sep = ""),
                  y = paste("scores on LV", options$LV_ind[2], " (",
                            toString(round(model@modelDF$R2X[options$LV_ind[2]] * 100)), "%)", sep = ""),
                  fill = y_name, color = y_name) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right",
                   aspect.ratio = 1,
                   axis.text = ggplot2::element_text(color = "black"))
  
  # in the case of a classification, add ellipses
  if (is.factor(y)) {
    plt_scores <- plt_scores +
      ggplot2::stat_ellipse(ggplot2::aes(color = y), level = options$level) +
      ggplot2::scale_fill_manual(values = options$colors[[y_name]]) +
      ggplot2::scale_color_manual(values = options$colors[[y_name]])
  } else {
    plt_scores <- plt_scores +
      ggplot2::scale_fill_gradient(low = options$colors$y[["low"]], high = options$colors$y[["high"]])
  }
  # ---------------------- END PLOT ---------------------- #
  
  return(plt_scores)
}


train_ropls <- function(X, y, options = list()) {
  # suppress annoying "error"s from ropls
  #sink(file = tempfile())
  if ("n_LV" %in% names(options)) {
    predI <- options$n_LV
  } else {
    predI <- NA
  }
  try_out <- try(model <- ropls::opls(X, y,
                                      permI = 0, # no permutation and other output to save computation time
                                      predI = predI,
                                      scaleC = 'none',
                                      info.txtC = "none",
                                      fig.pdfC = "none"
  )
  )
  if (is(try_out, "try-error")) {
    # No model could be build, provide a model with one latent variable
    # (even if this is not significant)
    model <- ropls::opls(X, y, predI = 1,
                         permI = 0,
                         scaleC = 'none',
                         info.txtC = "none",
                         fig.pdfC = "none"
    )
  }
  return(model)
}

validate_repeat <- function(X, y, method, options, n_trials = 100) {
  # we run validate() n_trials times, returning the vector of real scores
  # and the matrices of random feature scores/permutation test scores
  
  # ----------------- OPTIONS ----------------- #
  if (!("save" %in% names(options))) {
    options$save <- FALSE
  }
  if (options$save & !("save_file" %in% names(options))) {
    options$save_file <- paste0("model_validation_", Sys.time())
  }
  # ----------------- END OPTIONS ----------------- #
  
  
  return_vals <- list()
  
  for (trial in 1:n_trials) {
    message(paste("validate_repeat: trial", trial, "/", n_trials))
    val_scores <- validate(X, y, method, options)
    # remove the actual prediction from the validation
    val_scores["cv_y"] <- NULL
    return_vals[[trial]] <- val_scores
    if (options$save) {
      saveRDS(return_vals, file = paste0(options$save_file, ".RDS"))
    }
  }
  return(return_vals)
}


validate <- function(X, y, method, options) {
  # ----------------- INITIAL PROCESSING ----------------- #
  
  # see if a feature selector is passed, otherwise default
  # to selecting all features. in the latter case, also make
  # sure rf_trials = 0 since random features don't make sense
  if ("select" %in% names(method) ) {
    select <- method$select
  } else {
    select <- function(X,y) {return(colnames(X))}
    if (!("rf_trials" %in% names(options))) {
      options$rf_trials <- 0
    } else if (options$rf_trials != 0) {
      message("Warning in validate():")
      message("    no feature selector given but rf_trials != 0")
      message("    forcing rf_trials = 0")
      options$rf_trials <- 0
    }
  }
  
  # default to five-fold cross-validation
  if (!("n_folds" %in% names(options))) {
    options$n_folds <- 5
  }
  
  # for permuted labels, to which the predicted outcome should be compared to with the score function
  if (!("compare_pred" %in% names(options))) {
    options$compare_pred <- "y"
  } else {
    if (!(options$compare_pred %in% c("y", "y_perm"))) {
      stop("options$compare_pred needs to be \"y\" or \"y_perm\"")
    }
  }
  
  # for paired data, take structure into account
  if (!("paired" %in% names(options))) {
    options$paired <- FALSE
  }
  if (options$paired) {
    if (!("X_label" %in% names(options))) {
      stop("X_label needs to be provided to take into account the paired structure in cross-validation and permutation testing.")
    }
  }
  
  # also give these guys some shorter names
  train <- method$train
  predict <- method$predict
  score <- method$score
  
  # add return values to this list as we go along
  return_values <- list()
  
  # stores the number of features selected for each fold
  # during cross-validation. we need this later for the
  # random features test
  feats_per_fold <- list()
  
  # if score is a single function instead of a list of
  # functions, wrap it in a list
  if (class(score) != "list") {
    score <- list(score)
  }
  # ----------------- END INITIAL PROCESSING ----------------- #
  
  
  
  # ----------------- BEGIN CROSS-VALIDATION ----------------- #
  # split data into folds
  if (options$paired) {
    folds <- caret::createFolds(seq(1, nrow(X)/2), options$n_folds)
    fold_names <- names(folds)
    
    for (fname in fold_names) {
      folds[[fname]] <- which(options$X_label %in% options$X_label[folds[[fname]]])
    }
    
  } else {
    folds <- caret::createFolds(y, options$n_folds)
    fold_names <- names(folds)
  }
  
  # vector of cross-validation predictions
  y_pred <- y
  
  for (fname in fold_names) {
    indices <- folds[[fname]]
    X_train <- X[-indices, , drop = FALSE]
    y_train <- y[-indices]
    X_pred <- X[indices, , drop = FALSE]
    
    features <- select(X_train, y_train)
    
    # actually, check more for valid indices...
    if (length(features) == 0) {
      stop("method$select() did not return any features")
    }
    
    # store number of features selected in fold for later
    feats_per_fold[[fname]] <- length(features)
    
    model <- caret::train(as.matrix(X_train[, features, drop = FALSE]), y_train)
    
    y_pred[indices] <- predict(model, as.matrix(X_pred[, features, drop = FALSE]))
  }
  
  return_values$cv_y <- y_pred
  
  # apply the list of functions in score to y_pred, y
  f_star <- function(f) {f(y, y_pred)}
  return_values$cv_score <- lapply(score, f_star)
  # ----------------- END CROSS-VALIDATION ----------------- #
  
  
  
  # ----------------- BEGIN RANDOM FEATURES ----------------- #
  n_trials <- options$rf_trials
  if (n_trials > 0) {
    n_scores <- length(score)
    rf_scores <- list(vector(mode = "numeric", length = n_trials))
    rf_scores <- rep(rf_scores, n_scores)
    
    for (trial in 1:n_trials) {
      
      for (fname in fold_names) {
        indices <- folds[[fname]]
        X_train <- X[-indices, , drop = FALSE]
        y_train <- y[-indices]
        X_pred <- X[indices, , drop = FALSE]
        
        # careful with sample() pathology here...
        # select random features
        total_features <- ncol(X_train)
        features <- sample(1:total_features, feats_per_fold[[fname]])
        
        model <- train(as.matrix(X_train[, features, drop = FALSE]), y_train)
        
        y_pred[indices] <- predict(model, as.matrix(X_pred[, features, drop = FALSE]))
      }
      
      # compute list of scores
      score_list <- lapply(score, f_star)
      
      # assign them to vectors in the list
      rf_scores <- lv_assign(rf_scores, score_list, trial)
    }
    
    return_values$rf_scores <- rf_scores
  }
  # ----------------- END RANDOM FEATURES ----------------- #
  
  
  
  # ----------------- BEGIN PERMUTATION TESTING ----------------- #
  n_trials <- options$pt_trials
  
  if (n_trials > 0) {
    n_scores <- length(score)
    pt_scores <- list(vector(mode = "numeric", length = n_trials))
    pt_scores <- rep(pt_scores, n_scores)
    
    y_perm <- y
    
    for (trial in 1:n_trials) {
      if (options$paired) { # create permuted y, flip pairs with a 50% probability
        for (fname in fold_names) {
          indices <- folds[[fname]]
          tmp_rn <- runif(length(indices)/2)
          flip_pairs <- which(tmp_rn > 0.5)
          flip_pairs_labels <- unique(options$X_label[indices])[flip_pairs]
          
          y_perm[indices] <- y[indices]
          
          for (ind_pair in flip_pairs_labels) {
            y_perm[which(options$X_label == ind_pair)[2]] <- y[which(options$X_label == ind_pair)[1]]
            y_perm[which(options$X_label == ind_pair)[1]] <- y[which(options$X_label == ind_pair)[2]]
          }
        }
      } else {# create permuted y, but only permute inside each fold
        for (fname in fold_names) {
          indices <- folds[[fname]]
          perm <- sample(1:length(indices))
          y_perm[indices] <- y[indices[perm]]
        }
      }
      
      for (fname in fold_names) {
        indices <- folds[[fname]]
        X_train <- X[-indices, , drop = FALSE]
        y_train <- y_perm[-indices]
        X_pred <- X[indices, , drop = FALSE]
        
        features <- select(X_train, y_train)
        model <- train(as.matrix(X_train[, features, drop = FALSE]), y_train)
        
        y_pred[indices] <- predict(model, as.matrix(X_pred[, features, drop = FALSE]))
      }
      
      if (options$compare_pred == "y") {
        # compute list of scores
        score_list <- lapply(score, f_star)
      } else if (options$compare_pred == "y_perm") {
        f_star_perm <- function(f) {f(y_perm, y_pred)}
        score_list <- lapply(score, f_star_perm)
      }
      
      pt_scores <- lv_assign(pt_scores, score_list, trial)
    }
    
    return_values$pt_scores <- pt_scores
  }
  # ----------------- END PERMUTATION TESTING ----------------- #
  
  
  
  # ----------------- BEGIN FINAL PROCESSING ----------------- #
  # unpack the list if its length is 1 (score was just a function)
  if (length(score) == 1) {
    return_values$cv_score <- return_values$cv_score[[1]]
    
    if ("rf_scores" %in% names(return_values)) {
      return_values$rf_scores <- return_values$rf_scores[[1]]
    }
    
    if ("pt_scores" %in% names(return_values)) {
      return_values$pt_scores <- return_values$pt_scores[[1]]
    }
  }
  # ----------------- END FINAL PROCESSING ----------------- #
  
  return(return_values)
}


# for each index in vec_list, it sets
# vec_list[[ind]][v_index] = val_list[[ind]]
lv_assign <- function(vec_list, val_list, v_index) {
  list_indices <- c(1:length(vec_list))
  for (ind in list_indices) {
    vec <- vec_list[[ind]]
    vec[v_index] <- val_list[[ind]]
    vec_list[[ind]] <- vec
  }
  return(vec_list)
}

predict_ropls <- function(model, X) {
  y_pred <- ropls::predict(model, newdata = X)
  return(y_pred)
}

#' Cross-validating a PLSDA model
#' This code taken from systemsseRology package (https://github.com/LoosC/systemsseRology)
#' and modified for purpose

cross_validation_unpaired <- function(X, y, method, options, n_trials) {
  X <- as.matrix(X)
  y <- y
  
  vals_list <- list()
  
  if ("select" %in% names(method) ) {
    select <- method$select
  } else {
    select <- function(X,y) {return(colnames(X))}
    if (!("rf_trials" %in% names(options))) {
      options$rf_trials <- 0
    } else if (options$rf_trials != 0) {
      message("Warning in validate():")
      message("    no feature selector given but rf_trials != 0")
      message("    forcing rf_trials = 0")
      options$rf_trials <- 0
    }
  }
  
  # default to five-fold cross-validation
  if (!("n_folds" %in% names(options))) {
    options$n_folds <- 5
  }
  
  # also give these guys some shorter names
  train <- method$train
  predict <- method$predict
  score <- method$score
  
  
  for (trial in 1:n_trials) {
    message(paste("validate_repeat: trial", trial, "/", n_trials))
    
    # ----------------- INITIAL PROCESSING ----------------- #
    
    # see if a feature selector is passed, otherwise default
    # to selecting all features. in the latter case, also make
    # sure rf_trials = 0 since random features don't make sense
    
    
    
    # for permuted labels, to which the predicted outcome should be compared to with the score function
    if (!("compare_pred" %in% names(options))) {
      options$compare_pred <- "y"
    } else {
      if (!(options$compare_pred %in% c("y", "y_perm"))) {
        stop("options$compare_pred needs to be \"y\" or \"y_perm\"")
      }
    }
    
    # for paired data, take structure into account
    if (!("paired" %in% names(options))) {
      options$paired <- FALSE
    }
    if (options$paired) {
      if (!("X_label" %in% names(options))) {
        stop("X_label needs to be provided to take into account the paired structure in cross-validation and permutation testing.")
      }
    }
    
    # add return values to this list as we go along
    return_values <- list()
    
    # stores the number of features selected for each fold
    # during cross-validation. we need this later for the
    # random features test
    feats_per_fold <- list()
    
    # if score is a single function instead of a list of
    # functions, wrap it in a list
    if (class(score) != "list") {
      score <- list(score)
    }
    # ----------------- END INITIAL PROCESSING ----------------- #
    
    
    
    # ----------------- BEGIN CROSS-VALIDATION ----------------- #
    # split data into folds
    if (options$paired) {
      folds <- caret::createFolds(seq(1, nrow(X)/2), options$n_folds)
      fold_names <- names(folds)
      
      for (fname in fold_names) {
        folds[[fname]] <- which(options$X_label %in% options$X_label[folds[[fname]]])
      }
      
    } else {
      folds <- caret::createFolds(y, options$n_folds)
      fold_names <- names(folds)
    }
    
    # vector of cross-validation predictions
    y_pred <- y
    
    for (fname in fold_names) {
      indices <- folds[[fname]]
      X_train <- X[-indices, , drop = FALSE]
      y_train <- y[-indices]
      X_pred <- X[indices, , drop = FALSE]
      
      real_features <- select(X_train, y_train)
      all_features <- colnames(X)
      
      # actually, check more for valid indices...
      if (length(real_features) == 0) {
        stop("method$select() did not return any features")
      }
      
      # store number of features selected in fold for later
      feats_per_fold[[fname]] <- length(real_features)
      
      model <- train(as.matrix(X_train[, real_features, drop = FALSE]), y_train)
      try_out <- try(y_pred[indices] <- predict(model,
                                                as.matrix(X_pred[, real_features, drop = FALSE])),
                     silent = T)
      if (is(try_out, "try-error")) {
        # No model could be build, provide a model with two latent variables
        # (even if this is not significant)
        opts_model <- list(n_LV = 2)
        model <- train(as.matrix(X_train[, real_features, drop = FALSE]), y_train, opts_model)
        y_pred[indices] <- predict(model,
                                   as.matrix(X_pred[, real_features, drop = FALSE]))
      }
    }
    
    return_values$cv_y <- y_pred
    
    # apply the list of functions in score to y_pred, y
    f_star <- function(f) {f(y, y_pred)}
    return_values$cv_score <- lapply(score, f_star)
    print("end CV")
    # ----------------- END CROSS-VALIDATION ----------------- #
    
    
    
    # ----------------- BEGIN RANDOM FEATURES ----------------- #
    n_trials_r <- options$rf_trials
    if (n_trials_r > 0) {
      n_scores <- length(score)
      rf_scores <- list(vector(mode = "numeric", length = n_trials_r))
      rf_scores <- rep(rf_scores, n_scores)
      
      for (trial_r in 1:n_trials_r) {
        
        for (fname in fold_names) {
          indices <- folds[[fname]]
          X_train <- X[-indices, , drop = FALSE]
          y_train <- y[-indices]
          X_pred <- X[indices, , drop = FALSE]
          
          # careful with sample() pathology here...
          # select random features that are NOT ones already chosen
          total_features <- colnames(X)
          nonoverlap_features <- total_features[-which(total_features %in% real_features)]
          ##FEATS PER FOLD IS USING TRIAL_R FNAME AND NOT OVERALL FNAME (wait this is okay)
          random_features <- sample(nonoverlap_features, feats_per_fold[[fname]])
          opts_rf <- list(n_LV = 2)
          model <- train(as.matrix(X_train[, random_features, drop = FALSE]), y_train, opts_rf)
          
          y_pred[indices] <- predict(model, as.matrix(X_pred[, random_features, drop = FALSE]))
        }
        
        # compute list of scores
        score_list <- lapply(score, f_star)
        
        # assign them to vectors in the list
        rf_scores <- lv_assign(rf_scores, score_list, trial_r)
      }
      
      return_values$rf_scores <- rf_scores
    }
    print("end random features")
    # ----------------- END RANDOM FEATURES ----------------- #
    
    
    
    # ----------------- BEGIN PERMUTATION TESTING ----------------- #
    n_trials_p <- options$pt_trials
    
    if (n_trials_p > 0) {
      n_scores <- length(score)
      pt_scores <- list(vector(mode = "numeric", length = n_trials_p))
      pt_scores <- rep(pt_scores, n_scores)
      
      y_perm <- y
      
      for (trial_p in 1:n_trials_p) {
        print(paste0("STARTING_PERM_TRIAL_",trial_p))
        if (options$paired) { # create permuted y, flip pairs with a 50% probability
          for (fname in fold_names) {
            indices <- folds[[fname]]
            tmp_rn <- runif(length(indices)/2)
            flip_pairs <- which(tmp_rn > 0.5)
            flip_pairs_labels <- unique(options$X_label[indices])[flip_pairs]
            
            y_perm[indices] <- y[indices]
            
            for (ind_pair in flip_pairs_labels) {
              y_perm[which(options$X_label == ind_pair)[2]] <- y[which(options$X_label == ind_pair)[1]]
              y_perm[which(options$X_label == ind_pair)[1]] <- y[which(options$X_label == ind_pair)[2]]
            }
          }
        } else {# create permuted y, but only permute inside each fold
          for (fname in fold_names) {
            indices <- folds[[fname]]
            perm <- sample(1:length(indices))
            y_perm[indices] <- y[indices[perm]]
          }
        }
        
        for (fname in fold_names) {
          indices <- folds[[fname]]
          X_train <- X[-indices, , drop = FALSE]
          y_train <- y_perm[-indices]
          X_pred <- X[indices, , drop = FALSE]
          
          features <- select(X_train, y_train)
          model <- train(as.matrix(X_train[, features, drop = FALSE]), y_train)
          
          try_out <- try(y_pred[indices] <- predict(model,
                                                    as.matrix(X_pred[, features, drop = FALSE])),
                         silent = T)
          if (is(try_out, "try-error")) {
            # No model could be build, provide a model with two latent variables
            # (even if this is not significant)
            opts_model <- list(n_LV = 2)
            model <- train(as.matrix(X_train[, features, drop = FALSE]), y_train, opts_model)
            y_pred[indices] <- predict(model,
                                       as.matrix(X_pred[, features, drop = FALSE]))
          }
        }
        
        if (options$compare_pred == "y") {
          # compute list of scores
          score_list <- lapply(score, f_star)
        } else if (options$compare_pred == "y_perm") {
          f_star_perm <- function(f) {f(y_perm, y_pred)}
          score_list <- lapply(score, f_star_perm)
        }
        
        pt_scores <- lv_assign(pt_scores, score_list, trial_p)
        print(paste0("ENDING_PERM_TRIAL_",trial_p))
      }
      
      return_values$pt_scores <- pt_scores
    }
    # ----------------- END PERMUTATION TESTING ----------------- #
    print("end permutation testing")
    
    
    # ----------------- BEGIN FINAL PROCESSING ----------------- #
    # unpack the list if its length is 1 (score was just a function)
    if (length(score) == 1) {
      return_values$cv_score <- return_values$cv_score[[1]]
      
      if ("rf_scores" %in% names(return_values)) {
        return_values$rf_scores <- return_values$rf_scores[[1]]
      }
      
      if ("pt_scores" %in% names(return_values)) {
        return_values$pt_scores <- return_values$pt_scores[[1]]
      }
    }
    # ----------------- END FINAL PROCESSING ----------------- #
    
    # remove the actual prediction from the validation
    return_values["cv_y"] <- NULL
    vals_list[[trial]] <- return_values
  }
  return(vals_list)
}

# for each index in vec_list, it sets
# vec_list[[ind]][v_index] = val_list[[ind]]
lv_assign <- function(vec_list, val_list, v_index) {
  list_indices <- c(1:length(vec_list))
  for (ind in list_indices) {
    vec <- vec_list[[ind]]
    vec[v_index] <- val_list[[ind]]
    vec_list[[ind]] <- vec
  }
  return(vec_list)
}


suppressMessages({
  suppressWarnings({
    library(glmnet, quietly = TRUE)
    library(RhpcBLASctl, quietly = TRUE)
    library(foreach, quietly = TRUE)
    library(doSNOW, quietly = TRUE)
    library(parallelly,quietly=TRUE)
  })
})

feature_select_reg <- function(X, y, options = list(alpha = 0.6, subfolds = 5, error.measure = 'mse'), ...) {
  # select features based on non-zero coefs returns from CV glm
  # documentation: https://glmnet.stanford.edu/reference/glmnet.html, https://glmnet.stanford.edu/reference/cv.glmnet.html
  
  # GLM type
  if (is.factor(y)) {
    if (nlevels(y) == 1) {
      stop("y is a factor with only one level")
    } else if (nlevels(y) == 2) {
      fam <- "binomial"
    } else {
      fam <- "multinomial"
    }
  } else if (is.numeric(y)) {
    fam <- "gaussian"
  } else {
    stop("y must be of type factor or numeric vector")
  }
  
  if (!("alpha" %in% names(options))) {
    options$alpha <- 1 # lasso
  }
  if (!("error.measure" %in% names(options))) {
    options$error.measure <- 'mse' # lasso
  }
  
  # cv.glmnet needs at least 3 folds, so we need at least three features
  n_samples <- nrow(X)
  if (n_samples < 3) {
    stop("select_lasso() requires more than three samples for internal cross-validation")
  }
  
  if (!("subfolds" %in% names(options))) {
    if (n_samples >= 5) {
      options$subfolds <- 5
    } else {
      options$subfolds <- 3
    }
  } else {
    if (options$subfolds > n_samples) {
      message("Warning in select_lasso():")
      message("    options$subfolds greater than number of samples")
      message("    setting options$subfolds = number of samples")
      options$subfolds <- n_samples
    }
    if (options$subfolds < 3) {
      message("Warning in select_lasso():")
      message("    options$subfolds was less than 3")
      message("    setting options$subfolds = 3")
      options$subfolds <- 3
    }
  }
  
  # fit an appropriate lasso model with cv validation
  mod <- glmnet::cv.glmnet(X, y, type.measure = options$error.measure, alpha = options$alpha,
                           family = fam, type.multinomial = "grouped",
                           nfolds = options$subfolds, ...)
  
  # ID best fit model
  indices <- which(mod$nzero > 0) # lambdas that generated atleast one non-zero coefficients
  lambdas <- mod$lambda[indices] # lambda value 
  scores <- mod$cvm[indices] # mean cv score at that lambda
  best <- which.min(scores) # if a tie, chooses the one with fewer features
  mod_coeffs <- coef(mod, s = lambdas[best]) # ths can maybe be replaced with s = 'lambda.min'
  
  if (fam == "multinomial") {
    # if the data has multiple responses, the coefficients are a matrix
    # that is returned as a list of columns. type.multinomial = "grouped"
    # forced features to be selected for all responses or for none, so we
    # can get the selected features also by only considering the first
    # column. we just replace lasso_coeffs by this column and proceed as usual
    mod_coeffs <- mod_coeffs[[1]]
  }
  
  mod_coeffs <- mod_coeffs[-1,] # remove intercept
  mod_coeffs <- mod_coeffs[which(mod_coeffs != 0)] # remove non-zero coefs 
  
  return (names(mod_coeffs))
}


feature_select_iter <- function(X, y, options = list(n_trials = 1000, threshold = 0.5, force_select = TRUE, 
                                                     alpha = 0.6, subfolds = 5, error.measure = 'mse'),
                                par = T, n.cores = 16, ...) {
  
  # ----------------- BEGIN OPTIONS ----------------- #
  # How often it should be repeated
  if (!("n_trials" %in% names(options))) {
    options$n_trials <- 100
  }
  if (!("threshold" %in% names(options))) {
    options$threshold <- 0.8
  }
  # returns the whole data frame of how often a feature was selected
  if (!("force_select" %in% names(options))) {
    options$force_select <- TRUE
  }
  if (!("alpha" %in% names(options))) {
    options$alpha <- 1
  }
  if (!("subfolds" %in% names(options))) {
    options$subfolds <- 3
  }
  if (!("error.measure" %in% names(options))) {
    options$error.measure <- 'mse'
  }
  # ----------------- END OPTIONS ----------------- #
  
  # vector counting how often each feature is selected
  if (is.null(colnames(X))){colnames(X)<-sapply(1:dim(X)[[2]], function(x) paste0('V', x))}
  feature_count <- rep(0, ncol(X))
  names(feature_count) <- colnames(X)
  
  # run the feature selector trials times and increment the counters
  # for the features that are selected
  if (!par){
    for (trial in 1:options$n_trials) {
      features <- feature_select_reg(X, y, options = list(alpha = options$alpha, subfolds = options$subfolds, 
                                                          error.measure = options$error.measure), ...)
      feature_count[features] <- 1 + feature_count[features]
    }
  }
  else{
    n.cores =availableCores()# number of cores
    RhpcBLASctl::blas_set_num_threads(n.cores) # limit core usage
    
    
    #RhpcBLASctl::blas_set_num_threads(round(n.cores/2)) # limit core usage
    cl <- makeCluster(n.cores)
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = options$n_trials, style = 3)
    progress <- function(n_) setTxtProgressBar(pb, n_)
    opts <- list(progress = progress)
    
    features.all = foreach(trial = 1:options$n_trials, .combine = 'c', .packages = c('glmnet'), .export = c('feature_select_reg'),  
                           .verbose = TRUE, .options.snow = opts) %dopar% {
                             features<-feature_select_reg(X, y, options = list(alpha = options$alpha, subfolds = options$subfolds, 
                                                                               error.measure = options$error.measure), ...)
                           }
    close(pb)
    stopCluster(cl)
    
    for (features in features.all) {
      feature_count[features] <- 1 + feature_count[features]
    }
  }
  
  selected <- feature_count[-which(feature_count <= options$threshold * options$n_trials)] # those above thresh
  if (length(selected) == 0 & options$force_select) {
    selected <- feature_count[which(feature_count == max(feature_count))]
  }
  #return(selected)
  return(names(selected))
}


visualize_validate <- function(vals, options = list()) {
  
  # ----------------- OPTIONS ----------------- #
  if (!("y_label" %in% names(options))) {
    options$y_label <- "score"
  }
  # ----------------- END OPTIONS ----------------- #
  
  
  tmp_vals <-  unlist(vals)
  tmp_vals <- tmp_vals[which(grepl("score", names(tmp_vals)))]
  df_val <- data.frame(score = tmp_vals, model = gsub("_.*", "", names(tmp_vals)))
  df_val$model <- factor(df_val$model, levels = unique(df_val$model))
  
  # assign x-labels
  x_labels <- levels(df_val$model)
  x_labels[which(x_labels == "cv")] <- "model"
  x_labels[which(x_labels == "rf")] <- "random \nfeatures"
  x_labels[which(x_labels == "pt")] <- "permuted \nlabels"
  
  
  # Calculate p-values and generate label
  if (is.null(names(vals))) {
    n_repl <- length(vals)
  } else {
    n_repl <- 1
  }
  
  if (n_repl > 1) {
    # random features
    if ("rf" %in% levels(df_val$model)) {
      pval_rf <- rep(NA, length = length(vals))
      
      for (ind in 1:length(vals)) {
        pval_rf[ind] <- length(which(vals[[ind]]$rf_scores > vals[[ind]]$cv_score))/length(vals[[ind]]$rf_scores)
      }
      
      if (median(pval_rf) < 1/length(vals[[ind]]$rf_scores)) {
        label_rf <- paste0("p<", 1/length(vals[[ind]]$rf_scores))
      } else {
        label_rf <- paste0("p=", median(pval_rf))
      }
    }
    
    # permuted labels
    if ("pt" %in% levels(df_val$model)) {
      pval_pt <- rep(NA, length = length(vals))
      
      for (ind in 1:length(vals)) {
        pval_pt[ind] <- length(which(vals[[ind]]$pt_scores > vals[[ind]]$cv_score))/length(vals[[ind]]$pt_scores)
      }
      
      if (median(pval_pt) < 1/length(vals[[ind]]$pt_scores)) {
        label_pt <- paste0("p<", 1/length(vals[[ind]]$pt_scores))
      } else {
        label_pt <-  paste0("p=", median(pval_pt))
      }
    }
  } else {
    # random features
    if ("rf" %in% levels(df_val$model)) {
      pval_rf <- length(which(vals$rf_scores > vals$cv_score))/length(vals$rf_scores)
      if (pval_rf == 0) {
        label_rf <- paste0("p<", 1/length(vals$rf_scores))
      } else {
        label_rf <- paste0("p=", pval_rf)
      }
    }
    
    # permuted labels
    if ("pt" %in% levels(df_val$model)) {
      pval_pt <- length(which(vals$pt_scores > vals$cv_score))/length(vals$pt_scores)
      if (pval_pt == 0) {
        label_pt <- paste0("p<", 1/length(vals$pt_scores))
      } else {
        label_pt <-  paste0("p=", pval_pt)
      }
    }
  }
  
  
  y_pos <- max(df_val$score) + 0.05
  
  data_summary <- function(x) {
    m <- mean(x)
    ymin <- m - sd(x)
    ymax <- m + sd(x)
    return(c(y = m, ymin = ymin, ymax = ymax))
  }
  
  plt <- ggplot2::ggplot(df_val, ggplot2::aes(x = model, y = score), fill = "gray") +
    ggplot2::geom_violin(fill = "gray", color = "gray") +
    ggplot2::stat_summary(fun.data = data_summary, geom = "pointrange", size = 0.6, fatten = .8, color = "black") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.text = ggplot2::element_text(color = "black", size = 6),
                   axis.title = ggplot2::element_text(color = "black", size = 8),
                   axis.text.x =  ggplot2::element_text(size = 8, angle = 0, hjust = 0.5)) +
    ggplot2::ylab(options$y_label) +
    ggplot2::scale_x_discrete("", labels = x_labels)
  
  if ("rf" %in% levels(df_val$model)) {
    plt <- plt +  ggpubr::geom_bracket(xmin = 1, xmax = which(levels(df_val$model) == "rf"),
                                       inherit.aes = FALSE, label.size = 2.5,
                                       y.position = y_pos, label = label_rf)
  }
  
  if ("pt" %in% levels(df_val$model)) {
    plt <- plt +  ggpubr::geom_bracket(xmin = 1, xmax = which(levels(df_val$model) == "pt"),
                                       inherit.aes = FALSE, label.size = 2.5,
                                       y.position = y_pos + 0.12, label = label_pt)
  }
  
  if (grepl("ccuracy", options$y_label)) {
    plt <- plt + ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1),
                                             labels = c("0", "0.5", "1"),
                                             limits = c(0, max(1, max(df_val$score) + 0.22)))
  }
  return(plt)
  
}