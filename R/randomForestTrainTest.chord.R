#' Random forest training and testing procedure 
#'
#' @param train Training set
#' @param test Test set
#' @param colname.response Name of the response column
#' @param do.feature.prefilter Do t-test filtering?
#' @param feature.prefilter.p.value.n Negatively correlated features with p-values higher than this 
#' are removed
#' @param feature.prefilter.p.value.p Positively correlated features with p-values higher than this 
#' are removed
#' @param do.boruta Do boruta feature selection?
#' @param do.balance Do class scaling?
#' @param inner.cv.repeats Number of inner cross-validation repeats
#' @param train.test.seed Set seed
#' @param train.test.verbose Print messages?
#' @param export.dir Dir to output results
#' @param return.output Return output in R?
#' @param ... Other arguments that can be passed to randomForest.chord()
#'
#' @return A random forest and the test predictions
#' @export
#'
randomForestTrainTest.chord <- function(
  train, test = NULL, colname.response = 'response',
  do.feature.prefilter = T, feature.prefilter.p.value.n = 1, feature.prefilter.p.value.p = 1,
  do.boruta = T, do.balance = T, balance.aucpr.tol = NULL,
  inner.cv.repeats = 10,
  train.test.seed = NULL, train.test.verbose = T,
  export.dir = NULL, return.output = F, ...
){
  
  #train=hmf_training_set
  
  if(!is.null(train.test.seed)){ set.seed(train.test.seed) }
  
  if(!is.null(export.dir)){ mkdir(export.dir) }
  
  #--------- Feature pre filter ---------#
  if(do.feature.prefilter){
    if(train.test.verbose){ message('Performing univariate feature filtering...') }
    
    features <- colnames(train)[colnames(train) != colname.response]
    
    feature_by_class_cor <- lapply(features, function(i){
      #i = 'del.mh'
      featureByClassCor(
        df = train, feature = i, 
        p.value.p = feature.prefilter.p.value.p,
        p.value.n = feature.prefilter.p.value.n
      )
    })
    names(feature_by_class_cor) <- features
    
    ## Reorder list
    feature_by_class_cor_2 <- list(
      stats = lapply(feature_by_class_cor, function(i){ i$stats }),
      keep = sort( unlist(lapply(feature_by_class_cor, function(i){ unname(i[['keep']]) })) ,decreasing = T)
    )
    
    if(!is.null(export.dir)){ 
      saveRDS(feature_by_class_cor_2, file.path(export.dir, 'feature_by_class_cor.rds')) 
    }
    
    features_keep <- names(feature_by_class_cor_2$keep[feature_by_class_cor_2$keep == T])
    train <- train[,c(features_keep, colname.response)]
  }
  
  #--------- Boruta feature filter ---------#
  if(do.boruta){
    if(train.test.verbose){ message('Performing boruta feature selection...') }
    boruta_out <- Boruta(
      x = train[,colnames(train) != colname.response],
      y = as.factor(train[,colname.response]), maxRuns=20, doTrace=T ## Low max runs makes feature selection more stringent
    )
    
    if(!is.null(export.dir)){
      saveRDS(boruta_out, file.path(export.dir, '/boruta_out.rds'))
    }
    
    boruta_whitelist <- boruta_out$finalDecision %>% .[. == 'Confirmed'] %>% names()
    train <- train[,c(boruta_whitelist, colname.response)]
  }
  
  #--------- Grid search class scaling ---------#
  if(do.balance){
    if(train.test.verbose){ message('Performing grid search for class balancing...') }
    
    ## Prep grid
    search_grid <- expand.grid(list(none=c(1, 0.5, 0.25), BRCA1 = c(1, 1.5, 2)))
    balance_method <- 'simple'
    
    balance_options <- lapply(1:nrow(search_grid), function(i){
      i = unlist(search_grid[i,])
      list(
        c(names(i[1]), format(round(i[[1]], 2), nsmall = 2), balance_method),
        c(names(i[2]), format(round(i[[2]], 2), nsmall = 2), balance_method)
      )
    })
    
    names(balance_options) <- lapply(balance_options, function(i){
      v <- unlist(i)
      paste(v[v != balance_method], collapse = '_')
    })
    
    ## Do CV
    l_train_balanced <- lapply(balance_options, function(i){
      ## Resampling
      train_new <- balanceClassesMulti(train, colname.response, balance.options = i)
      return(train_new)
    })
    
    # l_rf_cv <- list()
    # for(i in 1:length(l_train_balanced)){
    #   #i = 1
    #   iter_name <- names(balance_options)[[i]]
    #   if(train.test.verbose){ message(' Executing repeated CV for: ', iter_name) }
    #   l_rf_cv[[iter_name]] <- randomForestCv.chord(
    #     train, cv.repeats = inner.cv.repeats,
    #     cv.keep.forests = F, cv.verbose = F
    #   )
    # }
    
    counter <- 0
    l_rf_cv <- lapply(l_train_balanced, function(i){
      counter <<- counter+1
      #i=l_train_balanced[[1]]
      if(train.test.verbose){ message(' Executing repeated CV for: ', names(balance_options)[[counter]]) }
      randomForestCv.chord(
        train, cv.repeats = inner.cv.repeats,
        cv.keep.forests = F, cv.verbose = F
      )
    })
    
    if(!is.null(export.dir)){ saveRDS(l_rf_cv, file.path(export.dir, 'l_rf_cv.rds')) }
    
    ## Get performance
    #l <- l_rf_cv$none_1.00_BRCA1_1.00
    perfCvRep <- function(l){
      perfs <- lapply(l, function(i){
        #i <- l_rf_cv$none_1.00_BRCA1_1.00[[1]]
        confusion <- confusionMatrix(
          i$BRCA1 + i$BRCA2,
          toBinaryResponse(i$response, c('BRCA1','BRCA2'), 1, 'none', 0)
        )
        
        m_pr <- calcPerfCompound(confusion,'pr', metric.names.as.x.y=T)
        m_roc <- calcPerfCompound(confusion,'roc', metric.names.as.x.y=T)
        
        c(
          aucpr = calcAUC(m_pr[,'x'],m_pr[,'y']),
          aucroc = calcAUC(m_roc[,'x'],m_roc[,'y'])
        )
      })
      
      perfs <- as.data.frame(do.call(cbind, perfs))
      perfs$mean <- rowMeans(perfs)
      
      return(perfs)
    }
    
    # l_balance_perfs <- list()
    # for(i in 1:length(l_rf_cv)){
    #   iter_name <- names(l_rf_cv)[[i]]
    #   if(train.test.verbose){ message('  Calculating performance for: ', iter_name) }
    #   l_balance_perfs[[iter_name]] <- perfCvRep(l_rf_cv[[i]])
    # }
    
    counter <- 0
    l_balance_perfs <- lapply(l_rf_cv, function(i){
      counter <<- counter+1
      if(train.test.verbose){ message('  Calculating performance for: ', names(l_rf_cv)[[counter]]) }
      perfCvRep(i)
    })
    
    balance_perfs <- t(as.data.frame( lapply(l_balance_perfs, function(i){ i['mean'] }) ))
    rownames(balance_perfs) <- names(l_balance_perfs)
    balance_perfs <- balance_perfs[order(rownames(balance_perfs)),]
    
    if(!is.null(export.dir)){
      saveRDS(l_balance_perfs, file.path(export.dir, 'l_balance_perfs.rds'))
      saveRDS(balance_perfs, file.path(export.dir, 'balance_perfs.rds'))
    }
    #saveRDS(balance_perfs, paste0(out_dir, 'balance_perfs.rds'))
    
    ## Determine best balance options
    best_balance_option <- names(which.max(balance_perfs[,'aucpr']))
    
    if(!is.null(export.dir)){ saveRDS(best_balance_option, file.path(export.dir, 'best_balance_option.rds')) }
    
    train <- balanceClassesMulti(
      train, colname.response,
      balance.options = balance_options[[best_balance_option]]
    )
  }
  
  #--------- Train 'final' model ---------#
  if(!is.null(train.test.seed)){ set.seed(train.test.seed) }
  
  if(train.test.verbose){ message('Training model with selected features and params...') }
  
  RF <- randomForest.chord(df = train, ...)
  #RF <- randomForest.chord(df = train)
  if(!is.null(test)){
    test[is.na(test)] <- 0
    pred <- as.data.frame(predict(
      object = RF,
      newdata = test[,colnames(test) != colname.response],
      type = "prob"
    ))
    pred$response <- test$response
    
  } else {
    pred <- NA
  }
  
  rf_out <- list(RF = RF, pred = pred)
  
  if(!is.null(export.dir)){ saveRDS(rf_out, file.path(export.dir, 'rf_out.rds')) }
  if(return.output){ return(rf_out) }
}