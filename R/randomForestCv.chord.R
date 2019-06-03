#' Cross-validation wrapper for randomForest.chord()
#'
#' @param df A dataframe containing the features and a response (class label) column
#' @param colname.response Name of the response column
#' @param cv.folds Number of cross-validation folds
#' @param cv.repeats Number of times to repeat the cross-validation
#' @param cv.verbose Print messages?
#' @param cv.keep.forests Keep the random forests or only keep the predictions?
#' @param ... Other arguments that can be passed to randomForest.chord()
#'
#' @return A list of random forests and/or the predictions
#' @export
#'
randomForestCv.chord <- function(
  df, colname.response = 'response', 
  cv.folds = 10, cv.repeats = 1, cv.verbose = T, cv.keep.forests = T, ...
){
  
  folds <- createCvTrainTestSets(df, k = cv.folds, stratify.by.col = colname.response)
  
  rf_cv_rep <- lapply(1:cv.repeats, function(j){
    
    if(cv.verbose){ message('   CV repeat: ', j) }
    
    rf_cv <- lapply(folds, function(i){
      RF <- randomForest.chord(df = i$train, ...)
      
      i$test[is.na(i$test)] <- 0
      pred <- as.data.frame(predict(
        object = RF,
        newdata = i$test[,colnames(i$test) != colname.response],
        type = "prob"
      ))
      pred$response <- i$test$response
      
      return(list(RF=RF, pred=pred))
    })
    
    ## Reorder list
    list(
      RF = lapply(rf_cv, function(i){ i$RF }),
      pred = do.call(rbind, lapply(rf_cv, function(i){ i$pred }))
    )
  })
  
  if(cv.keep.forests){
    rf_cv_rep
  } else {
    lapply(rf_cv_rep, function(i){ i$pred })
  }
}