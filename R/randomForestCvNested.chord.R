#' Outer cross-validation wrapper for  the training and testing procedure 
#'
#' @param df A dataframe containing the features and a response (class label) column
#' @param colname.response Name of the response column
#' @param k.outer Number of outer cross-validation folds
#' @param parent.export.dir Dir to output results
#' @param rfcv.nested.seed Set seed
#' @param rfcv.nested.verbose Print messages?
#' @param ... Other arguments that can be passed to randomForestTrainTest.chord()
#'
#' @return NULL
#' @export
randomForestCvNested.chord <- function(
  df = hmf_training_set, colname.response = 'response', k.outer = 10,
  parent.export.dir = NULL, rfcv.nested.seed = NULL, rfcv.nested.verbose = T, ...
){
  
  if(!is.null(rfcv.nested.seed)){
    ## Make seed vectors for nested CV. Ensures reproducibility with multithreading.
    train_test_seeds <- as.numeric(paste0(rfcv.nested.seed, 1:k.outer))
  }
  
  folds_outer <- createCvTrainTestSets(
    df, k = k.outer, stratify.by.col = colname.response,
    seed = rfcv.nested.seed
  )
  
  n_cores <- detectCores()
  if(.Platform$GUI == 'RStudio'){
    n_cores <- n_cores - 1
  }
  registerDoMC(n_cores)
  
  foreach(i = 1:k.outer) %dopar% {
    fold <- folds_outer[[i]]
    
    if(rfcv.nested.verbose){ message('################  ', 'NESTED CV, OUTER FOLD ', i, '  ################') }
    
    test_outer <- fold$test
    train_outer <- fold$train
    
    export_dir_fold <- paste0(parent.export.dir, '/fold_', str_pad(i, width = 2, pad = '0'), '/')
    mkdir(export_dir_fold)
    
    randomForestTrainTest.chord(
      train = train_outer,
      test = test_outer,
      train.test.seed = train_test_seeds[i], ## Seed
      export.dir = export_dir_fold,
      ...)
    
    return(NULL)
  }
}