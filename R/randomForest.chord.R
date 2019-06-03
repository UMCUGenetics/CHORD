#' Default randomForest() settings for training CHORD
#'
#' @param df A dataframe containing the features and a response (class label) column
#' @param colname.response Name of the response column
#' @param ... Other arguments that can be passed to randomForest()
#'
#' @return A random forest object
#' @export
#'
randomForest.chord <- function(df=hmf_training_set, colname.response='response', ...){
  randomForest(
    x = df[,colnames(df) != colname.response],
    y = df[,colname.response],
    strata = y, importance = T, proximity = F, ...)
}