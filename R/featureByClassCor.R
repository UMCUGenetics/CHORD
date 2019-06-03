#' T-test between positive vs. negative class
#'
#' @param df A dataframe containing the features and a response (class label) column
#' @param feature Name of the feature column 
#' @param colname.response Name of the response column
#' @param p.value.p Positively correlated features with p-values higher than this are removed
#' @param p.value.n Negatively correlated features with p-values higher than this are removed
#' @param p.responses Name(s) of the positive class
#' @param n.response Name of the negative class
#'
#' @return A list containing a character vector indicating which features are to be kept; and the
#' t-test results
#' @export

featureByClassCor <- function(
  df, feature, colname.response = 'response',
  p.value.p = 0.01, p.value.n = 0.01,
  p.responses = c('BRCA1','BRCA2'), n.response = 'none'
){
  #feature = 'del.mh'
  
  feature_values <- list(
    p = sapply(p.responses, function(i){ df[df[,colname.response] == i, feature] }),
    n = df[df[,colname.response] == n.response, feature]
  )
  
  stats <- do.call(rbind,lapply(feature_values$p, function(i){
    list(
      t_test = t.test(i, feature_values$n)$p.value,
      diff = median(i) - median(feature_values$n)
    )
  }))
  
  if( all(stats[,'diff'] < 0) & all(stats[,'t_test'] < p.value.n) ){ ## remove neg cor features; !redundant
    keep <- 0
  } else if( any(apply(stats, 1, function(i){ i$diff > 0 & i$t_test < p.value.p })) ){
    keep <- 1
  } else {
    keep <- 0
  }
  
  list(keep = structure(keep, names=feature),
       stats = stats)
}