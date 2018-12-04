#' Predict the probability of homogolous recombination deficiency using mutational signatures
#'
#' @description A wrapper for predict.randomForest() from the randomForest package
#'
#' @param newdata A data frame or matrix containing new data
#' @param raw.output If TRUE, the probabilities of 'BRCA1', 'BRCA2' and 'none' (no predicted
#' mutation) will also be returned
#'
#' @return A dataframe containing the prediction probabilities
#' @export
#'
chordPredict <- function(newdata, raw.output = F){
   pred <- as.data.frame(predict(newdata=newdata, object=CHORD, type='prob'))
   pred$hrd_prob <- pred$BRCA1 + pred$BRCA2

   if(raw.output){ return(pred) }
   else { return(pred['hrd_prob']) }
}

