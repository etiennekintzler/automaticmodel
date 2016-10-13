# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# devtools::missing_s3()



#rajouter l'option offset, si NULL laisser dans les explicatives mais ne pas la #mettre en offset,
#si pas NULL l'ajouter en offset

#' Finding the best lasso model
#'
#' Finding the best model using 5-fold cross-validation using Lasso-GLM model.
#'
#' @param data the dataset used. It must include the target. The offset must be removed if presents.
#' @param type.measure loss to use for the cross-validation. Available values are 'mse', 'mae' (for regression), 'auc' (for classification) and 'deviance'. Default is 'deviance'.
#' @param family distribution of the target. Available family are 'gaussian', 'gamma', 'binomial', 'multinomial', 'cox', 'mgaussian'.
#' @param offset (optionnal) : the vector of the offset values. It should be the same length as the number of rows of the dataset supplied.
#' @param lambda.type lambda used to use the best model. 'lambda.min' will choose lambda that optimize cross-validation performances.
#' 'lambda.1se' gives the most regularized model such that error is within one standard error of the minimum. Default is 'lambda.min'.
#' @param seq.lambda (optionnal) user supplied lambda sequence
#' @param nlambda the number of \code{lambda} values. Default is 50.
#' @param target the name of the target.
#' @param nfolds The number of folds used for k-folds validation. Default is 5.
#' @return \code{lassoSelection} returns an object of class "LassoGLM", a list consisting
#' \item{cvfit}{the returns object of the function \code{cv.glmnet} or \code{cv.HDtweedie} }
#' \item{formula}{the formula of the best model found}
#' \item{other}{to fill}
#' \item{prediction}{the predicted value by the best model}
#' \item{y}{true value of the target}
#' \item{selected_features}{vector of character of the selected features}
#' \item{non_selected_features}{vector of character of the non selected features}
#' @seealso \code{\link[HDtweedie]{cv.HDtweedie}}, \code{\link[glmnet]{cv.glmnet}}
#' @examples
#' data <- cbind(subset(x = mtcars, select = -mpg), y = mtcars$mpg)
#' lassoSelection(data, type.measure = 'mse', family = 'gaussian', lambda = 'lambda.min')
#' @export
lassoSelection <- function(data         = data,
                           target       = NULL,
                           family       = "gaussian",
                           type.measure = 'deviance',
                           offset       = NULL,
                           lambda.type  = "lambda.min",
                           seq.lambda   = NULL,
                           nlambda      = 50,
                           nfolds       = 5){
  data  <- na.omit(data)
  Xy    <- model.matrix( ~ . - 1 , data = data)

  if(family == 'gamma'){
    cv.lasso      <- HDtweedie::cv.HDtweedie
    plot.lasso    <- HDtweedie:::plot.cv.HDtweedie
    predict.lasso <- HDtweedie:::predict.cv.HDtweedie
    cv.control    <- list(pred.loss = type.measure, p = 2, nfolds = nfolds, nlambda = nlambda)
  } else {
    cv.lasso      <- glmnet::cv.glmnet
    plot.lasso    <- glmnet::plot.cv.glmnet
    predict.lasso <- glmnet::predict.cv.glmnet
    cv.control    <- list(type.measure = type.measure, family=family, nfolds = nfolds, nlambda = nlambda)
  }

  if(family == 'gamma' && !is.null(offset) ){
    stop('offset not allowed for gamma family')
  }
  if(!is.null(offset)){
    if( typeof(offset) != 'double' && typeof(offset) != 'integer'){
      stop("offset must be of type 'double' or 'integer'")
    }
  }

  if(! is.null(seq.lambda)) {
    cv.control <- c(cv.control, lambda = list(seq.lambda))
  }
  ifelse(family == 'gamma',
         print('Using function HDtweedie::cv.HDtweedie to perform LASSO-GLM procedure'),
         print('Using function glmnet::cv.glmnet to perform LASSO-GLM procedure'))

  if (is.null(offset)) {

    ### Launching Cross validation LASSO ###
    cvfit <- do.call(cv.lasso, c(list(x = subset(Xy, select = -get(target)),
                                      y = Xy[, target]), cv.control))
    #plot(cvfit)
    best.model <- coef(cvfit, s = lambda.type)

    #retirer l'exposure des var expl quand offset == NULL

    ### string work ###
    if (family == 'gamma'){
      i.best         <- which(best.model[, 1]!=0)
      matching.index <- !is.na(charmatch(colnames(data), names(i.best)))
    } else {
      i.best          <- best.model@i
      matching.index  <- !is.na(charmatch(colnames(data),
                                          best.model@Dimnames[[1]][i.best + 1]))
    }

    best.features         <- colnames(data)[matching.index]
    non.selected.features <- colnames(data)[!matching.index]
    paste0('#features', length(best.features))
    if (length(best.features) < 1) {
      stop('Lasso has not selected any feature')
    }
  } else {
    if(length(offset) != nrow(data)){
      stop("'offset' must have the same length as the number of rows in the dataset")
    }
    ### Launching Cross validation LASSO ###
    cvfit <- do.call(what = cv.lasso, args = c(list(x = subset(Xy, select = -get(target)),
                                                    y = Xy[, target],
                                                    offset = offset),
                                               cv.control))
    best.model <- coef(cvfit, s = lambda.type)

    ### string work ###
    if (family == 'gamma'){
      i.best <- which(best.model[, 1]!=0)
      matching.index <- !is.na(charmatch(colnames(data), names(i.best)))
    } else {
      i.best         <- best.model@i
      matching.index <- !is.na(charmatch(colnames(data),
                                         best.model@Dimnames[[1]][i.best + 1]))
    }
    best.features         <- colnames(data)[matching.index]
    non.selected.features <- colnames(data)[!matching.index]

    paste0('#features', length(best.features))

    if (length(best.features) < 1) {
      stop('Lasso has not selected any feature')
    }
  }

  formula.best <- paste(target, paste(best.features, collapse = '+'), sep = '~')
  if (!is.null(offset)){
    pred <- predict.lasso(cvfit, newx = subset(Xy, select = -get(target)),
                          offset = offset, s = lambda.type)
  } else {
    pred <- predict.lasso(cvfit, newx = subset(Xy, select = -get(target)), s = lambda.type)
  }
  true <- Xy[, target]
  output <- list(formula = formula.best, cvfit = cvfit, prediction = pred,
                 y = true, selected_features = best.features,
                 non_selected_features = non.selected.features)
  class(output) <- 'LassoGLM'
  return(output)
}

#' Plotting of the best Lasso selected model
#'
#' Plots the Gini gain curve the best model selected by Lasso-GLM procedure
#' @param self a LassoGLM object e.g the output of lassoSelection procedure
#'
#' @export
plot.LassoGLM <- function(self) {
  axaml::plotgain(predrisk = self$prediction, truerisk = self$y, return_val = T)$plot
}
#' Selected features of Lasso GLM procedure
#'
#' Give the selected and the non selected features
#' @export
summary.LassoGLM <- function(self){
  max.length                         <- max(length(self$non_selected_features),
                                            length(self$selected_features))
  length(self$selected_features)     <- max.length
  length(self$non_selected_features) <- max.length
  self$selected_features[is.na(self$selected_features)] <- ''
  self$non_selected_features[is.na(self$non_selected_features)] <- ''
  knitr::kable(data.frame('Selected Features'     = self$selected_features,
                          'Not Selected Features' = self$non_selected_features),
               caption = 'Features Selection', format = 'pandoc')
}



#' Give a summary of the performances of the Lasso-GLM procedure.
#'
#' @param format the format of the knitr table returns
#' @return a table with the RMSE and the Gini index of the selected model.
#' @export
perf.LassoGLM <- function(self, format = 'pandoc')
{
  rmse <- rmse(pred = self$prediction, true = self$y)
  gini <- axaml::kpi_gini(predrisk = self$prediction, truerisk = self$y, significance = 4 )
  knitr::kable(data.frame('Root Mean Squared error' = rmse,
                          'Gini index' = as.numeric(gini)),
               caption = 'Performances of Lasso-GLM model',
               format = format)
  #invisible(list('Root_mean_squared_error' = rmse, 'Gini_index' = as.numeric(gini)))
}


#' plot the cross-validation curve of Lasso-GLM procedure
#'
#' Plotting of the performances of the model according to lambda (LASSO regularization term)
#' @param self a LassoGLM object e.g the output of lassoSelection procedure
#' @export
plotLambda <- function(x, ...) UseMethod('plotLambda')

#' Plotting lambda for Lasso
#' @export
plotLambda.LassoGLM <- function(self) {
  plot(self$cvfit)
}

