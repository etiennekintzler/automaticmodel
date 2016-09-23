# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'



#rajouter l'option offset, si NULL laisser dans les explicatives mais ne pas la #mettre en offset,
#si pas NULL l'ajouter en offset

#' Finding the best lasso model
#'
#' Finding the best model using 5-fold cross-validation using Lasso-GLM model.
#'
#'#' @param data the dataset
#' @param type.measure the performance criterion use to perform cross-validation. Available values are 'mse', 'mae' (for regression), 'auc' (for classification) and deviance.
#' @param family distribution of the target. Available family are 'gaussian', 'gamma', 'binomial', 'multinomial', 'cox', 'mgaussian'.
#' @param offset name of the offset (if used)
#' @param lambda lambda used to use the best model. 'lambda.min' will choose lambda that optimize cross-validation performances. 'lambda.1se' will choose the lambda that
#' @param train.fraction fraction of data used to find the best model
#' @param seq.lambda user provided sequences of lambda on which to performs the cross-validation
#' @return \code{lassoSelection} returns an object of class "LassoGLM", a list consisting
#' \item{cvfit}{the returns object of the function \code{cv.glmnet} or \code{cv.HDtweedie} }
#' \item{formula}{the formula of the best model found}
#' \item{other}{to fill}
#' \item{prediction}{the predicted value by the best model}
#' \item{y}{true value of the target}
#' \item{selected_features}{vector of character of the selected features}
#' \item{non_selected_features}{vector of character of the non selected features}
#' @importFrom magrittr "%>%"
#' @seealso \code{\link[HDtweedie]{cv.HDtweedie}}, \code{\link[glmnet]{cv.glmnet}}
#' @examples
#' data = cbind(subset(x = mtcars, select = -mpg),
#'              y = mtcars$mpg)
#' lassoSelection(data, type.measure = 'mse', family = 'gaussian', lambda = 'lambda.min', train.fraction = 1)
#' @export
lassoSelection <- function(data = data, type.measure = "mse",
                           family = "gaussian", offset = NULL,
                           lambda = "lambda.min", train.fraction = 1,
                           seq.lambda = NULL){
  data <- na.omit(data)
  Xy <- model.matrix( ~ . - 1 , data = data)
  train <- sample(nrow(data), floor(train.fraction * nrow(data)))

  if(family == 'gamma'){
    cv.lasso   <- HDtweedie::cv.HDtweedie
    cv.control <- list(pred.loss = type.measure, p = 2)
  } else {
    cv.lasso   <- glmnet::cv.glmnet
    cv.control <- list(type.measure = type.measure, family=family)
  }

  if(! is.null(seq.lambda)) {
    cv.control <- c(cv.control, lambda = list(seq.lambda))
  }
  if (is.null(offset)) {

    ### Launching Cross validation LASSO ###
    cvfit <- do.call(cv.lasso, c(list(x = subset(Xy[train, ], select = -c(y)),
                                      y = Xy[train, "y"]), cv.control))
    #plot(cvfit)
    best.model <- coef(cvfit, s = lambda)

    #retirer l'exposure des var expl quand offset == NULL

    ### string work ###
    if (family == 'gamma'){
      i.best <- which(best.model[, 1]!=0)
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
    #cat(paste0("\t", best.features, "\n"))
    concatenate.features <- paste(best.features, collapse = '+')

  } else {

    ### Launching Cross validation LASSO ###
    cvfit <- do.call(what = cv.lasso, args = c(list(x = subset(Xy[train, ],
                                                               select = -c(get(offset), y)),
                                                    y = Xy[train, "y"],
                                                    offset = Xy[train, offset]),
                                               cv.control))

    best.model <- coef(cvfit, s = lambda)

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
    #cat(paste0("\t", best.features, "\n"))
    concatenate.features <- paste(paste(best.features, collapse = '+'),
                                  paste0("offset(", offset, ")"), sep = "+")
  }

  formula.best <- paste('y', concatenate.features, sep = '~')
  pred <- predict(cvfit, newx = subset(Xy, select = -y))
  true <- Xy[, 'y']
  output <- list(formula = formula.best, cvfit = cvfit, prediction = pred,
                 y = true, selected_features = best.features,
                 non_selected_features = non.selected.features)

  class(output) <- 'LassoGLM'
  return(output)
}

#' plot the cross-validation curve of Lasso-GLM procedure
#'
#' Plotting of the performances of the model according to lambda (LASSO regularization term)
#' @param self a LassoGLM object e.g the output of lassoSelection procedure
#'
#' @export
plotLambda <- function(x, ...) UseMethod('plotLambda')
plotLambda.LassoGLM <- function(self) {
  plot(self$cvfit)
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

summary.LassoGLM <- function(self){
  max.length                         <- max(length(self$non_selected_features),
                                            length(self$selected_features))
  length(self$selected_features)     <- max.length
  length(self$non_selected_features) <- max.length
  self$selected_features[is.na(self$selected_features)] <- ''
  self$non_selected_features[is.na(self$non_selected_features)] <- ''
  knitr::kable(data.frame('Selected Features' = self$selected_features,
                          'Not Selected Features' = self$non_selected_features),
               caption = 'Features Selection', format = 'pandoc')
}



perf <- function(x, ...) UseMethod('perf')
perf.LassoGLM <- function(self)
{
  rmse <- rmse(pred = self$prediction, true = self$y)
  gini <- axaml::kpi_gini(predrisk = self$prediction, truerisk = self$y )
  knitr::kable(data.frame('Root Mean Squared error' = rmse,
                          'Gini index' = as.numeric(gini)),
               caption = 'Performances of Lasso-GLM model')
  #invisible(list('Root_mean_squared_error' = rmse, 'Gini_index' = as.numeric(gini)))
}

