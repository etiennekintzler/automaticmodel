# si n.sample est NULL, sinon faire du sampling
#' Compute H-statistic for given pairs of variables
#'
#' Compute Friedman H-statistic for every pair of variable
#' @param gbm.model The \code{gbm.model} object used to compute the prediction
#' @param data the data used
#' @param importance.threshold the threshold from which to select the variables which interaction value will be computed, default is 0
#' @param max.select maximum number of variables selected for the computing of interactions
#' @param interact.threshold the threshold from which to keep/retain the interactions, default is 0
#' @param n.sample the number of observations used to compute H-statistic
#'
#' @seealso \code{\link[gbm]{interact.gbm}} for the computation of the H-statistic
#' @references Friedman, J. and Popescu B. (2005) \emph{Predictive Learning via Rule Ensembles} \url{http://statweb.stanford.edu/~jhf/ftp/RuleFit.pdf}
#' @return A dataframe which first two columns corresponds to the interacted variables and the third one corresponds to the H-statistic values for these interactions

#' @export
computeInteractions <- function(gbm.model, data, importance.threshold = 0, max.select = NULL,
                                n.sample = 1e4, interact.threshold = 0)
{
  # Selecting only the most important variables according to a threshold
  summary.gbm   <- summary(gbm.model, plot_it = F)
  important.var <- rownames(summary.gbm[summary.gbm$rel_inf> importance.threshold,])
  if (!is.null(max.select) & (length(summary.gbm) > max.select) ) {
    important.var <- summary.gbm[1:max.select, "var"]
  }
  # Computing all the pairs of interactions
  product.names   <- t(combn(important.var, 2))
  list.var        <- split(product.names, seq(nrow(product.names)))
  best.iter       <- gbm::gbm.perf(gbm.model, plot.it = F)
  training.sample <- trainingSample(data, train.fraction = (n.sample/nrow(data)))
  H               <- OpenRepGrid::sapply_pb(list.var, function(i){
    gbm::interact(gbm_fit_obj = gbm.model,
                  data        = data[training.sample, ],
                  var         = charmatch(i, gbm.model$variables$var_names),
                  num_tree    = best.iter)
  })
  # Formating the output in a dataframe and selected only the pairs which value are superior to interact.threshold
  df.temp  <- data.frame(product.names, H)
  mat.name <- subset(x = df.temp, subset = (H > interact.threshold) & !is.nan(H))
  if (nrow(mat.name) == 0) {
    print("There is no interactions at this threshold")
    return(NULL)
  } else {
  mat.name            <- mat.name[order(mat.name[, 3], decreasing = T), ]
  colnames(mat.name)  <- c('Variable 1', 'Variable 2', 'H-statistic')
  row.names(mat.name) <- NULL
  return(mat.name)
  }
}

glmPerformance <- function(data, gbm.model = gbm.model, formula.best, family,
                           interact.threshold = 0, train.fraction = .8, mat.name,
                           speed = F, include = F)
{
  if (!is.null(speed)) {
    if(speed == T)  glm <- speedglm::speedglm
  }
    mat.name      <- mat.name[mat.name[, 3] > interact.threshold, ]
    if (nrow(mat.name) == 0) {
      print("There is no interactions at this threshold")
      return(NULL)
    } else {
      # Computing the formula with interactions
      if (include) {
      interact.name <- paste(paste(mat.name[, 1], mat.name[, 2], sep = '*'),
                             collapse = "+")
    } else {
      interact.name <- paste(paste(mat.name[, 1], mat.name[, 2], sep = ':'),
                             collapse = "+")
    }

    print(paste0('Interaction(s): ', interact.name, ', number of interactions: ', nrow(mat.name)))

    # Training the GLMs
    data  <- na.omit(data)
    #train <- trainingSample(data, train.fraction = train.fraction)
    train <- sample(nrow(data), floor(train.fraction*nrow(data)))

    cat("\t Computing Basemodel GLM \n")
    glm.model <- glm(formula = as.formula(formula.best),
                     data    = data[train, ],
                     family  = family)

    cat("\t Computing Interaction Model GLM \n")
    formula.plus.interact <- paste(formula.best, interact.name, sep = "+")
    glm.model.plus        <- glm(formula = as.formula(formula.plus.interact),
                                 data    = data[train, ],
                                 family  = family)

    cat("\t Computing the performance measure (RMSE and Gini) \n")

    # Computing the performance criteria
    pred.base <- predict(glm.model, newdata = data[ - train, ], type = "response")
    pred.plus <- predict(glm.model.plus, newdata = data[ - train, ], type = "response")

    true <- data[ - train, "y"]
    rmse.base <- rmse(pred.base, true)
    rmse.plus <- rmse(pred.plus, true)
    res.rmse  <- c("RMSE base" = round(rmse.base, 4), "RMSE plus" = round(rmse.plus, 4),
                   "diff plus-base (%)" = paste0(round((rmse.plus / rmse.base - 1) * 100, 4), " %"))
    options <- c("Number of interactions" = dim(mat.name)[1], "H-statistic threshold" = interact.threshold,
                 "Train fraction" = train.fraction)
    gini <- axaml::kpi_gini(predrisk = list(gini.base = pred.base, gini.plus = pred.plus),
                            truerisk = true,
                            significance = 4)
    res.gini <- c('Gini Base' = gini[1], 'Gini Plus' = gini[2],
                  'Diff plus - base' = as.numeric(gini[2] - gini[1]) )

    lr.results        <- lmtest::lrtest(glm.model, glm.model.plus)
    anova.results     <- NULL #anova(glm.model.plus, test = 'Chisq')
    print(options);print(res.rmse);print(res.gini)
    pred              <- list('base' = pred.base, 'plus' = pred.plus)
    plot_gain.control <- list(predrisk = pred, truerisk = true)

    return(list(RMSE = res.rmse, Gini = res.gini, options = options, likelihood.ratio.test = lr.results,
                anova.results = anova.results, formula.plus = formula.plus.interact, glm.model.plus = glm.model.plus,
                predictions = pred, plot_gain.control = plot_gain.control, mat = mat.name))
    }
}

# Voir pour le calcul de la deviance hors-ech.

allGlmPerformances <- function(threshold.values = c(0.01, seq(0.05, 1, by = 0.05)),
                              data, mat.name, gbm.model, formula.best, family,
                              train.fraction = .8, max.interactions = 50, include = F, speed = F, seed = 2016)
{
  l  <- sapply(1:length(threshold.values),
               function(i, x = threshold.values) {
                 tmp <-sum(mat.name[, 3] > x[i])
                 return((tmp > 0) & (tmp < max.interactions))
                 })
  ll <- c(T, sapply(2:length(threshold.values),
                    function(i, x = threshold.values) sum(mat.name[, 3] > x[i]) < sum(mat.name[, 3] > x[i-1])))
  stored.results <- lapply(which(l * ll > 0),
                           function(i){
                             set.seed(seed)
                             glmPerformance(interact.threshold = threshold.values[i], train.fraction = train.fraction,
                                            gbm.model = gbm.model, formula.best = formula.best,
                                            mat.name = mat.name, data = data, family = family,
                                            include = include, speed = speed)
                             })
  return(stored.results)
}


#' Detect interactions using GBM modeling
#'
#' Using Gradient Boosting Model and Friedman H-statistic to find interactions in the model.
#' Comparison of the performances between the base model and the interactions-augmented model using root mean squared error and gini index.
#' @param data the dataset
#' @param formula.best the formula of the base model provided by the user or found in \code{\link{lassoSelection}}
#' @param gbm.control a list containing the parameters passed to the Gradient Boosting model \itemize{
#' \item{n.sample}{Number of points used to compute H-statistic}
#' \item{train.fraction}{Fraction of the data used to compute OOB error}
#' \item{shrinkage}{Shrinkage or learning rate of \code{gbm}}
#' \item{bag.fraction}{Subsampling rate}
#' \item{family}{distribution of the model, to see supported distribution please see \code{gbm} function}
#' \item{n.trees}{Number of trees of \code{gbm}}
#' \item{n.sample}{Number of observations drawn to compute the H-statistic, higher value slows down the computations}
#' \item{n.sample}{Correponds to \code{interaction.depth} parameter in \code{gbm}}
#' \item{importance.threshold}{The importance threshold at which the variables are selected to be tested as interactions}
#' \item{max.select}{Maximum number of variables selected to be tested as interactions -related to previous item}
#' }
#' @param glm.control a list containing the parameters for the computation of the base model and augmented model using \code{glm}\itemize{
#' \item train.fraction The fraction of the sample on which to train the GLM models
#' \item family The family of the data passed as \code{glm} function
#' \item include whether the user wants to include the univariate component of the interaction, default is False
#' \item speed whether the user wants to include \code{\link[speedglm]{speedglm}}
#' \item threshold.values The threshold values for the H-statistics
#' }
#'
#' @seealso \code{\link[gbm]{glm}} for Generalized Linear Models modeling
#' @seealso \code{\link[gbm]{gbm}} for Gradient Boosting Machine modeling
#' @seealso \code{\link[gbm]{interact.gbm}} to compute H-statitic,
#' @return A list of three object. The first object is a list corresponds to the performances of the regression objects of all
#'  the GLMs. The second is the dataframe with the interacted variables names and the corresponding H-statistic. The last object
#'  corresponds to the \code{function.call}.
#' @references Friedman, J. and Popescu B. (2005) \emph{Predictive Learning via Rule Ensembles} \url{http://statweb.stanford.edu/~jhf/ftp/RuleFit.pdf}
#' @references Firedman, J. (2009) \emph{Greedy Function Approximation : A Gradient Boosting Machine} \url{https://statweb.stanford.edu/~jhf/ftp/trebst.pdf}
#'
#' @export
detectInteractions <- function(data, max.interactions, formula.best, shuffle = F,
                               gbm.control = list(train.fraction = .75, family, shrinkage = .01,
                                                  bag.fraction = .5, n.trees = 100, n.sample = 5e4, depth = 5,
                                                  importance.threshold = 0, max.select = 30),
                               glm.control = list(train.fraction = .75, family, include = F, speed = F,
                                                  threshold.values = c(0.01, seq(0.05, 1, by = 0.05)), seed = 2016))
{
  if (sum(sapply(data, function(x) length(unique(x)) <= 1)) > 0 ) {
    stop('No variation in at least one of the explanatory variable')
  }
  print("Fitting GBM Model")
  # Fitting GBM model
  gbm.model <- gbm::gbm(y ~ . , distribution = gbm.control$family, data = data, n.trees = gbm.control$n.trees,
                        interaction.depth = gbm.control$depth, shrinkage = gbm.control$shrinkage,
                        bag.fraction = gbm.control$bag.fraction, verbose = T,
                        train.fraction = gbm.control$train.fraction)

  print("Computing interaction matrix")

  # Computing H-statistic matrix
  mat.name <- computeInteractions(gbm.model, data = data, importance.threshold = gbm.control$importance.threshold,
                                  n.sample = gbm.control$n.sample, max.select = gbm.control$max.select)
  cat('\n')

  if (shuffle == T) {
    mat.name[, 3] <- mat.name[sample(nrow(mat.name)), 3]
  }
  print('Fitting GLMs model')
  GLMs.perf <- allGlmPerformances(threshold.values = glm.control$threshold.values , data = data,
                                  mat.name = mat.name, gbm.model = gbm.model, formula.best = formula.best,
                                  family = glm.control$family, train.fraction = glm.control$train.fraction,
                                  max.interactions = max.interactions, include = glm.control$include,
                                  speed = glm.control$speed, seed = glm.control$seed)
  function.call <- as.list(sys.call())

  output        <- list(GLMs.perf         = GLMs.perf,
                        best.interactions = mat.name,
                        function.call     = function.call)
  class(output) <- 'InteractionsDetection'
  invisible(output)
}


#############    Definition of methods for InteractionsDetection class   #########################################

#' LR
#'
#'
LR <- function(x, ...) UseMethod('LR')
#' @export
LR.InteractionsDetection <- function(self)
{
  p.val <- rep(NA, length(self$GLMs.perf))
  for (i in 1:length(p.val)) {
    if (!is.null(self$GLMs.perf[[i]]$likelihood.ratio.test$`Pr(>Chisq)`[2])) {
      p.val[i] <- self$GLMs.perf[[i]]$likelihood.ratio.test$`Pr(>Chisq)`[2]
    }
  }
  names(p.val) <- 'p-values'
  knitr::kable(na.omit(p.val))
}


#' @export
returnBest <- function(x, ...) UseMethod('returnBest')
#' @export
returnBest.InteractionsDetection <- function(self, criterion = 'Gini') {
  index <- sapply(seq_along(self$GLMs.perf), function(i, x = self$GLMs.perf) x[[i]][[criterion]][3])
  index.best <- switch(criterion,
                       'Gini' = which.max(index),
                       'RMSE' = which.min(index),
                       stop("Criterion must be 'Gini' or 'RMSE'"))
  best.model <- self$GLMs.perf[[index.best]]
  return(best.model)
}


#' Return interactions of the best model
#'
#' Return all the interactions of the best model. The best model is selected according to a criterion ('Gini' or 'RMSE')
#' @param criterion Criterion used to select the best model. Values can be 'Gini' for gini index or 'RMSE' for Root Mean Squared Error
#' @return Printing a table with the best interactions and their correspond H-statistic
#' @export
summary.InteractionsDetection <- function(self, criterion = 'Gini', format = 'pandoc')
{
  #This method return the best model according to criterion
  best.model <- returnBest(self, criterion)
  mat.name   <- best.model$mat
  knitr::kable(mat.name, caption = paste0('Best model includes ', best.model$options[1],
                                          ' interactions, H-statistic treshold : ',
                                          best.model$options[2], sep = ''),
               format = format)
}


#' Performances of the best selected model
#'
#' Returning the performances (RMSE and Gini) of the best model selected according to criterion
#' @seealso \code{\link[axamodel]{perf.InteractionsDetection}} for the comparison of the performances of the GLM base model and GLM
#' augmented model.
#' @seealso \code{\link[axamodel]{perf.LassoGLM}} for the performances of the GLM-LASSO selected model
#' @export
perf <- function(x, ...) UseMethod('perf')
#' Performances of the best selected model
#'
#' Returning the performances (RMSE and Gini) of the best model selected according to criterion
#' @param criterion Criterion used to select the best model. Values can be 'Gini' or 'RMSE'
#' @param format the format of the knitr table returns
#' @return Table with the performances of the best model
#' @export
perf.InteractionsDetection <- function(self, criterion = 'Gini', format = 'pandoc')
{
  best.model <- returnBest(self, criterion)
  dataframe  <- data.frame('Root\ Mean Squared error' = best.model$RMSE, 'Gini index' = best.model$Gini)
  rownames(dataframe) <- c('Base Model', 'Augmented Model', 'Difference augmented - base ')
  paste1 <- paste0("Performances of the best augmented model vs base model according to '",
                 criterion, "' criterion,", sep = '')
  knitr::kable(dataframe, caption = paste(paste1, paste0('Train fraction:', best.model$options[3]), sep = '\n'),
               format = format)
}


#plot.myclass <- function(self, ...) UseMethod('plot') -- not needed for already defined methods like plot, summary
#' Interaction plot of interactions found
#'
#' Plot the interactions found in \code{\link{detectInteractions}} procedure
#' @param self an \code{InteractionsDetection}.
#' @param i the i-th best interaction according to its \emph{H-statistic} value. Default value is 1.
#' @param data The dataset used to perform the realize interaction plot. Default value is the dataset used in \code{detectInteractions} call.
#' @return An \code{interaction.plot} of the i-th best interaction.
#' @seealso \code{\link[stats]{interaction.plot}}
#' @export
plot.InteractionsDetection <- function(self, i = 1, rev = NULL, data = eval(self$function.call$data),
                                       cut.x = NULL, cut.trace = NULL,
                                       simplify = NULL, criterion = 'Gini')
{
  # This method plot the i-th interactions, by default it plot the best interaction
  best.model <- returnBest(self, criterion)
  plotInteract(x = self$best.interactions[i, 1:3],
               y = predict(best.model$glm.model.plus, data, type = 'response'),
               data = data, rev=rev, cut.x=cut.x, cut.trace=cut.trace, simplify = simplify)
}

#' @export
bestInteractions <- function(x, ...) UseMethod('bestInteractions')
#' @export
bestInteractions.InteractionsDetection <- function(self, i = 5, format = 'pandoc') {
  return(knitr::kable(self$best.interactions[1:i, ], format = format, caption = 'Best interactions'))
}


#' Plot Gini index of the augmented model vs base model
#'
#' This method performs a plot gain (gini index) of the best augmented-model (according to \code{criterion}) vs base model
#'
#' The function plot  all the interactions included in the best augmented-model
#' @return a \code{plot_gain} plot of the gini index of both the best augmented-model and the base model
#' @seealso \code{\link[axaml]{plot_gain}} in \pkg{axaml} package
#' @export
plotBest <- function(x, ...) UseMethod('plotBest')
#' @export
plotBest.InteractionsDetection <- function(self, criterion = 'Gini') {
  do.call(what = axaml::plotgain,
          args = c(returnBest(self, criterion)$plot_gain.control, return_val = T))$plot
}


#' Plot all the interactions included in the best augmented-model
#'
#' The function plot  all the interactions included in the best augmented-model
#'
#' @seealso \code{\link{plot.InteractionsDetection}}
#' @export
plotAll <- function(x, ...) UseMethod('plotAll')
#' @export
plotAll.InteractionsDetection <- function(self, simplify = 8, criterion = 'Gini') {
  for(i in 1:returnBest(self, criterion)$options[1]){
    plot(self, i = i, simplify = simplify)
  }
}
