# si n.sample est NULL, sinon faire du sampling
#' @export
computeInteractions2 <- function(gbm.model, data, importance.threshold = 0,
                                n.sample = 1e4, max.select = NULL, interact.threshold = 0)
{
  ### Selecting only the most important variables according to a threshold ###
  summary.gbm   <- summary(gbm.model, plotit = F)
  #print(summary.gbm)
  important.var <- rownames(summary.gbm[summary.gbm$rel.inf> importance.threshold,])
  if (!is.null(max.select) & (length(summary.gbm) > max.select) ) {
    important.var <- summary.gbm[1:max.select, "var"]
  }
  #Computing interaction
  product.names   <- t(combn(important.var, 2))
  list.var        <- split(product.names, seq(nrow(product.names)))
  best.iter       <- gbm::gbm.perf(gbm.model, plot.it = F)
  training.sample <- trainingSample(data, train.fraction = (n.sample/nrow(data)))
  H <- OpenRepGrid::sapply_pb(list.var,
                              function(i) gbm::interact.gbm(x       = gbm.model,
                                                            data    = data[training.sample, ],
                                                            i.var   = i,
                                                            n.trees = best.iter))
  df.temp  <- data.frame(product.names, H)
  mat.name <- df.temp[df.temp[, 3] > interact.threshold, ]
  if (nrow(mat.name) == 0) {
    print("There is no interactions at this threshold")
    return(NULL)
  } else {
  mat.name           <- mat.name[order(mat.name[, 3], decreasing = T), ]
  colnames(mat.name) <- c('Variable 1', 'Variable 2', 'H-statistic')
  row.names(mat.name) <- NULL
  return(mat.name)
  }
}


getInteractions <- function(mat.name, interact.threshold = 0, gbm.model)
{
  ### Selection of the highest interactions according to a threshold ###
  mat.name      <- mat.name[mat.name[, 3] > interact.threshold, ]
  if (nrow(mat.high) == 0) {
    print("There is no interactions at this threshold")
  } else {
    mat.name
  }
  return(mat.name)
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
      if (include) {
      interact.name <- paste(paste(mat.name[, 1], mat.name[, 2], sep = '*'),
                             collapse = "+")
    } else {
      interact.name <- paste(paste(mat.name[, 1], mat.name[, 2], sep = ':'),
                             collapse = "+")
    }

    print(paste0('Interaction(s): ', interact.name, ', number of interactions: ', nrow(mat.name)))

    ### Training the GLMs ###

    #TRAIN index with all factros
    data  <- na.omit(data)
    #train <- trainingSample(data, train.fraction = train.fraction)
    train <- sample(nrow(data), floor(train.fraction*nrow(data)))

    cat("\t Computing Basemodel GLM \n")
    #print(formula.best)
    glm.model <- glm(formula = as.formula(formula.best), data = data[train, ],
                     family = family)

    cat("\t Computing Interaction Model GLM \n")
    formula.plus.interact <- paste(formula.best, interact.name, sep = "+")
    glm.model.plus        <- glm(formula = as.formula(formula.plus.interact),
                                 data = data[train, ], family = family)

    cat("\t Computing the performance measure (RMSE and Gini) \n")

    ### Computing the performance measure ###
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
                            truerisk = true)
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
# mettre SPEED GLM
# Mettre indicateur GINI dans les if else. si package pour plotter/calculer gini facilement, le prendre. Difference simple en terme de GINI
#faire juse GINI et RMSE à chaque fois. SI count.
# Voir pour le calcul de la deviance hors-ech.








getGlmPerformance2 <- function(threshold.values = c(0.01, seq(0.05, 1, by = 0.05)),
                              data, mat.name, gbm.model, formula.best, family,
                              train.fraction = .8, max.interactions = 50, include = F, speed = F)
{
  l  <- sapply(1:length(threshold.values), function(i, x = threshold.values) {tmp <-sum(mat.name[, 3] > x[i]); (tmp > 0) & (tmp < max.interactions)})
  ll <- c(T, sapply(2:length(threshold.values), function(i, x = threshold.values) sum(mat.name[, 3] > x[i]) < sum(mat.name[, 3] > x[i-1])))
  stored.results <- lapply(which(l * ll > 0), function(i) {set.seed(1991); glmPerformance(interact.threshold = threshold.values[i],
                                                                                          train.fraction = train.fraction,
                                                                                          gbm.model = gbm.model, formula.best = formula.best,
                                                                                          mat.name = mat.name, data = data,
                                                                                          family = family, include = include, speed = speed)})
  return(stored.results)
}



#' Detect interactions
#'
#' Using Gradient Boosting Model and Friedman H-statistic to find interactions in the model.
#' Comparison of the performances between the base model and the interactions-augmented model using root mean squared error and gini index.
#' @param data the dataset
#' @param formula.best the formula of the base model provided by the user or found in \code{\link{lassoSelection}}
#' @param gbm.control a list containing the parameters passed to the Gradient Boosting model \itemize{
#' \item{n.sample}{Number of points used to compute H-statistic}
#' }
#'
#' @seealso \code{\link[gbm]{gbm}} for GBM modeling
#' @seealso \code{\link[gbm]{interact.gbm}} to compute H-statitic,
#'
#' @references Friedman, J. and Popescu B. (2005) \emph{Predictive Learning via Rule Ensembles} \url{http://statweb.stanford.edu/~jhf/ftp/RuleFit.pdf}
#' @references Firedman, J. (2009) \emph{Greedy Function Approximation : A Gradient Boosting Machine} \url{https://statweb.stanford.edu/~jhf/ftp/trebst.pdf}
#'
#' @export
detectInteractions <- function(data, max.interactions, formula.best, shuffle = F,
                               gbm.control = list(train.fraction = .8, family, shrinkage = .01, bag.fraction = .5,
                                                  n.trees = 100, n.sample = 5e4, depth = 5, importance.threshold = 0, max.select = 30),
                               glm.control = list(train.fraction = .8, family, include = F, speed = F,
                                                  threshold.values = c(0.01, seq(0.05, 1, by = 0.05))))
{
  if (sum(sapply(data, function(x) length(unique(x)) <= 1)) > 0 ) {
    stop('No variation in at least one of the explanatory variable')
  }
  print("Fitting GBM Model")
  #Fitting GBM model
  gbm.model <- gbm::gbm(y ~ . , distribution = gbm.control$family, data = data, n.trees = gbm.control$n.trees, interaction.depth = gbm.control$depth,
                        shrinkage = gbm.control$shrinkage, bag.fraction = gbm.control$bag.fraction, verbose = T, train.fraction = gbm.control$train.fraction)

  print("Computing interaction matrix")
  #computing friedman H matrix
  mat.name <- computeInteractions2(gbm.model, data = data, importance.threshold = gbm.control$importance.threshold,
                                  n.sample = gbm.control$n.sample, max.select = gbm.control$max.select)
  cat('\n')

  #print(getInteractions(mat.name[, 3] = mat.inter, gbm.model = gbm.model, interact.threshold = 0))
  if (shuffle == T) {
    mat.name[, 3] <- mat.name[sample(nrow(mat.name)), 3]
  }
  print('Fitting GLMs model')
  GLMs.perf <- getGlmPerformance2(threshold.values = glm.control$threshold.values , data = data, mat.name = mat.name,
                               gbm.model = gbm.model, formula.best = formula.best, family = glm.control$family,
                               train.fraction = glm.control$train.fraction, max.interactions = max.interactions,
                               include = glm.control$include, speed = glm.control$speed) #0-1 la famille doit être binomial

  function.call <- as.list(sys.call())
  output        <- list(GLMs.perf = GLMs.perf, best.interactions = mat.name,
                        function.call = function.call)
  class(output) <- 'InteractionsDetection'
  invisible(output)
}

#modifier la fonction pour faire apparaître bernoulli ou binomial


#definition of methods for InteractionsDetection object

#' LR
#'
#' @export
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


#' Return interactions of the best model
#'
#' Return all the interactions of the best model. The best model is selected according to a criterion ('Gini' or 'RMSE')
#' @param criterion Criterion used to select the best model. Values can be 'Gini' for gini index or 'RMSE' for Root Mean Squared Error
#' @return Printing a table with the best interactions and their correspond H-statistic
#' @export
summary.InteractionsDetection <- function(self, criterion = 'Gini', format = 'pandoc')
{
  #This method return the best model according to GINI
  index <- rep(NA, length(self$GLMs.perf))
  for (i in 1:length(index)) {
    if (!is.null(self$GLMs.perf[[i]][[criterion]][3])) {
      index[i] <- self$GLMs.perf[[i]][[criterion]][3]
    }
  }
  index.best <- which.max(index)
  best.model <- self$GLMs.perf[[index.best]]
  interaction.matrix <- best.model$mat
  colnames(interaction.matrix) <- c('Variable 1', 'Variable 2', 'H-statistic')
  knitr::kable(interaction.matrix, caption = paste0('Best model includes ', best.model$options[1],
                                                    ' interactions, H-statistic treshold : ', best.model$options[2],
                                                    sep = ''), format = format)
}

#' Performances of the best selected model
#'
#' Returning the performances (RMSE and Gini) of the best model selected according to criterion
#' @param criterion Criterion used to select the best model. Values can be 'Gini' or 'RMSE'
#' @return Table with the performances of the best model
#' @export
perf <- function(x, ...) UseMethod('perf')
#' @export
perf.InteractionsDetection <- function(self, criterion = 'Gini', format = 'pandoc')
{
  index <- sapply(seq_along(self$GLMs.perf), function(i, x = self$GLMs.perf) x[[i]][[criterion]][3])
  index.best <- which.max(index)
  best.model <- self$GLMs.perf[[index.best]]
  dataframe  <- data.frame('Root\ Mean Squared error' = best.model$RMSE, 'Gini index' = best.model$Gini)
  rownames(dataframe) <- c('Base Model', 'Augmented Model', 'Difference augmented - base ')
  cat1 <- paste0("Performances of the best augmented model vs base model according to '", criterion, "' criterion,", sep = '')
  knitr::kable(dataframe, caption = paste(cat1,paste0('Train fraction:', best.model$options[3]), sep = '\n'),
               format = format)
  #invisible(list('Root_mean_squared_error' = best.model$RMSE, 'Gini_index' = as.numeric(best.model$Gini)))
}


#    return(list(RMSE = res.rmse, Gini = res.gini, likelihood.ratio.test = lr.results, anova.results = anova.results, formula.plus = formula.plus.interact, predictions = pred, plot_gain.control = plot_gain.control))


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
plot.InteractionsDetection <- function(self, i = 1, rev = NULL,data = eval(self$function.call$data),
                                       cut.x = NULL, cut.trace = NULL,
                                       simplify = NULL, criterion = 'Gini')
{
  # This method plot the i-th interactions, by default it plot the best interaction
  index <- sapply(seq_along(self$GLMs.perf), function(i, x = self$GLMs.perf) x[[i]][[criterion]][3])
  index.best <- which.max(index)
  best.model <- self$GLMs.perf[[index.best]]
  plotInteract(x = self$best.interactions[i, 1:3],
               y = predict(best.model$glm.model.plus, eval(self$function.call$data), type = 'response'),
               data = data, rev=rev, cut.x=cut.x, cut.trace=cut.trace, simplify = simplify)
}

# Ne pas oublier les TROIS PETIT POINTS si on utilise plus d'un paramètre
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
  index <- sapply(seq_along(self$GLMs.perf), function(i, x = self$GLMs.perf) x[[i]][[criterion]][3])
  index.best <- which.max(index)
  do.call(what = axaml::plotgain, args = c(self$GLMs.perf[[index.best]]$plot_gain.control, return_val = T))$plot
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
  index <- sapply(seq_along(self$GLMs.perf), function(i, x = self$GLMs.perf) x[[i]][[criterion]][3])
  index.best <- which.max(index)
  best.model <- self$GLMs.perf[[index.best]]
  for(i in 1:best.model$options[1]){
    plot(self, i = i, simplify = simplify)
  }
}
