cleanData <- function(data, n.levels = 20)
{
  novar.col      <- sapply(data, function(x) length(unique(x)) <= 1)
  high.levels    <- sapply(data, function(x) length(levels(x)) > n.levels) #justifier dans rapport
  na.columns     <- sapply(data, function(x) sum(is.na(x)) > 0.2 * nrow(data))
  date.col       <- grepl(pattern = 'date', x = colnames(data))
  keep.columns   <- !(high.levels | na.columns | novar.col | date.col)
  return(data[, keep.columns])
}

#rajouter l'option offset, si NULL laisser dans les explicatives mais ne pas la #mettre en offset,
#si pas NULL l'ajouter en offset

#' Finding the best lasso model
#'
#' this is the description
#'
#' @param data the dataset
#'
#' @return explain the return
#' @importFrom magrittr "%>%"
#' @export
lassoSelection <- function(data, type.measure = "mse", family = "gaussian", offset, lambda = "lambda.min",
                    train.fraction = .8, seq.lambda = NULL){
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
    cvfit      <- do.call(cv.lasso, c(list(x = subset(Xy[train, ], select = -c(y)),
                                           y = Xy[train, "y"]), cv.control))
    #plot(cvfit)
    best.model <- coef(cvfit, s = lambda)

    #retirer l'exposure des var expl quand offset == NULL

    ### string work ###
    if (family == 'gamma'){
      i.best <- which(best.model[, 1]!=0)
      matching.index <- !is.na(charmatch(colnames(data), names(i.best)))
    } else {
      i.best               <- best.model@i
      matching.index       <- !is.na(charmatch(colnames(data), best.model@Dimnames[[1]][i.best + 1]))
    }

    best.features        <- colnames(data)[matching.index]
    non.selected.features <- colnames(data)[!matching.index]

    paste0('#features', length(best.features))
    if (length(best.features) < 1) {
      stop('Lasso has not selected any feature')
    }
    #cat(paste0("\t", best.features, "\n"))
    concatenate.features <- paste(best.features, collapse = '+')

  } else {

    ### Launching Cross validation LASSO ###
    cvfit     <- do.call(cv.lasso, c(list(x = subset(Xy[train, ], select = -c(get(offset), y)),
                                          y = Xy[train, "y"], offset = Xy[train, offset]), cv.control))

    plot(cvfit)
    best.model <- coef(cvfit, s = lambda)

    ### string work ###
    if (family == 'gamma'){
      i.best <- which(best.model[, 1]!=0)
      matching.index <- !is.na(charmatch(colnames(data), names(i.best)))
    } else {
      i.best               <- best.model@i
      matching.index       <- !is.na(charmatch(colnames(data), best.model@Dimnames[[1]][i.best + 1]))
    }
    best.features        <- colnames(data)[matching.index]
    non.selected.features <- colnames(data)[!matching.index]

    paste0('#features', length(best.features))

    if (length(best.features) < 1) {
      stop('Lasso has not selected any feature')
    }
    #cat(paste0("\t", best.features, "\n"))
    concatenate.features <- paste(paste(best.features, collapse = '+'), paste0("offset(", offset, ")"), sep = "+")
  }

  formula.best <- paste('y', concatenate.features, sep = '~')
  pred <- predict(cvfit, newx = subset(Xy, select = -y))
  true <- Xy[, 'y']
  output <- list(formula = formula.best, cvfit = cvfit, prediction = pred, y = true, selected_features = best.features, non_selected_features = non.selected.features)

  class(output) <- 'LassoGLM'
  return(output)
}


plotLambda <- function(x, ...) UseMethod('plotLambda')
plotLambda.LassoGLM <- function(self) {
  plot(self$cvfit)
}

plot.LassoGLM <- function(self) {
  axaml::plotgain(predrisk = self$prediction, truerisk = self$y, return_val = T)$plot
}

summary.LassoGLM <- function(self){
  max.length                         <- max(length(self$non_selected_features), length(self$selected_features))
  length(self$selected_features)     <- max.length
  length(self$non_selected_features) <- max.length
  self$selected_features[is.na(self$selected_features)] <- ''
  self$non_selected_features[is.na(self$non_selected_features)] <- ''
  knitr::kable(data.frame('Selected Features' = self$selected_features, 'Not Selected Features' = self$non_selected_features), caption = 'Features Selection', format = 'pandoc')
}

rmse <- function(pred, true)
{
  return(sqrt(mean((pred-true)**2)))
}

perf <- function(x, ...) UseMethod('perf')
perf.LassoGLM <- function(self)
{
  rmse <- rmse(pred = self$prediction, true = self$y)
  gini <- axaml::kpi_gini(predrisk = self$prediction, truerisk = self$y )
  (knitr::kable(data.frame('Root Mean Squared error' = rmse, 'Gini index' = as.numeric(gini)), caption = 'Performances of Lasso-GLM model'))
  #invisible(list('Root_mean_squared_error' = rmse, 'Gini_index' = as.numeric(gini)))
}




# si n.sample est NULL, sinon faire du sampling
#Entrer dans computeInteractions pour voir ce qu'il se passe niveau  du pmatch
computeInteractions <- function(gbm.model, data, importance.threshold = 3, n.sample = 5e4, max.select = NULL)
{
  ### Selecting only the most important variables according to a threshold ###
  summary.gbm   <- summary(gbm.model, plotit = F)
  #print(summary.gbm)
  important.var <- summary.gbm[summary.gbm$rel.inf > importance.threshold, "var"]
  if (!is.null(max.select)) {
    important.var <- summary.gbm[1:max.select, "var"]
  }
  #Computing interaction
  index     <- na.omit(pmatch(important.var, colnames(data)))
  n.var     <- length(index)
  mat.inter <- matrix(ncol = n.var, nrow = n.var)
  pb        <- txtProgressBar(0, n.var, style=3)
  best.iter <- gbm::gbm.perf(gbm.model, plot.it = F)

  for (i in 1:n.var) {
    for (j in 1:n.var) {
      setTxtProgressBar(pb, i)

      #print(paste(i, j))
      if (i < j) {
        mat.inter[i, j] = gbm::interact.gbm(x     = gbm.model, data=data[sample(n.sample), ],
                                            i.var = c(colnames(data)[index[i]],
                                                      colnames(data)[index[j]]),
                                            n.trees = best.iter)
      }
    }
  }

  return(mat.inter)
}

getInteractions <- function(mat.inter, interact.threshold = 0, gbm.model)
{
  ### Selection of the highest interactions according to a threshold ###
  mat.high <- which(mat.inter > interact.threshold, arr.ind = T)

  if (nrow(mat.high) == 0) {
    print("There is no interactions at this threshold")
  } else {
    #colnames(gbm.model$data$x.order)[which(mat.inter>0.3, arr.ind=T)]
    mat.name  <- matrix(NA, nrow = nrow(mat.high), ncol = 3)
    for (row in 1:nrow(mat.high)) {
      mat.name[row, 1:2] <- colnames(gbm.model$data$x.order)[mat.high[row, ]]
      mat.name[row, 3]   <- round(mat.inter[ mat.high[row, ][1], mat.high[row, ][2] ], 3)
    }

  }
  mat.name <- matrix(mat.name[order(mat.name[, 3], decreasing = T), ], ncol=3)
  colnames(mat.name) <- c('Variable 1', 'Variable 2', 'H-statistic')

  return(mat.name)
}




glmPerformance <- function(data, gbm.model = gbm.model, formula.best, family,
                    interact.threshold = 0, train.fraction = .8, mat.inter,
                    seed = 1991, speed = F, include = F)
{
  if (speed == T) {
    glm <- speedglm::speedglm
  }
  ### Selection of the highest interactions according to a threshold ###
  mat.high          <- which(mat.inter > interact.threshold, arr.ind = T)

  if (nrow(mat.high) == 0) {
    print("There is no interactions at this threshold")
    return (NULL)
  } else {
    #colnames(gbm.model$data$x.order)[which(mat.inter>0.3, arr.ind=T)]
    mat.name  <- matrix(NA, nrow = nrow(mat.high), ncol = 3)
    for (row in 1:nrow(mat.high)) {
      mat.name[row, 1:2] <- colnames(gbm.model$data$x.order)[mat.high[row, ]]
      mat.name[row, 3]   <- round(mat.inter[ mat.high[row, ][1], mat.high[row, ][2] ], 3)
    }

    mat.name      <- matrix(mat.name[(order(mat.name[,3], decreasing = T)), ], ncol=3)

    if (include) {
      interact.name <- paste(paste(mat.name[, 1], mat.name[, 2], sep = '*'),
                             collapse = "+")
    } else {
      interact.name <- paste(paste(mat.name[, 1], mat.name[, 2], sep = ':'),
                             collapse = "+")
    }

    print(paste0('Interaction(s): ', interact.name, ', number of interactions: ', nrow(mat.name)))

    ### Training the GLMs ###
    set.seed(seed)

    #train <- sample(nrow(data), floor(train.fraction * nrow(data)))
      #TRAIN index with all factros
    data <- na.omit(data)
    dots          <- lapply(X = colnames(data)[sapply(data, is.factor)], FUN = as.symbol)
    data$rownames <- as.numeric(rownames(data))
    dplyr.out <- data %>% dplyr::group_by_(.dots = dots) %>% sample_frac(train.fraction, replace = F)
    train     <- dplyr.out$rownames

    cat("\t Computing Basemodel GLM \n")
    #print(formula.best)
    glm.model <- glm(formula = as.formula(formula.best), data = data[train, ], family = family)

    cat("\t Computing Interaction Model GLM \n")
    formula.plus.interact <- paste(formula.best, interact.name, sep = "+")
    glm.model.plus        <- glm(formula = as.formula(formula.plus.interact), data = data[train, ], family = family)

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

    gini <- axaml::kpi_gini(predrisk = list(gini.base = pred.base, gini.plus = pred.plus), truerisk = true)
    res.gini <- c('Gini Base' = gini[1], 'Gini Plus' = gini[2], 'Diff plus - base' = as.numeric(gini[2] - gini[1]) )

    lr.results    <- lmtest::lrtest(glm.model, glm.model.plus)
    anova.results <- anova(glm.model.plus, test = 'Chisq')
    print(options)
    print(res.rmse)
    print(res.gini)
    pred <- list('base' = pred.base, 'plus' = pred.plus)
    plot_gain.control <- list(predrisk = pred, truerisk = true)
    return(list(RMSE = res.rmse, Gini = res.gini, options = options, likelihood.ratio.test = lr.results, anova.results = anova.results,
                formula.plus = formula.plus.interact, predictions = pred, plot_gain.control = plot_gain.control, mat = mat.name))
  }
}
# mettre SPEED GLM
# Mettre indicateur GINI dans les if else. si package pour plotter/calculer gini facilement, le prendre. Difference simple en terme de GINI
#faire juse GINI et RMSE à chaque fois. SI count.
# Voir pour le calcul de la deviance hors-ech.

getGlmPerformance <- function(threshold.values = c(0.01, seq(0.05, 1, by = 0.05)), data, mat.inter, gbm.model, formula.best,
                              family, train.fraction = .8, max.interactions = 50, include = F)
{
  min.interaction.value  <- min(mat.inter, na.rm = T)
  max.interaction.value  <- max(mat.inter, na.rm = T)
  threshold.values       <- threshold.values
  first.threshold        <- which(min.interaction.value <= threshold.values)[1]
  last.threshold         <- which(max.interaction.value <= threshold.values)[1]
  stored.results         <- list()
  length(stored.results) <- 20

  #Loop through all the values of the threshold

  for (i in first.threshold:last.threshold) {
    print(paste('H-statistic threshold value = ', threshold.values[i]))

    if (sum(mat.inter > threshold.values[i], na.rm = T) < max.interactions) {
      if ( (sum(mat.inter > threshold.values[i-1], na.rm = T) > sum(mat.inter > threshold.values[i], na.rm = T) ) & i > 1) {
        stored.results[[i]] <- glmPerformance(interact.threshold = threshold.values[i], train.fraction = train.fraction,
                                              gbm.model = gbm.model, formula.best = formula.best, mat.inter = mat.inter, data = data,
                                              family = family, include = include)
    } else {
      print('Same interactions as previous iteration, going to next iteration')
    }
      } else {
        print(paste("Number of interactions = ", sum(mat.inter > threshold.values[i], na.rm = T), '>', max.interactions, "= maximum number of interaction"))
      }
    cat('\n')
    }

  return(stored.results)
}


detectInteractions <- function(data, max.interactions, formula.best,
                               gbm.control = list(train.fraction = .8, family, shrinkage = .01, bag.fraction = .7, n.trees = 100, n.sample = 5e4, depth = 5, importance.threshold = 0, max.select = 30),
                               glm.control = list(train.fraction = .8, family, include = F, threshold.values = c(0.01, seq(0.05, 1, by = 0.05))))
{

  print("Fitting GBM Model")
  #Fitting GBM model
  gbm.model <- gbm::gbm(y ~ . , distribution = gbm.control$family, data = data, n.trees = gbm.control$n.trees, interaction.depth = gbm.control$depth,
                        shrinkage = gbm.control$shrinkage, bag.fraction = gbm.control$bag.fraction, verbose = T, train.fraction = gbm.control$train.fraction)

  print("Computing interaction matrix")
  #computing friedman H matrix
  mat.inter <- computeInteractions(gbm.model, data = data, importance.threshold = gbm.control$importance.threshold, n.sample = gbm.control$n.sample, max.select = gbm.control$max.select)
  cat('\n')

  print(getInteractions(mat.inter = mat.inter, gbm.model = gbm.model, interact.threshold = 0))

  print('Fitting GLMs model')
  results <- getGlmPerformance(threshold.values = glm.control$threshold.values , data = data, mat.inter = mat.inter,
                               gbm.model = gbm.model, formula.best = formula.best, family = glm.control$family,
                               train.fraction = glm.control$train.fraction, max.interactions = max.interactions, include = glm.control$include) #0-1 la famille doit être binomial

  function.call <- as.list(sys.call())
  output        <- list(GLMs.perf = results, best.interactions = getInteractions(mat.inter = mat.inter, gbm.model = gbm.model, interact.threshold = 0),
                        function.call = function.call)
  class(output) <- 'InteractionsDetection'
  invisible(output)

}

#modifier la fonction pour faire apparaître bernoulli ou binomial


#definition of methods for InteractionsDetection object

summary.InteractionsDetection <- function(self, criterion = 'Gini')
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
                                                     sep = ''), format = 'pandoc')
}


perf <- function(x, ...) UseMethod('perf')
perf.InteractionsDetection <- function(self, criterion = 'Gini')
{
  index <- rep(NA, length(self$GLMs.perf))
  for (i in 1:length(index)) {
    if (!is.null(self$GLMs.perf[[i]][[criterion]][3])) {
      index[i] <- self$GLMs.perf[[i]][[criterion]][3]
    }
  }
  index.best <- which.max(index)
  best.model <- self$GLMs.perf[[index.best]]
  dataframe <- data.frame('Root Mean Squared error' = best.model$RMSE, 'Gini index' = best.model$Gini)
  rownames(dataframe) <- c('Base Model', 'Augmented Model', 'Difference augmented - base ')
  cat1 <- paste0("Performances of augmented model vs base model according to '", criterion, "' criterion,", sep = '')
  knitr::kable(dataframe, caption = paste(cat1,paste0('Train fraction:', best.model$options[3]), sep = '\n'),
                     format = 'pandoc')
  #invisible(list('Root_mean_squared_error' = best.model$RMSE, 'Gini_index' = as.numeric(best.model$Gini)))
}


#    return(list(RMSE = res.rmse, Gini = res.gini, likelihood.ratio.test = lr.results, anova.results = anova.results, formula.plus = formula.plus.interact, predictions = pred, plot_gain.control = plot_gain.control))


#plot.myclass <- function(self, ...) UseMethod('plot') -- not needed for already defined methods like plot, summary
plot.InteractionsDetection <- function(self, i = 1, y, data = eval(self$function.call$data),
                          rev = NULL, cut.x = NULL, cut.trace = NULL, simplify = NULL)
{
  # This method plot the i-th interactions, by default it plot the best interaction

  plotInteract(self$best.interactions[i, 1:3], y=y, data=data, rev=rev, cut.x=cut.x, cut.trace=cut.trace, simplify = simplify)
}

# Ne pas oublier les TROIS PETIT POINTS si on utilise plus d'un paramètre
bestInteractions <- function(x, ...) UseMethod('bestInteractions')

bestInteractions.InteractionsDetection <- function(self, i = 5) {
  return(knitr::kable(self$best.interactions[1:i, ], format = 'pandoc', caption = 'Best interactions'))
}

plotBest <- function(x, ...) UseMethod('plotBest')
plotBest.InteractionsDetection <- function(self, criterion = 'Gini') {
  index <- rep(NA, length(self$GLMs.perf))
  for (i in 1:length(index)) {
    if (!is.null(self$GLMs.perf[[i]][[criterion]][3])) {
      index[i] <- self$GLMs.perf[[i]][[criterion]][3]
    }
  }
  index.best <- which.max(index)
  do.call(what = axaml::plotgain, args = c(self$GLMs.perf[[index.best]]$plot_gain.control, return_val = T))$plot
}

interactionsBestModelPlot <- function(x, ...) UseMethod('interactionsBestModelPlot')
interactionsBestModelPlot.InteractionsDetection <- function(self, criterion = 'Gini') {
  index <- rep(NA, length(self$GLMs.perf))
  for (i in 1:length(index)) {
    if (!is.null(self$GLMs.perf[[i]][[criterion]][3])) {
      index[i] <- self$GLMs.perf[[i]][[criterion]][3]
    }
  }
  index.best <- which.max(index)
  best.model <- self$GLMs.perf[[index.best]]
  for(i in 1:best.model$options[3]){
    plot(self, i = i, simplify= 8)
  }
}




#function plot interact
plotInteract <- function(x, y, data, rev = NULL, cut.x = NULL, cut.trace = NULL, simplify = NULL)
{
  if (is.null(rev)) {
    if(is.factor(data[, x[1]])){
      rev = T
    }
  } else {
    if (rev) {
      x[1:2] <- rev(x[1:2])
    }
  }


  xfactor     <- data[, x[1]]
  tracefactor <- data[, x[2]]

  if (!is.null(cut.x)) {
    xfactor <- Hmisc::cut2(x = xfactor, g = cut.x)
  }
  if (!is.null(cut.trace)) {
    tracefactor <- Hmisc::cut2(x = tracefactor, g = cut.trace)
  }
  if (!is.null(simplify)) {
    if(!is.factor(data[, x[1]])) {
      xfactor <- Hmisc::cut2(x = xfactor, g = simplify)
    }
    if(!is.factor(data[, x[2]])) {
      tracefactor <- Hmisc::cut2(x = tracefactor, g = simplify)
    }
  }
  interaction.plot(x.factor = xfactor, trace.factor = tracefactor,
                   response = data$y, xlab = x[1], ylab = "y", trace.label = x[2],
                   col = 1:length(levels(data[, x[2]])), main = paste('Interaction :', paste(x[1:2], collapse = '*'),
                                                                      '\n H-statistic :', x[3]))
}


##################   TRASH   ################################################


#
# automaticDetection <- function(data, family, max.interactions = 50,
#                                 gbm.control = list(train.fraction = .8, shrinkage = .01, bag.fraction = .7, n.trees = 100, n.sample = 5e4, depth = 5, importance.threshold = 0),
#                                 lasso.control = list(offset, lambda = "lambda.min", type.measure, train.fraction = .8),
#                                 glm.control = list(train.fraction = .8, include = F,  threshold.values = c(0.01, seq(0.05, 1, by = 0.05))))
# {
#   if (deparse(substitute(family)) == "binomial") {
#     family.gbm <- 'bernoulli'
#   } else if (deparse(substitute(family)) == "gaussian") {
#     family.gbm <- 'gaussian'
#   } else if (deparse(substitute(family)) == "gamma") {
#     family.gbm <- 'gamma'
#   }
#
#
#   print("Fitting GBM Model")
#   #Fitting GBM model
#   gbm.model <- gbm::gbm(y ~ . , distribution = family.gbm, data = data, n.trees = gbm.control$n.trees, interaction.depth = gbm.control$depth,
#                         shrinkage = gbm.control$shrinkage, bag.fraction = gbm.control$bag.fraction, verbose = T, train.fraction = gbm.control$train.fraction)
#
#   print("Computing interaction matrix")
#   #computing friedman H matrix
#   mat.inter <- computeInteractions(gbm.model = gbm.model, data = data, importance.threshold = gbm.control$importance, n.sample = gbm.control$n.sample)
#   #(gbm.model = gbm.model, data = data, importance.threshold = gbm.control$importance, n.sample = )
#   cat('\n')
#
#   print('Finding best model using LASSO')
#   #Finding the best base GLM model
#   best.glm <- bestLasso2(data = data, type.measure = lasso.control$type.measure, family = deparse(substitute(family)), offset = lasso.control$offset,
#                         lambda = lasso.control$lambda, train.fraction = lasso.control$train.fraction) # 0-1 family doit être "binomial"
#
#   print('Fitting GLMs model')
#   results <- getGlmPerformance(threshold.values = glm.control$threshold.values , data = data, mat.inter = mat.inter,
#                           gbm.model = gbm.model, formula.best = best.glm$formula, family = (family),
#                           train.fraction = glm.control$train.fraction, max.interactions = max.interactions, include = glm.control$include) #0-1 la famille doit être binomial
#
#   function.call <- as.list(sys.call())
#   output        <- list(GLMs.perf = results, best.interactions =  getInteractions(mat.inter, gbm.model = gbm.model), function.call = function.call)
#   class(output) <- 'InteractionsDetection'
#   return(output)
# }
#
#



#
# bestLasso <- function(data, type.measure = "mse", family = "gaussian", offset, lambda = "lambda.min", train.fraction = .8){
#
#   Xy    <- model.matrix( ~ . - 1 , data = na.omit(data))
#   train <- sample(nrow(data), floor(train.fraction * nrow(data)))
#
#   if (is.null(offset)) {
#
#     ### Launching Cross validation LASSO ###
#     cvfit      <- glmnet::cv.glmnet (x = subset(Xy[train, ], select = -c(y)),
#                                     y = Xy[train, "y"], type.measure = type.measure, family = family)
#     plot(cvfit)
#     best.model <- coef(cvfit, s = lambda)
#
#       #retirer l'exposure des var expl quand offset == NULL
#
#     ### string work ###
#     i.best               <- best.model@i
#     matching.index       <- na.omit(pmatch(best.model@Dimnames[[1]][i.best+1], colnames(data)))
#     best.features        <- colnames(data)[matching.index]
#     concatenate.features <- paste(best.features, collapse = '+')
#
#   } else {
#
#     ### Launching Cross validation LASSO ###
#     cvfit     <- glmnet::cv.glmnet(x = subset(Xy[train, ], select = -c(offset, y)),
#                                    y = Xy[train, "y"], type.measure = type.measure, family = family, offset = Xy[train, offset])
#     plot(cvfit)
#     best.model <- coef(cvfit, s = lambda)
#
#     ### string work ###
#     i.best               <- best.model@i
#     matching.index       <- pmatch(best.model@Dimnames[[1]][i.best + 1], colnames(data))
#     best.features        <- colnames(data)[na.omit(matching.index)]
#     concatenate.features <- paste(paste(best.features, collapse = '+'), paste0("offset(", offset, ")"), sep = "+")
#
#   }
#
#   formula.best <- paste('y', concatenate.features, sep = '~')
#
#   cat(paste0("\t", formula.best, "\n"))
#
#
#   return(list(formula = formula.best, cvfit = cvfit))
# }




#
#
# bestLasso2 <- function(data, type.measure = "mse", family = "gaussian", offset, lambda = "lambda.min",
#                        train.fraction = .8, seq.lambda = NULL){
#   data <- na.omit(data)
#   Xy <- model.matrix( ~ . - 1 , data = data)
#   train <- sample(nrow(data), floor(train.fraction * nrow(data)))
#
#   if(family == 'gamma'){
#     cv.lasso   <- HDtweedie::cv.HDtweedie
#     cv.control <- list(pred.loss = type.measure)
#   } else {
#     cv.lasso   <- glmnet::cv.glmnet
#     cv.control <- list(type.measure = type.measure, family=family)
#   }
#
#   if(! is.null(seq.lambda)) {
#     cv.control <- c(cv.control, lambda = list(seq.lambda))
#   }
#   if (is.null(offset)) {
#
#     ### Launching Cross validation LASSO ###
#     cvfit      <- do.call(cv.lasso, c(list(x = subset(Xy[train, ], select = -c(y)),
#                                            y = Xy[train, "y"]), cv.control))
#     #plot(cvfit)
#     best.model <- coef(cvfit, s = lambda)
#
#     #retirer l'exposure des var expl quand offset == NULL
#
#     ### string work ###
#     if (family == 'gamma'){
#       i.best <- which(best.model[, 1]!=0)
#       matching.index <- na.omit(pmatch(names(i.best), colnames(data)))
#     } else {
#       i.best               <- best.model@i
#       matching.index       <- na.omit(pmatch(best.model@Dimnames[[1]][i.best + 1], colnames(data)))
#     }
#
#     best.features        <- colnames(data)[matching.index]
#     paste0('#features', length(best.features))
#     if (length(best.features) < 1) {
#       stop('Lasso has not selected any feature')
#     }
#     cat(paste0("\t", best.features, "\n"))
#     concatenate.features <- paste(best.features, collapse = '+')
#
#   } else {
#
#     ### Launching Cross validation LASSO ###
#     cvfit     <- do.call(cv.lasso, c(list(x = subset(Xy[train, ], select = -c(get(offset), y)),
#                                           y = Xy[train, "y"], offset = Xy[train, offset]), cv.control))
#
#     plot(cvfit)
#     best.model <- coef(cvfit, s = lambda)
#
#     ### string work ###
#     if (family == 'gamma'){
#       i.best <- which(best.model[, 1]!=0)
#       matching.index <- na.omit(names(i.best), colnames(data))
#     } else {
#       i.best               <- best.model@i
#       matching.index       <- na.omit(pmatch(best.model@Dimnames[[1]][i.best + 1], colnames(data)))
#     }
#     best.features        <- colnames(data)[na.omit(matching.index)]
#
#     paste0('#features', length(best.features))
#
#     if (length(best.features) < 1) {
#       stop('Lasso has not selected any feature')
#     }
#     cat(paste0("\t", best.features, "\n"))
#     concatenate.features <- paste(paste(best.features, collapse = '+'), paste0("offset(", offset, ")"), sep = "+")
#   }
#   plot(cvfit)
#   formula.best <- paste('y', concatenate.features, sep = '~')
#
#   return(list(formula = formula.best, cvfit = cvfit))
# }
