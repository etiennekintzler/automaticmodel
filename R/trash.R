
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