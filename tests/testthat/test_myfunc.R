test_that("", {
  library(axamodel)
  library(magrittr)
  data(diabetes, package = 'lars')
  data <- cbind(data.frame(unclass(diabetes$x)), y = diabetes$y)
  df <- cbind(subset(x = data, select = -y),
              y = data$y)

  data <- axaml::portfolio_example
  df <- cbind(subset(x = data, select = -c(CalYear, Exposure, Numtpbi, Indtppd, Indtpbi, Numtppd, PolNum)),
              y = data$Numtpbi) %>% cleanData(n.levels = 10)

  n = nrow(df)

  best <- lassoSelection(data = df[1:n, ], train.fraction = 1, family = 'poisson')
  plot(best)

  debugonce(detectInteractions)
  set.seed(124)
  interactions <- axamodel::detectInteractions(data = df[1:n, ], max.interactions = 20,
                                               formula.best = best$formula, shuffle = F,
                                               gbm.control = list(train.fraction = 0.75, family = 'poisson', shrinkage = 0.01,
                                                                  bag.fraction = .5, n.trees = 200, n.sample = n/10, depth = 5,
                                                                  importance.threshold = 1, max.select = 100),
                                               glm.control = list(train.fraction = 0.25, family = poisson(link = 'log'), include = F,
                                                                  threshold.values = seq(0, 1, length.out = 50)))

  set.seed(124)
  m <- gbm::gbm(formula = y~ . , train.fraction = .75, data = df,  distribution = 'poisson',
                shrinkage = .01, bag.fraction = .5, n.trees = 500, interaction.depth = 5,
                verbose = T)
  gbm::gbm.perf(m)
  p <- computeInteractions2(gbm.model = m, data = df, importance.threshold = 1, n.sample = n/20, max.select = 100, interact.threshold = 0)
  p

  gbm::interact.gbm(x = m, data = df , i.var = c('Region', 'Value_Car')  , n.trees = gbm::gbm.perf(m, plot.it = F))
  perf(interactions)
  summary(interactions)


  interactionsBestModelPlot(interactions, simplify = 5)


  interactions$best.interactions[1:10,]
  #anova(glm(y ~ .*. , data  = data), test = 'Rao')




})
