test_that("", {
library(axamodel)
library(dplyr)

data(diabetes, package = 'lars')
data <- data.frame(unclass(diabetes$x), y = diabetes$y)
n    <- nrow(data)


data <- axaml::portfolio_example
df <- cbind(subset(x = data, select = -c(CalYear, Indtppd, Indtpbi, Numtppd, PolNum))) %>%
  axamodel::cleanData(n.levels = 10)


best <- lassoSelection(data = data, target ='y' ,family = 'poisson', nfolds = 3, nlambda = 20)
plotLambda(best)

plot(best)

set.seed(124)

interactions <- axamodel::detectInteractions(data = data, max.interactions = 20,
                                   formula.best = best$formula, shuffle = F, target = 'y',
                                   gbm.control = list(train.fraction = 0.75, distribution = 'poisson', shrinkage = 0.01,
                                                      bag.fraction = .5, n.trees = 20, n.sample = n/5,
                                                      interaction.depth = 5,
                                                      importance.threshold = 0, max.select = 100),
                                   glm.control = list(train.fraction = 0.75, family = poisson(link = 'log'), include = F,
                                                      threshold.values = seq(0, 1, length.out = 50)))

plot(interactions)
})
