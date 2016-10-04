test_that("", {
  library(axamodel)
  library(dplyr)

  data(diabetes, package = 'lars')
  data <- cbind(data.frame(unclass(diabetes$x)), y = diabetes$y)
  df <- cbind(subset(x = data, select = -y),
              y = data$y)

  data <- axaml::portfolio_example
  df <- cbind(subset(x = data, select = -c(CalYear, Exposure, Numtpbi, Indtppd, Indtpbi, Numtppd, PolNum)),
              y = data$Numtpbi) %>% cleanData(n.levels = 10)

  n = 1e4

  best <- lassoSelection(data = df[1:n, ], train.fraction = 1, family = 'poisson')
  plot(best)


  debugonce(detectInteractions)
  library(axamodel)

set.seed(124)
interactions <- detectInteractions(data = df[1:n, ], max.interactions = 20,
                                   formula.best = best$formula, shuffle = F,
                                   gbm.control = list(train.fraction = 0.75, family = 'poisson', shrinkage = 0.01,
                                                      bag.fraction = .5, n.trees = 2000, n.sample = n/5, depth = 5,
                                                      importance.threshold = 0, max.select = 100),
                                   glm.control = list(train.fraction = 0.75, family = poisson(link = 'log'), include = F,
                                                      threshold.values = seq(0, 1, length.out = 50)))

debugonce(plot)
plot(interactions, 1, rev = T)
plotAll(interactions)

interactions$best.interactions
Variable 1     Variable 2 H-statistic
1      Occupation         Region 0.364694663
2      Occupation      Value_Car 0.340063550
3  Density_Region     Occupation 0.302858277
4      Occupation       Type_Car 0.275163675
5          Region       Type_Car 0.188462309
6      Occupation         Gender 0.171255955
7      Occupation   Category_Car 0.168797644

Occupation_modified <- df$Occupation
#levels(Occupation_modified)  <- list(Active = c('Employed','Housewife','Self-employed'), Retired = 'Retired', Unemployed = 'Unemployed')
levels(Occupation_modified) <- c('active', 'active','retired' ,'active', 'unemployed')
df.modified <- df
df.modified$Occupation  <- Occupation_modified

Value_Car_modified <- Hmisc::cut2(x = df$Value_Car, g = 25)
levels(Value_Car_modified) <- c(rep('low_to_high_val', 24), 'really_high')
df.modified <- df
df.modified$Value_Car <- Value_Car_modified

set.seed(124)
a = detectInteractions(df.modified[1:n, ], max.interactions = 20,
                   formula.best = best$formula, shuffle = F,
                   gbm.control = list(train.fraction = 0.75, family = 'poisson', shrinkage = 0.01,
                                      bag.fraction = .5, n.trees = 2000, n.sample = n/5, depth = 5,
                                      importance.threshold = 0, max.select = 100),
                   glm.control = list(train.fraction = 0.75, family = poisson(link = 'log'), include = F,
                                      threshold.values = seq(0, 1, length.out = 50)))

a$best.interactions


debugonce(perf)
perf(interactions)

 set.seed(124)

 microbenchmark::microbenchmark(interactions <- detectInteractions(data = df[1:n, ], max.interactions = 20,
                                                                   formula.best = best$formula, shuffle = F,
                                                                   gbm.control = list(train.fraction = 0.75, family = 'poisson', shrinkage = 0.01,
                                                                                      bag.fraction = .5, n.trees = 200, n.sample = n/10, depth = 5,
                                                                                      importance.threshold = 1, max.select = 100),
                                                                   glm.control = list(train.fraction = 0.25, family = poisson(link = 'log'), include = F,
                                                                                      threshold.values = seq(0, 1, length.out = 50))),
                                times = 1L)


 debugonce(plot)
 plot(interactions, 1, simplify = 10)



  set.seed(124)
  m <- gbm::gbm(formula = y~ . , train.fraction = .75, data = df,  distribution = 'poisson',
                shrinkage = .01, bag.fraction = .5, n.trees = 500, interaction.depth = 5,
                verbose = T)
  gbm::gbm.perf(m)

  p <- computeInteractions2(gbm.model = m, data = df, importance.threshold = 1, n.sample = n/20, max.select = 100, interact.threshold = 0)
  p
  interactionsBestModelPlot(interactions, simplify = 5)


})
