data <- mtcars
df <- cbind(subset(x = data, select = -mpg),
            y = data$mpg)
#debugonce(lassoSelection)
best.glm <- lassoSelection(data = df[sample(x = nrow(df), size = 1e3, replace=T),], family = 'gaussian', train.fraction = .75)
perf(best.glm)


#debugonce(detectInteractions)

interactions <- detectInteractions(data = df[sample(x = nrow(df), size = 1e4, replace=T),], max.interactions = 10, formula.best = best.glm$formula,
                                   gbm.control = list(train.fraction = .5, family = 'gamma', shrinkage = .1, bag.fraction = .7,
                                                      n.trees = 100, n.sample = 1e3, depth = 5, importance.threshold = 0, max.select = 30),
                                   glm.control = list(train.fraction = .8, family  = gaussian(), include = F,
                                                      threshold.values = seq(0.01, 1, by = 0.02), speed = F))

perf(interactions)
plot(interactions)
#plot(interactions, 9, simplify = 8)
#interactionsBestModelPlot(interactions)


