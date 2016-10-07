test_that("", {
library(axamodel)
library(dplyr)

data(diabetes, package = 'lars')
data <- cbind(data.frame(unclass(diabetes$x)), y = diabetes$y)
df <- cbind(subset(x = data, select = -y),
            y = data$y)
n <- nrow(df)

#
data <- axaml::portfolio_example
df <- cbind(subset(x = data, select = -c(CalYear, Indtppd, Indtpbi, Numtppd, PolNum))) %>%
  axamodel::cleanData(n.levels = 10)

n = 5e4
debugonce(lassoSelection)

best <- lassoSelection(data = df[1:n, ], target ='Numtpbi' ,
                       train.fraction = 1, family = 'poisson', nlambda = 50)
plotLambda(best)

lassoSelection()
plot(best)

set.seed(124)
debugonce(axamodel::detectInteractions)
interactions <- axamodel::detectInteractions(data = df[1:n, ], max.interactions = 20,
                                   formula.best = best$formula, shuffle = F,
                                   gbm.control = list(train.fraction = 0.75, family = 'poisson', shrinkage = 0.01,
                                                      bag.fraction = .5, n.trees = 20, n.sample = n/5, depth = 5,
                                                      importance.threshold = 0, max.select = 100),
                                   glm.control = list(train.fraction = 0.75, family = poisson(link = 'log'), include = F,
                                                      threshold.values = seq(0, 1, length.out = 50)))


d <- data.frame(y = (x1 <- rnorm(1000))+ (x2 <- rpois(1000, 4)), x1 , x2)
set.seed(1)
microbenchmark::microbenchmark(m <- gbm::gbm(y ~ . , data = df, interaction.depth = 5, verbose = T,
                                             n.trees =1000, cv.folds = 1, train.fraction = .75, shrinkage = .1),
                               times = 1L)
set.seed(1)

microbenchmark::microbenchmark(mm <- gbmt(y ~ . , data = df,
           train_params = training_params(interaction_depth = 5, num_trees =500,
                                          bag_fraction = .5, num_train = floor(.75*nrow(df)),
                                          id = seq_len(nrow(df)), shrinkage = .1),
           is_verbose = T, cv_folds = 1, par_details = gbmParallel(num_threads = 1)), times = 1L)

gbm.perf(m)
plot(gbmt_performance(m, method = 'OOB'))
plot(gbmt_performance(m, method = 'cv')); gbmt_performance(m, method = 'cv')

plot(gbmt_performance(mm, method = 'OOB'))


plot(gbmt_performance(mm, method = 'OOB')); plot(gbmt_performance(mm, method = 'cv'))

o$variables$var_names
debug(interact)
o <- mm
interact(gbm_fit_obj = o, data = df ,var_indices = c(5, 7 ), num_trees = gbmt_performance(o))
undebug(interact)


anova(glm(y ~ .*., data = df[, (ncol(df)-6):ncol(df)]), test = 'Chisq')
(train <- sample(nrow(df), floor(.8*nrow(df))) )
g <- glm(y ~ . , data = df[train, ])
axaml::plotgain(truerisk = df$y[-train],
                 predrisk = list(base = predict(g, newdata= df[-train, ], type = 'response' ),
                                 plus = predict(g2, newdata= df[-train, ], type = 'response' )))

g2 <- glm(y ~ . + Occupation:Age + Age*Adind , data = df[train, ])

g3 <-  glm(y ~ .*. , data = df[train, ])
summary(g3)

returnBest(interactions)
plotBest(interactions)
summary(interactions)
plotAll(interactions)
interactions$best.interactions

plot(interactions, simplify = 25)

var_indices = c('x1', 'x2')
(
!is.atomic(var_indices) || any(is.infinite(var_indices)) ||
  any(is.na(var_indices)) || any(is.nan(var_indices)) ||
  !(all(var_indices == as.integer(var_indices)) ||
      all(var_indices == as.character(var_indices))) ||
  is.na(all(var_indices == as.integer(var_indices))) ||
  is.na(all(var_indices == as.character(var_indices)))
)

# Variable 1     Variable 2 H-statistic
# 1      Occupation         Region 0.364694663
# 2      Occupation      Value_Car 0.340063550
# 3  Density_Region     Occupation 0.302858277
# 4      Occupation       Type_Car 0.275163675
# 5          Region       Type_Car 0.188462309
# 6      Occupation         Gender 0.171255955
# 7      Occupation   Category_Car 0.168797644

Occupation_modified <- df$Occupation
#levels(Occupation_modified)  <- list(Active = c('Employed','Housewife','Self-employed'), Retired = 'Retired', Unemployed = 'Unemployed')
levels(Occupation_modified) <- c('active', 'active','retired' ,'active', 'unemployed')
df.modified <- df
df.modified$Occupation  <- Occupation_modified

Value_Car_modified <- Hmisc::cut2(x = df$Value_Car, g = 25)
levels(Value_Car_modified) <- c(rep('low_to_high_val', 23),'really_high', 'highest')
df.modified <- df
df.modified$Value_Car <- Value_Car_modified

set.seed(124)
a = detectInteractions(df.modified[1:nrow(df), ], max.interactions = 20,
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
              shrinkage = .01, bag.fraction = 1, n.trees = 500, interaction.depth = 5,
              verbose = T)
gbm::gbm.perf(m)

p <- computeInteractions(gbm.model = m, data = df, importance.threshold = 1, n.sample = n/20, max.select = 100, interact.threshold = 0)
p

})
