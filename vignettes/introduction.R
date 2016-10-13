## ---- eval = F-----------------------------------------------------------
#  devtools::install_git('https://github.com/etiennekintzler/automaticmodel.git')

## ---- eval = F-----------------------------------------------------------
#  devtools::install_git('https://github.com/etiennekinthttps://github.axa.com/ds4a/axamodel.git')

## ---- include = F--------------------------------------------------------
knitr::opts_chunk$set(cache = TRUE, message = FALSE)

## ---- include = T, echo = F, results = 'hide', message = F---------------
library(axamodel)
library(dplyr)

## ------------------------------------------------------------------------
data <- axaml::portfolio_example
df   <- subset(x = data, select = -c(CalYear, Indtppd, Indtpbi, Numtppd, PolNum)) %>% 
  cleanData(n.levels = 10)

## ---- echo = F-----------------------------------------------------------
hist(df$Numtpbi, main = 'Histogram of Numtpbi', xlab = 'Numtpbi')

