#function plot interact
plotInteract <- function(x, y, data, rev = NULL, cut.x = NULL, cut.trace = NULL,
                         simplify = NULL)
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

#' Root mean squared error
#'
#' Compute the RMSE
#'
#' @details
#' The root mean squared error is a an accuracy criterion which writes as follow:
#' \deqn{RMSE = \sqrt{\frac{1}{N}\sum_{i=1}^N ( y_i - \hat{y}_i)^2  } }
rmse <- function(pred, true)
{
  return(sqrt(mean((pred-true)**2)))
}

#' Cleaning the data
#'
#' Removing variables with too much levels, too much NA, too few variance
#' @param data the dataset
#' @param n.levels the maximum number of levels allowed
#' @param per.na the maximum percentage of missing data allowed for a variable
#' @return A clean and neat dataset
#' @export
cleanData <- function(data, n.levels = 20, perc.na = 0.2)
{
  novar.col      <- sapply(data, function(x) length(unique(x)) <= 1)
  high.levels    <- sapply(data, function(x) length(levels(x)) > n.levels) #justifier dans rapport
  na.columns     <- sapply(data, function(x) sum(is.na(x)) > perc.na * nrow(data))
  date.col       <- grepl(pattern = 'date', x = colnames(data))
  keep.columns   <- !(high.levels | na.columns | novar.col | date.col)
  return(data[, keep.columns])
}
