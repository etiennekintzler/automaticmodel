# faire une vignette avec portfolio example sur le README.
# dans le package folder - 'doc' , mettre le markdown dont avc le regroupement des interactions
# carrément inclure la data de portfolio_example ; regarder wickham rpackage
# mettre un target = 'name'
# permettre le log dans le string pour offset. Dire de mettre le vecteur dans les paramètres

plotInteractionsMatrix <- function(mat.inter)
{
  rgb.palette <- colorRampPalette(c("white", "blue"), space = "rgb")
  lattice::levelplot(matCost, col.regions = rgb.palette(40), main = "Interaction visualization")
}


#' Taking sample for each level of factors in a dataset
#'
#' Returns the row numbers of a training sample of a dataset making sure to draw in each level of the dataset
#'
#' @param data The dataset
#' @param train.fraction the fraction of the observations to be drawn in each level of the factors of the dataset
#' @seealso \code{\link[dplyr]{group_by}} and \code{\link[dplyr]{sample_frac}} in \pkg{dplyr} package
#' @export
trainingSample <- function(data, train.fraction = 0.75){
  dots         <- colnames(data)[sapply(data, is.factor)]
  data$row.num <- 1:nrow(data)
  dplyr.out    <- dplyr::sample_frac(dplyr::group_by_(data, .dots = dots), train.fraction, replace = F)
  train        <- dplyr.out$row.num
  return(train)
}


#function plot interact
#' plot interaction
#'
#' @export
plotInteract <- function(x, y, data, rev = NULL, cut.x = NULL, cut.trace = NULL,
                         simplify = NULL)
{
  x <- c(as.character(x[1, 1]), as.character(x[1, 2]), x[1, 3])
  if (!is.null(rev)) {
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
    if (length(levels(xfactor)) > length(levels(tracefactor))) {
      temp <- xfactor
      xfactor <- tracefactor
      tracefactor <- temp
      temp <- x[1]
      x[1] <- x[2]
      x[2] <- temp
    }
  }
  interaction.plot(x.factor = xfactor, trace.factor = tracefactor,
                   response = y, xlab = x[1], ylab = "y", trace.label = x[2],
                   col = 1:length(unique(tracefactor)), fixed = T,
                   main = paste('Interaction :', paste(x[1:2], collapse = '*'),
                                '\n H-statistic :', round(as.numeric(x[3]), 4)))
}

#' Root mean squared error
#'
#' Compute the RMSE
#'
#' @details
#' The root mean squared error is a an accuracy criterion which writes as follow:
#' \deqn{RMSE = \sqrt{\frac{1}{N}\sum_{i=1}^N ( y_i - \hat{y}_i)^2  } }
#' @export
rmse <- function(pred, true)
{
  return(sqrt(mean((pred-true)**2)))
}

#' Cleaning the data
#'
#' Removing variables with too much levels, too much NA, no variance, too much correlation.
#' @param data the dataset
#' @param n.levels the maximum number of levels allowed
#' @param per.na the maximum percentage of missing data allowed for a variable
#' @param na.string the name of the \code{NA} values in the data
#' @param remove.cor used to remove variables which correlation is superior to this value
#' @return A clean and neat dataset
#' @export
cleanData <- function(data, n.levels = 20, perc.na = 0.2, na.string = NULL, remove.cor = NULL, target = NULL)
{
  if (!is.null(remove.cor)){
    if (is.null(target)){
      "'target' name must be provided is 'remove.cor' is not 'NULL'"
    }
  }
  data <- data.frame(data)
  high.levels    <- sapply(data, function(x) length(levels(x)) > n.levels) #justifier dans rapport
  na.columns     <- sapply(data, function(x) sum(is.na(x)) > perc.na * nrow(data))
  date.col       <- grepl(pattern = 'date', x = colnames(data))
  keep.columns   <- !(high.levels | na.columns | date.col)
  data.temp      <- data[, keep.columns]
  if (!is.null(na.string)){
    data.temp[data.temp == na.string] <- NA
  }
  data.temp <- na.omit(data.temp)
  if (!is.null(remove.cor)) {
    X                   <- subset(data.temp, select = -y)
    tmp                 <- cor(X[, sapply(X, is.numeric)])
    tmp[upper.tri(tmp)] <- 0
    diag(tmp)           <- 0
    X.new               <- X[, ! apply(tmp, 2,function(x) any(x > remove.cor))]
    data.temp           <- cbind(X.new, y = data.temp$y)
  }
  data.final <- as.data.frame(lapply(data.temp, function(x) if(is.factor(x)) factor(x) else x))
  novar.col  <- sapply(data.final, function(x) length(unique(x)) <= 1)
  return(na.omit(data.final[, ! novar.col]))
}

