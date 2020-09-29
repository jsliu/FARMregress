#' farm regression
#' @param x predictors
#' @param y response
#' @import glmnet
#' @export
farmRegress <- function(x, y, robust=T,
                        type=c("regression","classification"), ...)
{
  family <- switch(type,
                   "regression"="gaussian",
                   "classification"="binomial")

  if (!is.matrix(x)) x <- as.matrix(x)

  xt <- t(x)


  P <- NROW(xt)
  N <- NCOL(xt)

  if (robust) {
    if (type == "regression") {
      # y has to be centralized in linear regression
      my <- meanHuber(y)
      y <- y-my
    }
    mx <- meanHuber(x)
    xt <- sweep(xt,1,mx)
    covx <- covHuber(x)
    eigen_fit <- eigen(covx)
    values <- eigen_fit$values
    vectors <- eigen_fit$vectors

  } else {
    if (type == "regression") {
      y <- y-mean(y,na.rm=T)
    }
    xt <- sweep(xt,1,rowMeans(xt))
    covx <- cov(t(x))
    pca_fit <- stats::prcomp(xt, center = TRUE, scale = TRUE)
    values <-  (pca_fit$sdev)^2
    vectors <- pca_fit$rotation
  }

  # choose number of factors
  values <- pmax(values,0)
  ratio <- c()
  K.factors <- NULL
  K.factors <- if (is.null(K.factors)) {
    for(i in 1:(floor(min(N, P)/2))){
      ratio <- append(ratio, values[i+1]/values[i])
    }
    ratio <- ratio[is.finite(ratio)]
    K.factors = which.min(ratio)
  } else {
    length(values)
  }
  if(K.factors>min(N, P)/2) warning('\n Warning: Number of factors supplied is > min(n,p)/2. May cause numerical inconsistencies \n')
  if(K.factors>max(N, P)) stop('\n Number of factors supplied cannot be larger than n or p \n')

  lambda.hat <- findLambdaClass(values, vectors, K.factors)
  f.hat <- findFactorClass(x, values, lambda.hat)

  if (type == "regression") {
    PF <- findPF(f.hat)
    x.res <- findXStar(xt, PF)
    y.res <- findYStar(y, PF)

    cv.model <- cv.glmnet(x.res, y.res, family=family, ...)
    coef.tmp <-  coef(cv.model,s = "lambda.min")
    betas <- rownames(coef.tmp)
    coef <- as.vector(coef.tmp)
    fitted.value <- predict(cv.model, x, s="lambda.min")
    residual <- y-fitted.value
  }

  if (type == "classification") {
    x.res <- findXStarClass(xt, lambda.hat, f.hat)
    y.res <- as.factor(y)

    cv.model <- cv.glmnet(cbind(f.hat,x.res), y.res, family=family, ...)
    F.factor <- NCOL(f.hat)
    coef.tmp <-  coef(cv.model,s = "lambda.min")
    betas <- rownames(coef.tmp)[-c(1:F.factor+1)]
    coef <- coef.tmp[-c(1:F.factor+1)]
    f.zero <- array(0,dim=dim(f.hat))
    prob <- predict(cv.model, cbind(f.zero,x), s="lambda.min", type="response")
    fitted.value <- ifelse(prob > 0.5,1,0)
    residual <- mean(y==fitted.value)
  }

  ic <- cor(y,fitted.value,use="pairwise")
  lambda <- cv.model$lambda.min
  model <- cv.model$glmnet.fit

  result <- list(coefficients=coef, betas=betas, lambda = lambda,
                 fitted.value=fitted.value, residual=residual,
                 ic=ic, model=model)

  class(result) <- "ic_farm"
  return(result)
}

#' find B matrix such that X=BxF+e, i.e. lambda.hat
#' @param eigen.value eigen values of covariance matrix
#' @param eigen.vector eigen vectors of covariance matrix
#' @param k.factors number of factors needed
findLambdaClass <- function(eigen.value, eigen.vector, k.factors)
{
  if (k.factors > 1) {
    lambda.hat <- eigen.vector[,1:k.factors]%*%diag(sqrt(eigen.value[1:k.factors]))
  } else {
    lambda.hat <- as.matrix(eigen.vector[,1]*sqrt(eigen.value[1]))
  }
  return(lambda.hat)
}

#' find F such that X = BxF + e, i.e F.hat
#' @param x predictors
#' @param eigen.value covraince matrix of x
#' @param lambda.hat the estimation of B
findFactorClass <- function(x, eigen.value, lambda.hat)
{
  k.factors <- NCOL(lambda.hat)
  if (k.factors > 1) {
    F.hat <- x %*% lambda.hat %*% diag(1/eigen.value[1:k.factors])
  } else {
    F.hat <- x * lambda.hat / eigen.value[1]
  }

  return(F.hat)
}

#' find U such that U = X - BF
#' @param x predictors
#' @param lambda.hat estimation of B
#' @param F.hat estimation of F
findXStarClass <- function(x, lambda.hat, F.hat)
{
  U.star <- x - lambda.hat %*% t(F.hat)
  U.hat <- t(U.star)
  return(U.hat)
}

#' find x star such that U = X - BF under linear regression
#' @param x predictors
#' @param PF transformation of F
findXStar <- function(x, PF)
{
  return(t(x%*%PF))
}

#' find y star
#' @param y response
#' @param PF transformation of F
findYStar <- function(y, PF)
{
  return(PF%*%y)
}

#' find PF, i.e. PF = I - F x (F^TF)^-1 x F^T
#' @param F.hat estimation of F
findPF <- function(F.hat)
{
  N <- nrow(F.hat)
  if (N > 1)
    I <- diag(N)
  else
    I <- 1

  PF <- I - F.hat %*% solve(t(F.hat)%*%F.hat) %*% t(F.hat)
  return(PF)
}

#' calculate covariance matrix based on huber loss
#' @param x predictors
#' @import sgd
#' @importFrom gtools combinations
covHuber <- function(x)
{
  P <- NCOL(x);
  index <- combinations(P, 2, repeats=T)
  x2 <- apply(index, 1, function(l) x[,l[1]]*x[,l[2]])
  data <- data.frame(y=1, x2)
  out <- sgd(y ~ .-1, data = data, model = "m",
             sgd.control=list(method="sgd", lr="adagrad", npass=10, pass=T))
  sgd.theta <- matrix(nrow = P, ncol = P)
  sgd.theta[index] <- out$coefficients
  sgd.theta[lower.tri(sgd.theta)] <- t(sgd.theta)[lower.tri(t(sgd.theta))]
  return(sgd.theta)
}


#' calculate mean based on huber loss
#' @param x predicts
#' @import sgd
meanHuber <- function(x)
{
  data <- data.frame(y=1, x)
  out <- sgd(y ~ .-1, data = data, model = "m",
             sgd.control=list(method="sgd", lr="adagrad", npass=10, pass=T))
  sgd.theta <- as.vector(out$coefficients)
  return(sgd.theta)
}
