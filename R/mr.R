#' Multiply Robust GEE
#'
#' Fits a multiply robust generalized estimating equation (MR-GEE) using
#' empirical maximum likelihood to handle missing outcomes.
#'
#' This function inputs a model formula as well as lists of propensity score
#' models and outcome regression models to carry out MR-GEE. It outputs model
#' fit information including empirical MLEs for the coefficients. The weighted
#' GEE is fit using \code{\link[geeM]{geem}}, which ensures consistency of this
#' weighted estimator, rather than the conventional  implementation in PROC
#' GENMOD for example.
#'
#' @param formula an object of class "\code{\link[stats]{formula}}"
#'     (or one that can be coerced to that class): a symbolic description of the
#'     model to be fitted. This argument is passed to
#'     \code{\link[geeM]{geem}}
#' @param family a description of the error distribution and link function to
#'     be used in the model. This can be a character string naming a family
#'     function, a family function or the result of a call to a family function.
#'     This argument is passed to \code{\link[geeM]{geem}}. See
#'     \code{\link[stats]{family}} for details.
#' @param data a data frame containing the variables in the model. This argument
#'     is passed to \code{\link[geeM]{geem}}
#' @param id a vector which identifies the clusters. The length of \code{id}
#'     should be the same as the number of observations. Data are assumed to be
#'     sorted so that observations on a cluster are contiguous rows for all
#'     entities in the formula.
#' @param corstr a character string specifying the correlation structure. The
#'     following are permitted: "independence", "exchangeable", "ar1",
#'     "unstructured", and "userdefined".
#' @param ... further arguments passed to or from other methods
#' @param pmodels a list of propensity score models for \eqn{Pr(R=1|X, A)}, where
#'     \eqn{R} indicates outcome observance, X are fully observed baseline
#'     covariates, and A is the treatment indicator. List elements must be
#'     objects of class "glm" from \code{\link[stats]{glm}}.
#' @param bmodels a list of outcome models for \eqn{E(Y|X, A)}, where \eqn{Y}
#'     is the outcome, X are fully observed baseline covariates, and A is the
#'     treatment indicator. List elements must be objects of class "gee" from
#'     \code{\link[gee]{gee}}.
#' @export
mrgee <- function(formula, family = gaussian, data, id, corstr, pmodels,
                  bmodels, ...) {
  # errors/warnings
  if (corstr != 'exchangeable')
    stop('Currently only works with corstr = "exchangeable"')

  # set up cluster IDs
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$R <- m$b <- m$tol <- m$maxiter <- m$link <- m$varfun <- m$corstr <- m$Mv <-
    m$silent <- m$contrasts <- m$family <- m$scale.fix <- m$scale.value <-
    m$v4.4compat <- m$pmodels <- m$bmodels <- m$... <- NULL
  m$na.action <- as.name("na.pass")
  if (is.null(m$id))
    m$id <- as.name("id")
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  id <- model.extract(m, id)
  if ('id' %in% names(data))
    data <- data[, -which(names(data) == 'id')]
  data$id <- id

  N <- nrow(data) # total number of observations
  K <- length(unique(id)) # total number of clusters
  mf <- model.frame(formula, data, na.action = na.pass)
  y <- mf[, 1] # outcome vector
  r <- as.numeric(!is.na(y)) # observance indicator
  M <- sum(r) # total number of observed outcomes
  X <- model.matrix.lm(mf) # design matrix
  p <- ncol(X)

  # PS models
  J <- length(pmodels)
  pihats_c <- matrix(NA, nrow = J, ncol = N)
  for (j in seq_along(pmodels)) {
    pmodel <- pmodels[[j]]
    pihats <- predict(pmodel, type = 'response')
    thetahat <- mean(pihats)
    pihats_c[j, ] <- pihats - thetahat
  }

  # OR models
  L <- length(bmodels)
  Uhats_c <- matrix(NA, nrow = p * L, ncol = N)
  for (l in seq_along(bmodels)) {
    rowStart <- (l - 1) * p + 1
    rowEnd <- l * p
    bmodel <- bmodels[[l]]

    gammahats <- coef(bmodel)
    bmf <- model.frame(bmodel$terms, data = data, na.action = na.pass)
    bX <- model.matrix.lm(bmf)

    ahat <- bX %*% as.matrix(gammahats)
    data$yor <- y
    data$yor[r == 0] <- ahat[r == 0]

    formulatmp <- update(formula, 'yor ~ .')

    sink(tempfile()) # suppress gee output
      or <- suppressMessages( # suppress gee messages
        gee(formulatmp, family = family, data = data, id = id, corstr = corstr)
        )
    sink()

    betahats <- coef(or)
    orhat <- as.vector(fitted.values(or))
    Vk <- or$scale * or$working.correlation # need to generalize this
    zetahat <- matrix(NA, nrow = length(betahats), ncol = K)
    Uhat <- matrix(NA, nrow = length(betahats), ncol = N)
    for (k in 1:K) {
      ind <- which(id == k)
      datatmp <- data[ind, ]
      difftmp <- (ahat - orhat)[ind]
      Dk <- model.matrix.lm(formulatmp, data = datatmp)
      zetahat[, k] <- t(Dk) %*% solve(Vk) %*% difftmp
      Uhat[, ind] <- t(Dk) %*% solve(Vk) %*% diag(difftmp)
    }
    # etahat <- (1 / N) * apply(zetahat, 1, sum) # if N^{-1} for eta terms
    etahat <- apply(zetahat, 1, mean) # if K^{-1} for eta terms
    Uhats_c[rowStart:rowEnd, ] <- Uhat - etahat
  }

  ghat <- rbind(pihats_c, Uhats_c)

  # modified Newton-Raphson
  modNR_results <- modNR(ghat, r, ...)
  rhohat <- modNR_results$rhohat

  # weighted GEE
  w <- data$w <- as.vector((1 / M) * (1 / (1 + t(rhohat) %*% ghat))) * r
  sumw <- sum(w) # should sum to 1

  fit_mr <- geeM::geem(formula = formula, family = family, data = data, id = id,
                       corstr = corstr, weights = w)
  coefmr <- coef(fit_mr)

  out <- list(
    call = call,
    ghat = ghat,
    rhohat = rhohat,
    rhohat_converged = modNR_results$converged,
    weights = w,
    sumw = sumw,
    fit = fit_mr,
    coefficients = coefmr
  )
  class(out) <- 'mrgee'
  return(out)
}

#' Multiply Robust GLM
#'
#' Fits a multiply robust generalized linear model (MR-GLM) for independent data
#' using empirical maximum likelihood to handle missing outcomes. See Han
#' (2014).
#'
#' This function inputs a model formula as well as lists of propensity score
#' models and outcome regression models to carry out multiply robust regression.
#' It outputs model fit information including empirical MLEs for the
#' coefficients.
#'
#' @param formula an object of class "\code{\link[stats]{formula}}"
#'     (or one that can be coerced to that class): a symbolic description of the
#'     model to be fitted. This argument is passed to
#'     \code{\link[stats]{glm}}
#' @param family a description of the error distribution and link function to
#'     be used in the model. This can be a character string naming a family
#'     function, a family function or the result of a call to a family function.
#'     This argument is passed to \code{\link[stats]{glm}}. See
#'     \code{\link[stats]{family}} for details.
#' @param data a data frame containing the variables in the model. This argument
#'     is passed to \code{\link[stats]{glm}}
#' @param pmodels a list of propensity score models for \eqn{Pr(R=1|X, S)}, where
#'     \eqn{R} indicates outcome observance, X are covariates, and S are
#'     auxiliary variables. List elements must be objects of class "glm" from
#'     \code{\link[stats]{glm}}.
#' @param bmodels a list of outcome models for \eqn{E(Y|X, A)}, where \eqn{Y}
#'     is the outcome, X are fully observed baseline covariates, and A is the
#'     treatment indicator. List elements must be objects of class "glm" from
#'     \code{\link[stats]{glm}}.
#' @param ... further arguments passed to or from other methods
#' @export
mrglm <- function(formula, family = gaussian, data, pmodels, bmodels, ...) {
  # set up cluster IDs
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$pmodels <- m$bmodels <- m$... <- NULL
  m$na.action <- as.name("na.pass")
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())

  N <- nrow(data) # total number of observations
  mf <- model.frame(formula, data, na.action = na.pass)
  y <- mf[, 1] # outcome vector
  r <- as.numeric(!is.na(y)) # observance indicator
  M <- sum(r) # total number of observed outcomes
  X <- model.matrix.lm(mf) # design matrix
  p <- ncol(X)

  # PS models
  J <- length(pmodels)
  pihats_c <- matrix(NA, nrow = J, ncol = N)
  for (j in seq_along(pmodels)) {
    pmodel <- pmodels[[j]]
    pihats <- predict(pmodel, type = 'response')
    thetahat <- mean(pihats)
    pihats_c[j, ] <- pihats - thetahat
  }

  # OR models
  L <- length(bmodels)
  Uhats_c <- matrix(NA, nrow = p * L, ncol = N)
  for (l in seq_along(bmodels)) {
    rowStart <- (l - 1) * p + 1
    rowEnd <- l * p
    bmodel <- bmodels[[l]]

    gammahats <- coef(bmodel)
    bmf <- model.frame(bmodel$terms, data = data, na.action = na.pass)
    bX <- model.matrix.lm(bmf)

    ahat <- as.vector(bX %*% as.matrix(gammahats))
    data$yor <- y
    data$yor[r == 0] <- ahat[r == 0]

    formulatmp <- update(formula, 'yor ~ .')

    or <- glm(formulatmp, family = family, data = data)

    betahats <- coef(or)
    orhat <- predict(or)
    Uhat <- diag(ahat - orhat) %*% X
    etahat <- apply(Uhat, 2, mean)
    Uhats_c[rowStart:rowEnd, ] <- t(Uhat - etahat)
  }

  ghat <- rbind(pihats_c, Uhats_c)

  # modified Newton-Raphson
  modNR_results <- modNR(ghat, r, ...)
  rhohat <- modNR_results$rhohat

  # weighted linear model
  w <- data$w <- as.vector((1 / M) * (1 / (1 + t(rhohat) %*% ghat))) * r
  sumw <- sum(w) # should sum to 1

  fit_mr <- glm(formula, family = family, data = data, weights = w)
  coefmr <- coef(fit_mr)

  out <- list(
    call = call,
    ghat = ghat,
    rhohat = rhohat,
    rhohat_converged = modNR_results$converged,
    weights = w,
    sumw = sumw,
    fit = fit_mr,
    coefficients = coefmr
  )
  class(out) <- 'mrglm'
  return(out)
}
