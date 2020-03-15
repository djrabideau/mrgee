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
#' @param corstr  character string specifying the correlation structure. Allowed
#'     structures are: "independence", "exchangeable", "ar1", "m-dependent",
#'     "unstructured", "fixed", and "userdefined". Any unique substring may be
#'     supplied. If "fixed" or "userdefined", then corr.mat must be specified. If
#'     "m-dependent", then Mv is relevant.
#' @param Mv for "m-dependent", the value for m.
#' @param ... further arguments passed to \code{\link[mrgee]{modNR}}.
#' @param pmodels a list of propensity score models for \eqn{Pr(R=1|X, A)}, where
#'     \eqn{R} indicates outcome observance, X are fully observed baseline
#'     covariates, and A is the treatment indicator. List elements must be
#'     objects of class "glm" from \code{\link[stats]{glm}}.
#' @param bmodels a list of outcome models for \eqn{E(Y|X, A)}, where \eqn{Y}
#'     is the outcome, X are fully observed baseline covariates, and A is the
#'     treatment indicator. List elements must be objects of class "gee" from
#'     \code{\link[gee]{gee}}.
#' @export
mrgee <- function(formula, family = gaussian, data, id,
                  corstr = 'independence', corr.mat = NULL, init.phi = 1, scale.fix = F,
                  pmodels, bmodels, ...) {
  if (is.character(family)) {
      famret <- get(family, mode = "function", envir = parent.frame())
  } else if (is.function(family)) {
      famret <- family()
  } else if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
  }

  # set up cluster IDs
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$R <- m$b <- m$tol <- m$maxiter <- m$link <- m$varfun <- m$corstr <- m$Mv <-
    m$silent <- m$contrasts <- m$family <- m$scale.fix <- m$scale.value <-
    m$v4.4compat <- m$pmodels <- m$bmodels <- m$... <- m$corr.mat <- m$init.phi <- NULL
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
  ids <- unique(id) # unique cluster ids
  K <- length(ids) # total number of clusters
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

    # get gammahat
    gammahats <- coef(bmodel)
    bmf <- model.frame(bmodel$terms, data = data, na.action = na.pass)
    bX <- model.matrix.lm(bmf)

    # get b(gammahat), which is estimated E(Y|X,A)
    b_linkinv <- bmodel$family$linkinv
    bhat <- b_linkinv(bX %*% as.matrix(gammahats))

    # get betahat (with missing outcomes replaced by bhat)
    data$yor <- y
    data$yor[r == 0] <- bhat[r == 0]
    formulatmp <- update(formula, 'yor ~ .')
    or <- geem_Vinv(formulatmp, family = family, data = data, id = id, corstr = corstr,
                    corr.mat = corr.mat, init.phi = init.phi, scale.fix = scale.fix)
    betahats <- coef(or)

    # get mu(betahat), which is estimated E(Y|A)
    mu_linkinv <- famret$linkinv
    muhat <- mu_linkinv(or$X %*% as.matrix(betahats))

    # transformed cond/marg diffs for each cluster
    Vinv <- or$Vinv
    Uhat <- matrix(NA, nrow = length(betahats), ncol = N)
    for (k in 1:K) {
      ind <- which(id == ids[k])
      datatmp <- data[ind, ]
      difftmp <- as.vector(bhat - muhat)[ind]

      # Get Vinvk (sub)matrix
      Vinvk <- Vinv[ind, ind]

      # Get Dk matrix, now generalized beyond Dk=Xk
      InvLinkDeriv <- famret$mu.eta
      Xtmp <- model.matrix.lm(formulatmp, data = datatmp)
      etatmp <- as.vector(Xtmp %*% betahats) # maybe need to change this to work with an offset term
      Dk <- diag(InvLinkDeriv(etatmp)) %*% Xtmp

      Uhat[, ind] <- as.matrix(t(Dk) %*% Vinvk %*% diag(difftmp))
    }
    etahat <- apply(Uhat, 1, mean)
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
                       corstr = corstr, corr.mat = corr.mat, init.phi = init.phi,
                       scale.fix = scale.fix, weights = w)
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

#' Multiply Robust GEE (no V)
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
#' @param corstr  character string specifying the correlation structure. Allowed
#'     structures are: "independence", "exchangeable", "ar1", "m-dependent",
#'     "unstructured", "fixed", and "userdefined". Any unique substring may be
#'     supplied. If "fixed" or "userdefined", then corr.mat must be specified. If
#'     "m-dependent", then Mv is relevant.
#' @param Mv for "m-dependent", the value for m.
#' @param ... further arguments passed to \code{\link[mrgee]{modNR}}.
#' @param pmodels a list of propensity score models for \eqn{Pr(R=1|X, A)}, where
#'     \eqn{R} indicates outcome observance, X are fully observed baseline
#'     covariates, and A is the treatment indicator. List elements must be
#'     objects of class "glm" from \code{\link[stats]{glm}}.
#' @param bmodels a list of outcome models for \eqn{E(Y|X, A)}, where \eqn{Y}
#'     is the outcome, X are fully observed baseline covariates, and A is the
#'     treatment indicator. List elements must be objects of class "gee" from
#'     \code{\link[gee]{gee}}.
#' @export
mrgee_noV <- function(formula, family = gaussian, data, id, corstr, pmodels,
                      bmodels, ...) {
  if (is.character(family)) {
      famret <- get(family, mode = "function", envir = parent.frame())
  } else if (is.function(family)) {
      famret <- family()
  } else if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
  }

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
  ids <- unique(id) # unique cluster ids
  K <- length(ids) # total number of clusters
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

    linkinv <- bmodel$family$linkinv
    ahat <- linkinv(bX %*% as.matrix(gammahats))
    data$yor <- y
    data$yor[r == 0] <- ahat[r == 0]

    formulatmp <- update(formula, 'yor ~ .')

    or <- geem(formulatmp, family = family, data = data, id = id, corstr = 'independence') # no Vinv in OM
    # if this works, I could probably use faster way to fit this, i.e. gls()

    betahats <- coef(or)
    orhat <- as.vector(fitted.values(or))
    Uhat <- matrix(NA, nrow = length(betahats), ncol = N)
    for (k in 1:K) {
      ind <- which(id == ids[k])
      datatmp <- data[ind, ]
      difftmp <- (ahat - orhat)[ind]

      # Get Dk matrix, now generalized beyond Dk=Xk
      InvLinkDeriv <- famret$mu.eta
      Xtmp <- model.matrix.lm(formulatmp, data = datatmp)
      etatmp <- as.vector(Xtmp %*% betahats) # maybe need to change this to work with an offset term
      Dk <- diag(InvLinkDeriv(etatmp)) %*% Xtmp

      Uhat[, ind] <- as.matrix(t(Dk) %*% diag(difftmp))
    }
    etahat <- apply(Uhat, 1, mean)
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
