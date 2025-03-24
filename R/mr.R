#' Multiply Robust GEE
#'
#' Fits a multiply robust generalized estimating equation (MR-GEE) using
#' empirical maximum likelihood to handle missing outcomes.
#'
#' This function inputs a model formula, a list of propensity score models, and
#' lists of outcome regression models for the marginal means, scales, and
#' correlations to carry out MR-GEE. It outputs model fit information including
#' empirical MLEs for first- and second-order parameters (beta, gamma, alpha).
#' This function iteratively calls \code{\link[mrgee]{mrgee.fit}} to solve
#' MR-GEEs for each set of parameters (beta, gamma, alpha), corresponding to
#' parameters for the marginal mean, scales, and correlations, respectively.
#'
#' Currently this function is limited to common scalar parameters for scale and
#' correlation. It also currently assumes propensity score model specified is
#' for individual probabilities and that joint probabilities are simply product
#' of individuals (i.e. assuming independence between missingness probabilities
#' ).
#'
#' @seealso \code{\link[mrgee]{mrgee.fit}} for the utility function used to
#'     solve each MR-GEE
#'
#' @param formula an object of class "\code{\link[stats]{formula}}"
#'     (or one that can be coerced to that class): a symbolic description of the
#'     model to be fitted.
#' @param family a description of the error distribution and link function to
#'     be used in the model. This can be a character string naming a family
#'     function, a family function or the result of a call to a family function.
#'     See \code{\link[stats]{family}} for details.
#' @param data a data frame containing the variables in the model.
#' @param id a vector which identifies the clusters. The length of \code{id}
#'     should be the same as the number of observations. Data are assumed to be
#'     sorted so that observations on a cluster are contiguous rows for all
#'     entities in the formula.
#' @param corstr  character string specifying the correlation structure. Allowed
#'     structures are: "independence", "exchangeable"
#' @param scale.fix if T, then the scale parameter is fixed at init.gamma
#' @param pmodels a list of propensity score models for \eqn{Pr(R=1|X, A)}, where
#'     \eqn{R} indicates outcome observance, X are fully observed baseline
#'     covariates, and A is the treatment indicator. List elements must be
#'     objects of class "glm" from \code{\link[stats]{glm}}.
#' @param bmodels a list of outcome models for \eqn{E(Y|X, A)}, where \eqn{Y}
#'     is the outcome, X are fully observed baseline covariates, and A is the
#'     treatment indicator. List elements must be objects of class "gee" from
#'     \code{\link[gee]{gee}}.
#' @param smodels a list of outcome models for \eqn{var(Y|X, A)}
#' @param cmodels a list of outcome models for \eqn{corr(Y, Y'|X, A)}
#' @param init.gamma initial scale parameter value
#' @param init.alpha initial correlation parameter value
#' @param maxiter maximum Newton-Raphson updates of MR-GEE parameters (beta,
#'      gamma, alpha).
#' @param tol convergence tolerance for final MR-GEE updates. Maximum absolute
#'      difference between step (s-1) and (s) across all parameter updates (
#'      i.e. all betas, gamma, alpha) is compared to tol.
#' @param printIter logical whether to print iteration information
#' @export
mrgee <- function(formula, family = gaussian, data, id,
                  corstr = 'independence', scale.fix = F,
                  pmodels = NULL, bmodels = NULL, smodels = NULL, cmodels = NULL,
                  init.gamma = 1, init.alpha = 0.01, maxiter = 1, tol = 0.001,
                  printIter = F) {
  # basic checks
  if (length(pmodels) == 0 &
      (length(bmodels) == 0 | length(smodels) == 0 | length(cmodels) == 0))
    warning("no pmodel, bmodel, smodel, or cmodel specified")

  if (is.character(family)) {
      famret <- get(family, mode = "function", envir = parent.frame())
  } else if (is.function(family)) {
      famret <- family()
  } else if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
  }
  LinkFun <- famret$linkfun
  InvLink <- famret$linkinv
  VarFun <- famret$variance
  InvLinkDeriv <- famret$mu.eta

  # set up cluster IDs
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$R <- m$b <- m$tol <- m$maxiter <- m$link <- m$varfun <- m$corstr <- m$Mv <-
    m$silent <- m$contrasts <- m$family <- m$scale.fix <- m$scale.value <-
    m$v4.4compat <- m$pmodels <- m$bmodels <- m$smodels <- m$cmodels <- m$... <-
    m$corr.mat <- m$init.phi <- m$init.gamma <- m$init.alpha <- m$maxiter <-
    m$tol <- m$printIter <- NULL
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

  # initialize beta using identity matrix for working covariance
  betaFit0 <- mrgee.fit(formula = formula, family = family, data = data, id = id,
                        corstr = 'independence',
                        pmodels = pmodels, bmodels = bmodels,
                        scale = 1, corr = 0)
  beta <- as.matrix(betaFit0$coefficients)

  # inits for gamma, alpha
  gamma <- init.gamma
  alpha <- init.alpha

  betas <- matrix(NA, nrow = nrow(beta), ncol = 0)
  gammas <- c()
  alphas <- c()
  for (it in 1:maxiter) {
    # previous values
    beta.last <- beta
    gamma.last <- gamma
    alpha.last <- alpha

    # scales
    eta <- as.matrix(X %*% beta)
    etavec <- as.vector(eta)
    mu <- InvLink(eta)
    v <- VarFun(mu)
    e <- (y - mu) / sqrt(v)
    phi <- rep(gamma, N)
    s <- data$s <- e^2
    if (!scale.fix) {
      gammaFit0 <- mrgee.fit(formula = s ~ 1, family = gaussian, id = id, data = data,
                             corstr = 'independence', pmodels,
                             scale = 1, corr = 0)
      gamma <- as.numeric(gammaFit0$coefficients)
    }
    gammas <- c(gammas, gamma)

    # correlations
    if (corstr == 'exchangeable') {
      alpha <- init.alpha
      nks <- as.vector(table(id))
      Nstar <- sum(choose(nks, 2))
      id2 <- rep(ids, choose(nks, 2))
      rho <- rep(alpha, Nstar)
      z <- c()
      for (i in seq_along(ids)) {
        idtmp <- ids[i]
        rows <- which(id == idtmp)
        etmp <- e[rows]
        ecomb <- combn(etmp, 2)
        eprods <- apply(ecomb, 2, prod)
        z <- c(z, eprods)
      }
      z <- z / gamma
      data2 <- data.frame(z, id2)
      alphaFit0 <- mrgee.fit(formula = z ~ 1, family = gaussian, id = id2, data = data2,
                             corstr = 'independence', pmodels = pmodels,
                             scale = 1, corr = 0, jointP = T, indiv_id = id)
      alpha <- as.numeric(alphaFit0$coefficients)
    } else {
      alpha <- 0
    }
    alphas <- c(alphas, alpha)

    # update beta given (gamma, alpha)
    betaFit <- mrgee.fit(formula = formula, family = family, data = data, id = id,
                         corstr = corstr, pmodels = pmodels, bmodels = bmodels,
                         scale = gamma, corr = alpha)
    beta <- betaFit$coefficients
    betas <- cbind(betas, beta)

    # check convergence
    converged <- NA
    if (maxiter > 1) {
      beta_diff <- max(abs(beta - beta.last))
      gamma_diff <- abs(gamma - gamma.last)
      alpha_diff <- abs(alpha - alpha.last)
      diffs <- c(beta_diff, gamma_diff, alpha_diff)
      if (printIter)
        print(round(c(beta, gamma, alpha, beta_diff, gamma_diff, alpha_diff), 5))
      if (max(diffs) < tol) {
        converged <- T
        break
      } else if (it == maxiter) {
        if (printIter) {
          print(paste0('beta_diff = ', round(beta_diff, 5)))
          print(paste0('gamma_diff = ', round(gamma_diff, 5)))
          print(paste0('alpha_diff = ', round(alpha_diff, 5)))
        }
        converged <- F
        warning(paste0('did not converge after ', it, ' iterations'))
      }
    }
  }

  # output
  betaOut <- beta
    names(betaOut) <- dimnames(X)[[2]]

  out <- list(
    call = call,
    ghat = betaFit$ghat,
    rhohat = betaFit$rhohat,
    rhohat_converged = betaFit$rhohat_converged,
    weights = betaFit$w,
    sumw = betaFit$sumw,
    fit = betaFit,
    coefficients = betaOut,
    gamma = gamma,
    alpha = alpha,
    converged = converged # NA means only 1 iter, F means didn't, T means did
  )
  class(out) <- 'mrgee'
  return(out)
}

#' Fit Single Level of Multiply Robust GEE
#'
#' A utility function used by \code{\link[mrgee]{mrgee}} to solve a single level
#' of MR-GEE.
#'
#' This function solves a single weighted GEE corresponding to either the set of
#' marginal mean, scale, or correlation parameters in MR-GEE. This function is
#' called iteratively by \code{\link[mrgee]{mrgee}}. If used for correlation
#' parameter GEE, \code{jointP} should be TRUE and \code{indiv_id} should be
#' additionally specified.
#'
#' @seealso \code{\link[mrgee]{mrgee}} for the general fitting function for
#' MR-GEE
#'
#' @inheritParams mrgee
#' @param bmodels a list of outcome models for \eqn{E(Y|X, A)}, where \eqn{Y}
#'     is the outcome, X are fully observed baseline covariates, and A is the
#'     treatment indicator. List elements must be objects of class "gee" from
#'     \code{\link[gee]{gee}}.
#' @param smodels a list of outcome models for \eqn{var(Y|X, A)}
#' @param cmodels a list of outcome models for \eqn{corr(Y, Y'|X, A)}
#' @param scale current value of scale parameter for mrgee.fit()
#' @param corr current value of correlation parameter for mrgee.fit()
#' @param jointP if jointP == T, we're requesting joint probabilities,
#'      currently setup as just the products. Should be T if fitting MR-GEE for
#'      correlation parameters.
#' @param indiv_id if provided, this corresponds to individual (rather than
#'      joint) vector of ids. Should be provided if fitting MR-GEE for
#'      correlation parameters.
#' @param ... further arguments passed to \code{\link[mrgee]{modNR}}.
#' @export
mrgee.fit <- function(formula, family, data, id, corstr = 'independence',
                      pmodels = NULL, bmodels = NULL, scale, corr, jointP = F,
                      indiv_id = NULL, ...) {

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
    m$v4.4compat <- m$pmodels <- m$bmodels <- m$... <- m$corr.mat <- m$init.phi
  m$scale <- m$corr <- m$jointP <- m$indiv_id <- NULL
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

  # PS models - for now, just individual probs, not joint
  J <- length(pmodels)
  if (J > 0) {
    pihats_c <- matrix(NA, nrow = J, ncol = N)
    if (!jointP) { # individual  probabilities
      for (j in seq_along(pmodels)) {
        pmodel <- pmodels[[j]]
        pihats <- predict(pmodel, type = 'response')
        thetahat <- mean(pihats)
        pihats_c[j, ] <- pihats - thetahat
      }
    } else { # joint probabilities = product, assuming independence R_ki R_ki'
      nks <- as.vector(table(indiv_id))
      for (j in seq_along(pmodels)) {
        pmodel <- pmodels[[j]]
        pihats <- predict(pmodel, type = 'response')
        pipihats <- c()
        for (i in seq_along(ids)) {
          idtmp <- ids[i]
          rows <- which(indiv_id == idtmp)
          pitmp <- pihats[rows]
          picomb <- combn(pitmp, 2)
          piprods <- apply(picomb, 2, prod)
          pipihats <- c(pipihats, piprods)
        }
        thetahat <- mean(pipihats)
        pihats_c[j, ] <- pipihats - thetahat
      }
    }
  } else {
    pihats_c <- NULL
  }

  # OR models
  L <- length(bmodels)
  if (L > 0) {
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
      yor <- data$yor <- y
      yor[r == 0] <- data$yor[r == 0] <- bhat[r == 0]
      formulatmp <- update(formula, 'yor ~ .')
      or <- gee.fit(X, yor, family = family, id = id, corstr = corstr,
                    scale = scale, corr = corr, calcVinv = T)
      betahats <- or$coefficients

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
        Xtmp <- model.matrix.lm(update(formula, 'yor ~ .'), data = datatmp)
        etatmp <- as.vector(Xtmp %*% betahats) # maybe need to change this to work with an offset term
        Dk <- diag(InvLinkDeriv(etatmp)) %*% Xtmp

        Uhat[, ind] <- as.matrix(t(Dk) %*% Vinvk %*% diag(difftmp))
      }
      etahat <- apply(Uhat, 1, mean)
      Uhats_c[rowStart:rowEnd, ] <- Uhat - etahat
    }
  } else {
    Uhats_c <- NULL
  }

  ghat <- rbind(pihats_c, Uhats_c)

  if (!is.null(ghat)) {
    # modified Newton-Raphson
    modNR_results <- modNR(ghat, r, ...)
    rhohat <- modNR_results$rhohat
    rhohat_converged <- modNR_results$converged

    # weighted GEE
    w <- data$w <- as.vector((1 / M) * (1 / (1 + t(rhohat) %*% ghat))) * r
  } else { # if no models provided, unconstrained maximization results in this:
    rhohat <- rhohat_converged <- NULL
    w <- data$w <- (1 / M) * r
  }
  sumw <- sum(w) # should sum to 1

  fit_mr <- gee.fit(X, y, family = family, id = id, corstr = corstr,
                    weights = w, scale = scale, corr = corr)
  coefmr <- fit_mr$coefficients

  out <- list(
    call = call,
    ghat = ghat,
    rhohat = rhohat,
    rhohat_converged = rhohat_converged,
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
