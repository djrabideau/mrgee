#' Modified Newton-Raphson
#'
#' Carries out modified Newton-Raphson algorithm to solve for \eqn{\rho} for
#' MR-GEE.
#'
#' This function inputs a \eqn{(J + 2L) x N} matrix \code{ghat} and outcome
#' observance \eqn{N}-vector \code{r}, performs modified Newton-Raphson, and
#' returns a \eqn{(J + 2L)}-vector \eqn{\hat{\rho}]}.
#'
#' @param ghat a \eqn{(J + 2L) x N} matrix, where \eqn{J} is the number of
#' propensity score models, \eqn{L} is the number of outcome regression models,
#' \eqn{N} is the total number of individuals in the data set (both observed and
#' unobserved). Row order must correspond to \code{r}
#' @param r an \eqn{N}-vector of 0s and 1s indicating which individuals had
#' their outcome observed (1) or missing (0). Order must correspond to rows in
#' \code{ghat}
#' @param maxit maximum number of iterations for algorithm; default
#' \code{maxit=100}
#' @param eps convergence criterion for algorithm; default \code{eps=10e-8}
#' @param quietly logical; if F, convergence info for each step is printed
#' @export
modNR <- function(ghat, r, maxit = 100, eps = 10e-8, quietly = T) {
  call <- match.call()
  N <- length(r)
  M <- sum(r)
  observed <- which(as.logical(r))

  # function to minimize
  f <- function(x) {
    -(1 / M) * sum(as.vector(log(1 + t(x) %*% ghat[, observed])))
  }

  # step 0
  rhohat <- matrix(rep(0, nrow(ghat)), ncol = 1)
  s <- 0
  tau <- 1

  converged <- F
  for (j in 1:maxit) {
    # step 1
    delta1 <- rep(0, nrow(ghat))
    for (i in observed) {
      tmp <- ghat[, i] / as.numeric(1 + t(rhohat) %*% ghat[, i])
      delta1 <- delta1 + tmp
    }
    delta1 <- as.matrix(delta1)

    sumouter <- matrix(0, nrow = nrow(ghat), ncol = nrow(ghat))
    for (i in observed) {
      tmp <- (ghat[, i] %*% t(ghat[, i])) /
        as.numeric(1 + t(rhohat) %*% ghat[, i])^2
      sumouter <- sumouter + tmp
    }
    delta2 <- solve(-sumouter) %*% delta1

    if (!quietly)
      print(paste0('step ', j, ', norm = ', base::norm(delta2, type = 'F')))
    if (base::norm(delta2, type = 'F') < eps) {
      converged <- T
      break
    }

    # step 2
    condA <- condB <- T
    maxit2 <- 0 # in case of a runaway while loop
    while ((condA | condB) & maxit2 <= 100) {
      delta <- tau * delta2
      rhotmp <- rhohat - delta
      condA <- sum((1 + t(rhotmp) %*% ghat[, observed]) <= 0) > 0
      condB <- F
      if (!condA) condB <- f(rhotmp) > f(rhohat)
      tau <- tau / 2
      maxit2 <- maxit2 + 1
      if (maxit2 == 100)
        stop("step 2 conditions not met after 100 restarts")
    }

    # step 3
    rhohat <- rhohat - delta
    s <- s + 1
    # tau <- 1 / sqrt(s) # need this for convergence proof
    tau <- 1 # "numerical behavior not affected with this instead" - Han (2014)
            # also seems to converge faster with tau <- 1
  }

  if (!converged)
    warning(paste0('Search for rhohat did not converge with eps = ', eps,
                   ', maxit = ', maxit))

  out <- list(
    call = call,
    rhohat = rhohat,
    converged = converged,
    args = list(
      maxit = maxit,
      eps = eps
    )
  )
  return(out)
}

#' Update parameters of GEE using weighted Newton-Raphson update
#'
#' Carries out one step of the usual GEE Newton-Raphson iterative fitting
#' procedure
#'
#' @param beta vector of parameters to update, e.g. \eqn{\beta} for marginal
#' means, \eqn{\gamma} for scales, \eqn{\alpha} for correlations
#' @param D matrix of partial derivatives
#' @param V working covariance matrix
#' @param W weight matrix
#' @param y vector of outcomes or residuals, e.g. \eqn{Y[k]}, \eqn{S[k]},
#' \eqn{Z[k]}
#' @param mu vector of marginal means, scales, or correlations, e.g.
#' \eqn{\mu[k]}, \eqn{\phi[k]}, \eqn{\rho[k]}
#' @param id vector of cluster identifiers
#' @export
gee_update <- function(beta, D, V, W, y, mu, id, noInv) {
  p <- length(beta)
  ids <- unique(id)
  sumHess <- matrix(0, nrow = p, ncol = p)
  sumGrad <- matrix(0, nrow = p, ncol = 1)
  for (i in seq_along(ids)) {
    idtmp <- ids[i]
    rows <- which(id == idtmp)
    Dk <- D[rows, ]
    Vk <- V[rows, rows]
    Wk <- W[rows, rows]
    yk <- y[rows]
    muk <- mu[rows]
    if (!noInv) {
      Vkinv <- solve(Vk)
    } else {
    Vkinv <- Vk
    }

    hessk <- t(Dk) %*% Vkinv %*% Wk %*% Dk
    gradk <- t(Dk) %*% Vkinv %*% Wk %*% (yk - muk)

    sumHess <- sumHess + hessk
    sumGrad <- sumGrad + gradk
  }
  betaNext <- beta + solve(sumHess) %*% sumGrad
  return(betaNext)
}

#' Update first-order working covariance
#'
#' Updates the entire working covariance matrix \eqn{V[1]} with most recent
#' second order parameter values, i.e. \eqn{(\gamma, \alpha)}
#'
#' @param mu vector of marginal means
#' @param gamma scale parameters
#' @param alpha correlation parameters
#' @param id vector of cluster identifiers
#' @param corstr a character string specifying the correlation structure.
#' Allowed structures are: "independence" and "exchangeable"
#' @param VarFun function: the variance as a function of the mean
#' @export
updateV <- function(mu, gamma, alpha, id, corstr, VarFun) {
  muvec <- as.vector(mu)
  ids <- unique(id)
  N <- length(id)
  blocks <- list()
  for (i in seq_along(ids)) {
    idtmp <- ids[i]
    rows <- which(id == idtmp)
    nk <- length(rows)
    muvectmp <- muvec[rows]
    Akhalf <- Matrix::Diagonal(nk, sqrt(VarFun(muvectmp)))
    Ck <- Matrix::Diagonal(nk)
    if (corstr == 'exchangeable') {
      Ck[upper.tri(Ck)] <- alpha
      Ck[lower.tri(Ck)] <- alpha
    }
    blocks[[i]] <- gamma * (Akhalf %*% Ck %*% Akhalf)
  }
  V <- Matrix::bdiag(blocks)
  return(V)
}

#' Fit weighted GEE (point estimates only)
#'
#' This utility function fits a weighted GEE and returns parameter point
#' estimates, design matrix, and inverse of the working covariance matrix. We
#' solve for the regression parameters only; scale and correlation parameters
#' are assumed fixed at their given values. Note, this function can be used to
#' solve for second-order regression parameters such as scales or correlations;
#' in that case, scale and corr arguments of this function correspond to even
#' higher-order parameters (e.g. third and fourth moment quantities).
#'
#' @param X design matrix
#' @param y outcome vector
#' @param family a description of the error distribution and link function to
#'     be used in the model. This can be a character string naming a family
#'     function, a family function or the result of a call to a family function.
#'     See \code{\link[stats]{family}} for details.
#' @param id a vector which identifies the clusters. The length of \code{id}
#'     should be the same as the number of observations. Data are assumed to be
#'     sorted so that observations on a cluster are contiguous rows for all
#'     entities in the formula.
#' @param corstr  character string specifying the correlation structure. Allowed
#'     structures are: "independence", "exchangeable"
#' @param maxiter maximum Newton-Raphson updates of GEE parameters.
#' @param tol convergence tolerance for final GEE updates. Maximum absolute
#'      difference between step (s-1) and (s) across all parameter values.
#' @param weights vector of weights corresponding to \eqn{y}
#' @param scale current value of scale parameter
#' @param corr current value of correlation parameter
#' @param calcVinv if T, calculates and outputs inverse of entire covariance
#'      matrix V (note, this calculation could be very slow for correlation
#'      parameters if cluster sizes are large, say >30 observations/cluster).
#' @export
gee.fit <- function(X, y, family, id, corstr,
                    maxiter = 25, tol = 0.001, weights = NULL,
                    scale, corr, calcVinv = F) {
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

  # check supported corstr
  if (!(corstr %in% c('independence', 'exchangeable'))) {
    stop(paste0('corstr = ', corstr, ' not supported'))
  }

  # weights
  if (is.null(weights))
    weights <- ifelse(is.na(y), 0, 1)
  w <- weights
  sumw <- sum(w)

  # drop 0 weights
  dropRows <- which(weights == 0)
  if (length(dropRows) > 0) {
    X <- X[-dropRows, ]
    y <- y[-dropRows]
    id <- id[-dropRows]
    w <- w[-dropRows]
  }

  # setup data
  N <- length(y) # total number of observations
  ids <- unique(id) # unique cluster ids
  K <- length(ids) # total number of clusters
  p <- ncol(X)
  W <- Matrix::Diagonal(N, w)

  # init betahat
  fit0 <- suppressWarnings(glm.fit(X, y, family = famret, weights = w))
  beta <- as.matrix(coef(fit0))
  # ...suppressWarnings when binary and non-integer predicted values for miss outcomes

  # iteratively solve GEE
  for (i in 1:maxiter) {
    # calculate mu, eta, v, e for most recent updates
    eta <- as.matrix(X %*% beta)
    etavec <- as.vector(eta)
    mu <- InvLink(eta)
    v <- VarFun(mu)

    # update matrices
    D <- as.matrix(Matrix::Diagonal(N, InvLinkDeriv(etavec)) %*% X)
    V <- updateV(mu, gamma = scale, alpha = corr, id, corstr, VarFun)
    # if V is identity matrix, no need to slow things down by taking inverse
    # when doing gee_update below (i.e. when fitting scale and corr params)
    if (corstr == 'independence' & scale == 1) {
      noInv <- T
    } else {
      noInv <- F
    }

    # update beta
    beta.last <- beta
    beta <- gee_update(beta, D, V, W, y, mu, id, noInv)

    # check convergence
    if (max(abs(beta - beta.last)) < tol) {
      break
    } else if (i == maxiter) {
      warning(paste0('init.betahat did not converge after ', i, ' iterations'))
    }
  }

  betahat <- as.vector(beta)
  if (calcVinv) {
    Vinv <- solve(V)
  } else {
    Vinv <- NULL
  }
  out <- list(
    coefficients = betahat,
    X = X,
    Vinv = Vinv
  )
  return(out)
}
