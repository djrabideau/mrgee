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
    tau <- 1 / sqrt(s)
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
