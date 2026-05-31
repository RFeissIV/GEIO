# ============================================================================
# GEIO -- Equilibrium Solver, Index Computation, and Fragility Diagnostics
# ============================================================================
# Author : Richard A. Feiss IV, Ph.D.
# Contact: feiss026@umn.edu
# Org    : Minnesota Center for Prion Research and Outreach (MNPRO)
#          University of Minnesota
# License: MIT
# ============================================================================
#
# This file implements the equilibrium-level machinery:
#
#   1. Newton's method for F_theta(x) = 0
#   2. Equilibrium index:  chi = sgn det D_x F_theta(x*)
#   3. Fragility:          sigma_min of D_x F_theta(x*)
#   4. Numerical Jacobian:  central-difference approximation
#
# The index theory follows Shapley (1974) and Govindan-Wilson (2003).
# For a regular real zero, the local degree reduces to sgn det J
# (the Eisenbud-Khimshiashvili-Levine signature), which provides the
# local-degree background for the determinant-sign orientation
# diagnostic.  GEIO uses this only as a diagnostic at recovered regular
# zeros, not as an exhaustive theorem for arbitrary custom game systems.
# ============================================================================


# ============================================================================
# Numerical Jacobian (central difference)
# ============================================================================

#' Numerical Jacobian via central differences
#'
#' Computes the Jacobian matrix of a vector-valued function \code{F} at
#' point \code{x} using central finite differences.  Each column
#' \eqn{j} is approximated as
#' \eqn{(F(x + e_j h) - F(x - e_j h)) / (2h)} where
#' \eqn{h = \text{eps} \cdot (|x_j| + 1)}.
#'
#' @param F Vector-valued function: \code{function(x) -> numeric vector}
#'   of length \eqn{m}.
#' @param x Numeric vector of length \eqn{n} at which to evaluate.
#' @param eps Step size scaling factor (default \code{1e-7}).
#' @return Numeric matrix of dimension \eqn{m \times n}.
#'
#' @examples
#' F <- function(x) c(x[1]^2 + x[2] - 1, x[1] + x[2]^2 - 1)
#' J <- geio_numjac(F, c(0.5, 0.5))
#' # Analytical: [[2*x1, 1], [1, 2*x2]] = [[1, 1], [1, 1]]
#'
#' @export
geio_numjac <- function(F, x, eps = 1e-7) {
  if (!is.function(F))
    stop("geio_numjac: F must be a function")
  x  <- as.numeric(x)
  if (!all(is.finite(x)))
    stop("geio_numjac: x must be a finite numeric vector")
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0)
    stop("geio_numjac: eps must be a positive finite scalar")
  n  <- length(x)
  if (n < 1L)
    stop("geio_numjac: x must have positive length")
  f0 <- F(x)
  if (!is.numeric(f0) || length(f0) < 1L || !all(is.finite(f0)))
    stop("geio_numjac: F(x) must return a finite numeric vector")
  m  <- length(f0)
  J  <- matrix(NA_real_, nrow = m, ncol = n)

  for (j in seq_len(n)) {
    h      <- eps * (abs(x[j]) + 1)
    x_p    <- x;  x_p[j] <- x[j] + h
    x_m    <- x;  x_m[j] <- x[j] - h
    fp <- F(x_p)
    fm <- F(x_m)
    if (!is.numeric(fp) || !is.numeric(fm) ||
        length(fp) != m || length(fm) != m ||
        !all(is.finite(fp)) || !all(is.finite(fm)))
      stop("geio_numjac: F must return finite numeric vectors of consistent length")
    J[, j] <- (fp - fm) / (2 * h)
  }

  if (!all(is.finite(J)))
    stop("geio_numjac: non-finite finite-difference derivative encountered")
  J
}


# ============================================================================
# Newton solver for F_theta(x) = 0
# ============================================================================

#' Newton's method for equilibrium systems
#'
#' Finds \eqn{x^*} such that \eqn{F(x^*) \approx 0} using damped
#' Newton iteration with backtracking.  If the Jacobian is singular or
#' near-singular, a Tikhonov regularization term is added.
#'
#' @param F Vector-valued function: \code{function(x) -> numeric vector}.
#' @param x0 Initial guess, numeric vector.
#' @param DF Optional Jacobian function: \code{function(x) -> matrix}.
#'   If \code{NULL}, \code{\link{geio_numjac}} is used.
#' @param max_iter Maximum Newton iterations (default 200).
#' @param tol Convergence tolerance on \eqn{\|F(x)\|_\infty}
#'   (default \code{1e-10}).
#' @param damping Initial damping factor for backtracking (default 1.0).
#' @param regularize Tikhonov regularization added to diagonal when
#'   the Jacobian condition number exceeds \code{cond_max}
#'   (default \code{1e-8}).
#' @param cond_max Condition number threshold triggering regularization
#'   (default \code{1e12}).
#'
#' @return A named list:
#' \describe{
#'   \item{\code{x}}{Numeric vector; the equilibrium point (or last iterate).}
#'   \item{\code{converged}}{Logical; \code{TRUE} if
#'     \eqn{\|F(x)\|_\infty < \text{tol}}.}
#'   \item{\code{residual}}{Scalar; \eqn{\|F(x)\|_\infty} at termination.}
#'   \item{\code{iterations}}{Integer; Newton iterations used.}
#'   \item{\code{jacobian}}{Matrix; Jacobian \eqn{D_x F} at the final
#'     iterate.}
#' }
#'
#' @examples
#' # Intersection of two parabolas: x^2 + y = 1, x + y^2 = 1
#' F <- function(x) c(x[1]^2 + x[2] - 1, x[1] + x[2]^2 - 1)
#' sol <- geio_newton(F, c(0.5, 0.5))
#' sol$x          # should be near (1, 0) or (0, 1)
#' sol$converged
#'
#' @export
geio_newton <- function(F, x0, DF = NULL, max_iter = 200L, tol = 1e-10,
                         damping = 1.0, regularize = 1e-8,
                         cond_max = 1e12) {

  if (!is.function(F))
    stop("geio_newton: F must be a function")
  x <- as.numeric(x0)
  if (!all(is.finite(x)))
    stop("geio_newton: x0 must be a finite numeric vector")
  n <- length(x)
  if (n < 1L)
    stop("geio_newton: x0 must have positive length")
  if (!is.numeric(max_iter) || length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1)
    stop("geio_newton: max_iter must be a positive integer-like scalar")
  max_iter <- as.integer(max_iter)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0)
    stop("geio_newton: tol must be a positive finite scalar")
  if (!is.null(DF) && !is.function(DF))
    stop("geio_newton: DF must be NULL or a function")
  if (is.null(DF)) DF <- function(xx) geio_numjac(F, xx)
  if (!is.numeric(damping) || length(damping) != 1L || !is.finite(damping) || damping <= 0)
    stop("geio_newton: damping must be a positive finite scalar")
  if (!is.numeric(regularize) || length(regularize) != 1L || !is.finite(regularize) || regularize < 0)
    stop("geio_newton: regularize must be a non-negative finite scalar")
  if (!is.numeric(cond_max) || length(cond_max) != 1L || !is.finite(cond_max) || cond_max <= 0)
    stop("geio_newton: cond_max must be a positive finite scalar")

  .safe_DF <- function(xx) {
    out <- tryCatch(DF(xx), error = function(e) matrix(NA_real_, n, n))
    as.matrix(out)
  }

  for (k in seq_len(max_iter)) {
    Fx <- F(x)
    if (!all(is.finite(Fx)))
      return(list(x = x, converged = FALSE, residual = Inf,
                  iterations = k, jacobian = .safe_DF(x)))

    res <- max(abs(Fx))
    if (res < tol)
      return(list(x = x, converged = TRUE, residual = res,
                  iterations = k, jacobian = .safe_DF(x)))

    J <- DF(x)
    if (!is.matrix(J) || nrow(J) != ncol(J) || nrow(J) != n)
      return(list(x = x, converged = FALSE, residual = res,
                  iterations = k, jacobian = J))
    if (!all(is.finite(J)))
      return(list(x = x, converged = FALSE, residual = res,
                  iterations = k, jacobian = J))

    # Regularize if near-singular
    cn <- tryCatch(kappa(J, exact = FALSE), error = function(e) Inf)
    if (!is.finite(cn) || cn > cond_max)
      J <- J + regularize * diag(n)

    # Solve J * dx = -F(x)
    dx <- tryCatch(
      solve(J, -Fx),
      error = function(e) {
        # Fallback: pseudoinverse via SVD
        sv <- svd(J)
        d_inv <- ifelse(sv$d > 1e-14, 1 / sv$d, 0)
        sv$v %*% (d_inv * crossprod(sv$u, -Fx))
      }
    )

    dx <- as.numeric(dx)
    if (length(dx) != n || !all(is.finite(dx)))
      return(list(x = x, converged = FALSE, residual = res,
                  iterations = k, jacobian = J))

    # Backtracking line search on ||F||^2
    alpha <- damping
    merit_old <- sum(Fx^2)
    accepted <- FALSE
    for (bt in seq_len(20L)) {
      x_try <- x + alpha * dx
      Fx_try <- F(x_try)
      if (is.numeric(Fx_try) && length(Fx_try) == n &&
          all(is.finite(Fx_try)) && sum(Fx_try^2) < merit_old) {
        x <- x_try
        accepted <- TRUE
        break
      }
      alpha <- 0.5 * alpha
    }
    if (!accepted) {
      return(list(x = x, converged = FALSE, residual = res,
                  iterations = k, jacobian = J))
    }
  }

  Fx_final <- F(x)
  final_res <- if (all(is.finite(Fx_final))) max(abs(Fx_final)) else Inf
  list(
    x          = x,
    converged  = is.finite(final_res) && final_res < tol,
    residual   = final_res,
    iterations = max_iter,
    jacobian   = .safe_DF(x)
  )
}


# ============================================================================
# Equilibrium index
# ============================================================================

#' Determinant-sign local equilibrium index
#'
#' Computes a determinant-sign local index
#' \eqn{\chi = \operatorname{sgn} \det D_x F_\theta(x^*)} at a given
#' regular equilibrium point. Under standard generic-game assumptions,
#' equilibrium-index sums are expected to equal \eqn{+1}.
#'
#' \describe{
#'   \item{\eqn{\chi = +1}}{Positive local orientation.}
#'   \item{\eqn{\chi = -1}}{Negative local orientation.}
#'   \item{\eqn{\chi = 0}}{Degenerate: Jacobian is singular, equilibrium
#'     is singular or numerically near singular.}
#' }
#'
#' @param J Jacobian matrix \eqn{D_x F_\theta(x^*)} at an equilibrium.
#'   Either a pre-computed matrix, or will be computed internally if
#'   \code{F} and \code{x} are supplied.
#' @param F Optional vector-valued function (used if \code{J} is
#'   \code{NULL}).
#' @param x Optional equilibrium point (used if \code{J} is \code{NULL}).
#' @param DF Optional Jacobian function (used if \code{J} is \code{NULL};
#'   falls back to \code{\link{geio_numjac}}).
#' @param tol Non-negative scalar tolerance.  The determinant is treated
#'   as zero (returning index \code{0}) when its absolute value is at or
#'   below \code{tol} times a singular-value scale factor (default
#'   \code{sqrt(.Machine$double.eps)}).
#'
#' @return Integer: \code{+1}, \code{-1}, or \code{0}.
#'
#' @references
#' Shapley, L. S. (1974). A note on the Lemke-Howson algorithm.
#' \emph{Mathematical Programming Study}, 1, 175--189.
#'
#' Govindan, S., & Wilson, R. (2003). A global Newton method to compute
#' Nash equilibria. \emph{Journal of Economic Theory}, 110(1), 65--86.
#' \doi{10.1016/S0022-0531(03)00005-X}
#'
#' For a regular real zero, the local degree reduces to
#' \eqn{\operatorname{sgn} \det J} (the Eisenbud-Khimshiashvili-Levine
#' signature), providing the local-degree background for the
#' determinant-sign diagnostic.
#'
#' @examples
#' # 1D: F(x) = x^2 - 1 at x* = 1.  F'(1) = 2 > 0, so chi = +1.
#' geio_index(J = matrix(2))
#'
#' # 1D: at x* = -1.  F'(-1) = -2 < 0, so chi = -1.
#' geio_index(J = matrix(-2))
#'
#' @export
geio_index <- function(J = NULL, F = NULL, x = NULL, DF = NULL, tol = sqrt(.Machine$double.eps)) {
  if (is.null(J)) {
    if (is.null(F) || is.null(x))
      stop("geio_index: supply either J, or both F and x")
    if (is.null(DF)) DF <- function(xx) geio_numjac(F, xx)
    J <- DF(x)
  }

  J <- as.matrix(J)
  if (!is.numeric(J)) stop("geio_index: Jacobian must be numeric")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("geio_index: tol must be a non-negative finite scalar")
  if (nrow(J) == 0L || ncol(J) == 0L || nrow(J) != ncol(J))
    stop("geio_index: Jacobian must be a non-empty square matrix")
  if (!all(is.finite(J)))
    return(0L)

  sv <- svd(J, nu = 0, nv = 0)$d
  scale <- max(1, prod(pmax(1, sv)))
  d <- det(J)
  if (!isTRUE(is.finite(d)) || abs(d) <= tol * scale) return(0L)
  as.integer(sign(d))
}


# ============================================================================
# Fragility diagnostic
# ============================================================================

#' Equilibrium fragility (minimum singular value)
#'
#' Computes \eqn{\sigma_{\min}(D_x F_\theta(x^*))}, the distance of the
#' Jacobian from singularity.  When \eqn{\sigma_{\min} \to 0}, the
#' computed equilibrium is near a singularity; this can indicate
#' loss of regularity or proximity to a bifurcation.
#'
#' @param J Jacobian matrix at an equilibrium.  If \code{NULL}, computed
#'   from \code{F} and \code{x}.
#' @param F Optional vector-valued function.
#' @param x Optional equilibrium point.
#' @param DF Optional Jacobian function.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{sigma_min}}{Scalar; minimum singular value.}
#'   \item{\code{sigma_max}}{Scalar; maximum singular value.}
#'   \item{\code{condition}}{Scalar; condition number
#'     \eqn{\sigma_{\max} / \sigma_{\min}}.}
#'   \item{\code{singular_direction}}{Numeric vector; the right singular
#'     vector corresponding to \eqn{\sigma_{\min}}, indicating the
#'     state-space direction associated with the smallest singular value.}
#' }
#'
#' @examples
#' J <- matrix(c(2, 0.1, 0.1, 0.001), 2, 2)
#' frag <- geio_fragility(J = J)
#' frag$sigma_min   # small => fragile
#'
#' @export
geio_fragility <- function(J = NULL, F = NULL, x = NULL, DF = NULL) {
  if (is.null(J)) {
    if (is.null(F) || is.null(x))
      stop("geio_fragility: supply either J, or both F and x")
    if (is.null(DF)) DF <- function(xx) geio_numjac(F, xx)
    J <- DF(x)
  }

  J <- as.matrix(J)
  if (!is.numeric(J)) stop("geio_fragility: Jacobian must be numeric")
  if (nrow(J) == 0L || ncol(J) == 0L || nrow(J) != ncol(J) || !all(is.finite(J)))
    stop("geio_fragility: Jacobian must be a non-empty finite square matrix")
  sv <- svd(J)

  sigma_min <- min(sv$d)
  sigma_max <- max(sv$d)
  cond_num  <- if (sigma_min > 0) sigma_max / sigma_min else Inf

  # Right singular vector for sigma_min
  idx_min <- which.min(sv$d)
  v_min   <- sv$v[, idx_min]

  list(
    sigma_min          = sigma_min,
    sigma_max          = sigma_max,
    condition          = cond_num,
    singular_direction = v_min
  )
}
