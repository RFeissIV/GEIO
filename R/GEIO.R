# ============================================================================
# GEIO 0.1.0 -- Game-Equilibrium Index Optimizer
# ============================================================================
# Author : Richard A. Feiss IV, Ph.D.
# Contact: feiss026@umn.edu
# Org    : Minnesota Center for Prion Research and Outreach (MNPRO)
#          University of Minnesota
# License: MIT
# ============================================================================
#
# Architecture overview
# ---------------------
# GEIO optimizes game-theoretic parameters theta while tracking the
# equilibrium structure at every iteration.  The optimization loop:
#
#   for each GALAHAD iteration k:
#     1. GALAHAD proposes theta_k
#     2. Attempt to solve F_{theta_k}(x) = 0 from multiple starts
#     3. At each equilibrium x*:
#        a. Compute chi = sgn det D_x F(x*)       [local index]
#        b. Compute sigma_min(D_x F(x*))           [fragility]
#     4. Check sum-of-indices = +1 when assumptions apply
#     5. Flag possible index changes or near-singular events
#     6. Log everything to eq_history
#
# GALAHAD (Feiss, 2026) handles the mixed-geometry optimization of theta.
# GEIO adds the game-theoretic layer on top.
#
# Key references
# --------------
#  - Shapley (1974): classical equilibrium-index background for finite games
#  - Govindan & Wilson (2003): global Newton for Nash equilibria
#  - Conn, Gould & Toint (2000): trust-region methods (via GALAHAD)
#  - Barzilai & Borwein (1988): step-size control (via GALAHAD)
# ============================================================================


#' GEIO: Game-Equilibrium Index Optimization for Finite Games
#'
#' Tracks local equilibrium-index diagnostics in simple finite-game
#' equilibrium systems.  Provides constructors for small games, a damped
#' Newton solver, determinant-sign local-index calculations,
#' singular-value near-singularity diagnostics, and conservative
#' consistency checks for selected generic examples.  Parameter
#' optimization is delegated to the \pkg{GALAHAD} package.
#'
#' @references
#' Shapley, L. S. (1974). A note on the Lemke-Howson algorithm.
#' \emph{Mathematical Programming Study}, 1, 175--189.
#'
#' Govindan, S., & Wilson, R. (2003). A global Newton method to compute
#' Nash equilibria. \emph{Journal of Economic Theory}, 110(1), 65--86.
#' \doi{10.1016/S0022-0531(03)00005-X}
#'
#' @keywords internal
"_PACKAGE"


# Item 6.5: Parts validator
#' @noRd
.geio_validate_parts <- function(parts, n_theta) {
  if (!is.list(parts)) stop("GEIO: parts must be a list")
  idx <- integer(0)
  for (nm in intersect(names(parts), c("positive", "euclidean"))) {
    val <- parts[[nm]]
    if (length(val) == 0L) next
    if (!is.numeric(val) || any(!is.finite(val)) || any(val != as.integer(val)))
      stop("GEIO: parts indices must be finite integer-like values")
    idx <- c(idx, as.integer(val))
  }
  if (length(idx) > 0L) {
    if (any(idx < 1L | idx > n_theta)) stop("GEIO: parts indices out of range")
    if (any(duplicated(idx))) stop("GEIO: parts indices must not be duplicated")
  }
  invisible(TRUE)
}

# Item 7: Safe GALAHAD wrapper
#' @noRd
.geio_run_galahad <- function(V, gradV, theta0, parts, control, callback) {
  if (!requireNamespace("GALAHAD", quietly = TRUE))
    stop("GEIO: package 'GALAHAD' is required")
  out <- GALAHAD::GALAHAD(V = V, gradV = gradV, theta0 = theta0,
                          parts = parts, control = control, callback = callback)
  required <- c("theta", "value", "converged", "reason", "iterations")
  missing <- setdiff(required, names(out))
  if (length(missing) > 0L)
    stop("GEIO: GALAHAD output missing field(s): ", paste(missing, collapse = ", "))
  if (!is.numeric(out$theta) || length(out$theta) != length(theta0) || !all(is.finite(out$theta)))
    stop("GEIO: GALAHAD returned invalid theta")
  out
}


#' GEIO: Game-Equilibrium Index Optimizer
#'
#' @description
#' Optimizes game parameters \eqn{\theta} by minimizing an objective
#' \eqn{V(\theta)} while tracking the equilibrium structure of a
#' parameterized game at every iteration.  At each parameter iterate,
#' GEIO attempts to solve the equilibrium system \eqn{F_\theta(x) = 0}
#' from multiple starts, computes the local index
#' \eqn{\chi = \operatorname{sgn} \det D_x F_\theta(x^*)} at each one,
#' and monitors fragility via \eqn{\sigma_{\min}(D_x F)}.
#'
#' The optimizer is powered by \code{\link[GALAHAD]{GALAHAD}}, a
#' geometry-adaptive trust-region method with softplus
#' reparameterization.  Parameters declared as \strong{positive} (e.g.,
#' payoffs that must be non-negative) are handled via softplus in
#' z-space; \strong{euclidean} parameters are unconstrained.
#'
#' @section Equilibrium tracking:
#' At each GALAHAD iteration, the equilibrium tracker:
#' \enumerate{
#'   \item Uses multi-start Newton from pure-strategy corners and a
#'     user-supplied or default interior guess to recover equilibria.
#'   \item Deduplicates solutions within tolerance \code{eq_tol}.
#'   \item Computes the determinant-sign local index at each equilibrium.
#'   \item Checks the sum-of-indices identity
#'     (\eqn{\sum \chi_i = +1}) when the generic-game assumptions apply.
#'   \item Computes \eqn{\sigma_{\min}} fragility at each equilibrium.
#'   \item Detects bifurcation events (index sign changes, equilibrium
#'     creation/destruction between iterations).
#' }
#'
#' @param game A \code{geio_game} object from \code{\link{geio_game}},
#'   or a named list with elements \code{F_theta}, \code{DF_theta},
#'   \code{n_x}, and optionally \code{x0_default}, \code{all_pure}.
#' @param V Objective function: \code{function(theta) -> scalar}.
#' @param gradV Gradient of \code{V}: \code{function(theta) -> numeric
#'   vector}.  Use \code{GALAHAD::galahad_numgrad} if analytical
#'   gradient is unavailable.
#' @param theta0 Initial parameter vector.
#' @param parts GALAHAD-style parameter partition.  See
#'   \code{\link[GALAHAD]{galahad_parts}}.
#' @param x0 Optional initial equilibrium guess.  If \code{NULL}, the
#'   game's \code{x0_default(theta0)} is used.
#' @param control Named list of control parameters.  Includes all
#'   GALAHAD control parameters plus GEIO-specific ones:
#'   \describe{
#'     \item{\code{eq_tol}}{Tolerance for declaring two equilibria
#'       identical (default \code{1e-6}).}
#'     \item{\code{newton_max}}{Maximum Newton iterations per
#'       equilibrium solve (default \code{200}).}
#'     \item{\code{newton_tol}}{Newton convergence tolerance
#'       (default \code{1e-10}).}
#'     \item{\code{bifurcation_tol}}{Threshold on
#'       \eqn{\sigma_{\min}} below which a bifurcation warning is
#'       issued (default \code{1e-4}).}
#'     \item{\code{track_every}}{Track equilibria every \code{n}
#'       iterations (default \code{1}).  Set higher to reduce
#'       computational cost.}
#'   }
#' @param callback Optional user callback, called after equilibrium
#'   tracking at each iteration.  Receives a list with \code{iter},
#'   \code{theta}, \code{value}, \code{equilibria}, \code{indices},
#'   \code{fragilities}, \code{bifurcation}.
#'
#' @return A named list with class \code{"geio_result"}:
#' \describe{
#'   \item{\code{theta}}{Final parameter estimate.}
#'   \item{\code{value}}{Final objective value.}
#'   \item{\code{converged}}{Logical; optimizer convergence.}
#'   \item{\code{reason}}{Convergence reason string.}
#'   \item{\code{iterations}}{Number of GALAHAD iterations.}
#'   \item{\code{galahad}}{Full GALAHAD output object.}
#'   \item{\code{equilibria}}{List of equilibrium points at the final
#'     \eqn{\theta}.}
#'   \item{\code{indices}}{Integer vector of determinant-sign local indices at
#'     final equilibria.}
#'   \item{\code{fragilities}}{List of fragility diagnostics at
#'     final equilibria.}
#'   \item{\code{index_sum}}{Sum of indices at final \eqn{\theta}
#'     (should be +1 for generic games).}
#'   \item{\code{bifurcations}}{Data frame of detected bifurcation
#'     events: iteration, type (creation/destruction/flip), and
#'     \eqn{\sigma_{\min}} at the event.}
#'   \item{\code{eq_history}}{List of per-iteration equilibrium
#'     tracking records.}
#'   \item{\code{certificate}}{Convergence and index-sum certificate.}
#' }
#'
#' @examples
#' # Stag Hunt: find the S,H that minimize deviation from target
#' # mixed-strategy equilibrium at p* = 0.6
#' game <- geio_game("stag_hunt")
#'
#' target_p <- 0.6
#' V <- function(theta) {
#'   S <- theta[1]; H <- theta[2]
#'   p_star <- H / S   # mixed eq for Stag Hunt
#'   (p_star - target_p)^2
#' }
#' gV <- function(theta) GALAHAD::galahad_numgrad(V, theta)
#'
#' \dontrun{
#' fit <- GEIO(game   = game,
#'             V      = V,
#'             gradV  = gV,
#'             theta0 = c(S = 4, H = 1),
#'             parts  = list(positive = c(1L, 2L),
#'                           euclidean = integer(0)),
#'             control = list(max_iter = 100))
#'
#' fit$theta          # optimized (S, H)
#' fit$indices        # determinant-sign local indices at final equilibria
#' fit$index_sum      # local-index sum when applicable
#' fit$bifurcations   # any bifurcation events during optimization
#' }
#'
#' @references
#' Feiss, R. A. (2026). GALAHAD: Geometry-Adaptive Lyapunov-Assured
#' Hybrid Optimizer with Softplus Reparameterization and Trust-Region
#' Control. R package version 2.0.0.
#' \url{https://CRAN.R-project.org/package=GALAHAD}
#'
#' Shapley, L. S. (1974). A note on the Lemke-Howson algorithm.
#' \emph{Mathematical Programming Study}, 1, 175--189.
#'
#' Govindan, S., & Wilson, R. (2003). A global Newton method to compute
#' Nash equilibria. \emph{Journal of Economic Theory}, 110(1), 65--86.
#' \doi{10.1016/S0022-0531(03)00005-X}
#'
#' Conn, A. R., Gould, N. I. M., & Toint, P. L. (2000).
#' \emph{Trust-Region Methods}. SIAM.
#' \doi{10.1137/1.9780898719857}
#'
#' @seealso \code{\link{geio_game}} for game constructors;
#'   \code{\link{geio_index}} for standalone index computation;
#'   \code{\link{geio_fragility}} for standalone fragility;
#'   \code{\link{geio_newton}} for the equilibrium solver;
#'   \code{\link[GALAHAD]{GALAHAD}} for the underlying optimizer.
#'
#' @export
GEIO <- function(game, V, gradV, theta0, parts, x0 = NULL,
                 control = list(), callback = NULL) {

  # --- Validate game object --------------------------------------------------
  if (!is.list(game) || is.null(game$F_theta))
    stop("GEIO: 'game' must be a geio_game object or list with F_theta")

  n_x <- game$n_x
  if (is.null(n_x) || length(n_x) != 1L || !is.finite(n_x) || n_x < 1)
    stop("GEIO: game$n_x must be a positive finite integer")
  n_x <- as.integer(n_x)

  if (!is.function(V))
    stop("GEIO: V must be a function")
  if (!is.function(gradV))
    stop("GEIO: gradV must be a function")
  if (!is.numeric(theta0) || !all(is.finite(theta0)))
    stop("GEIO: theta0 must be a finite numeric vector")
  n_theta <- game$n_theta
  if (!is.null(n_theta)) {
    if (length(n_theta) != 1L || !is.finite(n_theta) || n_theta < 1)
      stop("GEIO: game$n_theta must be a positive finite integer")
    n_theta <- as.integer(n_theta)
    if (length(theta0) != n_theta)
      stop("GEIO: length(theta0) must equal game$n_theta")
  }
  # Preflight V and gradV at theta0
  v0 <- V(theta0)
  if (!is.numeric(v0) || length(v0) != 1L || !is.finite(v0))
    stop("GEIO: V(theta0) must return a finite numeric scalar")
  g0 <- gradV(theta0)
  if (!is.numeric(g0) || length(g0) != length(theta0) || !all(is.finite(g0)))
    stop("GEIO: gradV(theta0) must return a finite numeric vector matching theta0")
  if (!is.list(parts))
    stop("GEIO: parts must be a GALAHAD partition list")
  .geio_validate_parts(parts, length(theta0))

  F_theta  <- game$F_theta
  DF_theta <- game$DF_theta
  is_valid <- game$is_valid
  if (is.null(is_valid)) {
    is_valid <- function(theta, x, tol = 1e-8) all(is.finite(x))
  }

  # If DF_theta is NULL, use numerical Jacobian

  if (is.null(DF_theta))
    DF_theta <- function(theta, x) geio_numjac(function(xx) F_theta(theta, xx), x)

  # --- GEIO-specific control defaults ----------------------------------------
  if (!is.list(control)) stop("GEIO: control must be a list")
  eq_tol          <- if (!is.null(control$eq_tol))          control$eq_tol          else 1e-6
  newton_max      <- if (!is.null(control$newton_max))      control$newton_max      else 200L
  newton_tol      <- if (!is.null(control$newton_tol))      control$newton_tol      else 1e-10
  bifurcation_tol <- if (!is.null(control$bifurcation_tol)) control$bifurcation_tol else 1e-4
  track_every     <- if (!is.null(control$track_every))     control$track_every     else 1L

  if (!is.numeric(eq_tol) || length(eq_tol) != 1L || !is.finite(eq_tol) || eq_tol <= 0)
    stop("GEIO: control$eq_tol must be a positive finite scalar")
  if (!is.numeric(newton_max) || length(newton_max) != 1L || !is.finite(newton_max) || newton_max < 1)
    stop("GEIO: control$newton_max must be a positive integer-like scalar")
  if (!is.numeric(newton_tol) || length(newton_tol) != 1L || !is.finite(newton_tol) || newton_tol <= 0)
    stop("GEIO: control$newton_tol must be a positive finite scalar")
  if (!is.numeric(bifurcation_tol) || length(bifurcation_tol) != 1L || !is.finite(bifurcation_tol) || bifurcation_tol < 0)
    stop("GEIO: control$bifurcation_tol must be a non-negative finite scalar")
  if (!is.numeric(track_every) || length(track_every) != 1L || !is.finite(track_every) || track_every < 1)
    stop("GEIO: control$track_every must be a positive integer-like scalar")
  newton_max  <- as.integer(newton_max)
  track_every <- as.integer(track_every)

  # Strip GEIO-specific keys before passing to GALAHAD
  if (!is.null(callback) && !is.function(callback)) stop("GEIO: callback must be NULL or a function")
  galahad_ctrl <- control
  galahad_ctrl$eq_tol          <- NULL
  galahad_ctrl$newton_max      <- NULL
  galahad_ctrl$newton_tol      <- NULL
  galahad_ctrl$bifurcation_tol <- NULL
  galahad_ctrl$track_every     <- NULL

  # --- Initialize equilibrium tracking state ---------------------------------
  eq_env <- new.env(parent = emptyenv())
  eq_env$history      <- list()
  eq_env$prev_n_eq    <- NA_integer_
  eq_env$prev_indices <- integer(0)
  eq_env$prev_equilibria <- list()
  eq_env$bifurcations <- data.frame(
    iter       = integer(0),
    type       = character(0),
    sigma_min  = numeric(0),
    n_eq_before = integer(0),
    n_eq_after  = integer(0),
    stringsAsFactors = FALSE
  )

  # --- Build the GALAHAD callback that does equilibrium tracking -------------
  # Safe fragility wrapper (item 8.6)
  .geio_safe_fragility <- function(J, nx) {
    tryCatch(
      geio_fragility(J = J),
      error = function(e) list(
        sigma_min = NA_real_, sigma_max = NA_real_,
        condition = NA_real_,
        singular_direction = rep(NA_real_, nx)
      )
    )
  }

  geio_callback <- function(info) {
    if (is.null(info$iter) || is.null(info$theta)) return(invisible(NULL))
    k     <- as.integer(info$iter)
    theta <- as.numeric(info$theta)
    if (!is.finite(k) || !all(is.finite(theta))) return(invisible(NULL))
    value <- if (!is.null(info$value)) info$value else NA_real_

    if (k %% track_every != 0L && k > 1L) return(invisible(NULL))

    # Multi-start Newton: try from pure corners and interior guess
    starts <- list()
    if (!is.null(game$all_pure)) starts <- game$all_pure(theta)
    if (!is.null(x0))           starts <- c(starts, list(x0))
    if (!is.null(game$x0_default))
      starts <- c(starts, list(game$x0_default(theta)))
    if (length(starts) == 0L && !is.null(n_x))
      starts <- list(rep(0.5, n_x))

    # Solve from each start
    raw_eqs <- list()
    for (s in starts) {
      sol <- geio_newton(
        F       = function(x) F_theta(theta, x),
        x0      = s,
        DF      = function(x) DF_theta(theta, x),
        max_iter = newton_max,
        tol      = newton_tol
      )
      if (sol$converged && is_valid(theta, sol$x, eq_tol)) {
        raw_eqs <- c(raw_eqs, list(sol))
      }
    }

    # Deduplicate interior equilibria found by Newton
    unique_eqs <- .geio_deduplicate(raw_eqs, eq_tol)

    # Compute index and fragility at each interior equilibrium
    indices     <- integer(length(unique_eqs))
    fragilities <- vector("list", length(unique_eqs))

    for (i in seq_along(unique_eqs)) {
      J <- unique_eqs[[i]]$jacobian
      indices[i]     <- geio_index(J = J)
      fragilities[[i]] <- .geio_safe_fragility(J, n_x)
    }

    # --- Merge verified pure equilibria (boundary) ---------------------------
    # The indifference system F_theta(x) = 0 only captures interior (fully
    # mixed) equilibria.  Pure NE are best-response fixed points, not zeros
    # of F.  We add them from game$pure_equilibria() with pre-computed
    # determinant-sign local indices (+1 for strict, 0 for degenerate).
    has_pure <- !is.null(game$pure_equilibria)
    if (has_pure) {
      pure_eqs <- game$pure_equilibria(theta)
      for (pe in pure_eqs) {
        # Skip if Newton already found something at this point
        is_dup <- FALSE
        for (ue in unique_eqs) {
          if (max(abs(pe$x - ue$x)) < eq_tol) { is_dup <- TRUE; break }
        }
        if (!is_dup) {
          unique_eqs <- c(unique_eqs, list(list(
            x = pe$x, converged = TRUE, residual = 0,
            iterations = 0L, jacobian = matrix(NA_real_, n_x, n_x)
          )))
          indices     <- c(indices, pe$index)
          fragilities <- c(fragilities, list(list(
            sigma_min = NA_real_, sigma_max = NA_real_,
            condition = NA_real_, singular_direction = rep(NA_real_, n_x)
          )))
        }
      }
    }

    idx_sum <- sum(indices)

    # --- Bifurcation/proximity diagnostics -----------------------------------------------
    n_eq <- length(unique_eqs)

    if (!is.na(eq_env$prev_n_eq)) {
      # Equilibrium count changed
      if (n_eq != eq_env$prev_n_eq) {
        bif_type <- if (n_eq > eq_env$prev_n_eq) "creation" else "destruction"
        sigma_vals <- if (length(fragilities) > 0)
          vapply(fragilities, function(f) f$sigma_min, numeric(1))
        else numeric(0)
        sigma_vals <- sigma_vals[is.finite(sigma_vals)]
        min_sigma <- if (length(sigma_vals) > 0L) min(sigma_vals) else NA_real_
        eq_env$bifurcations <- rbind(eq_env$bifurcations, data.frame(
          iter        = k,
          type        = bif_type,
          sigma_min   = min_sigma,
          n_eq_before = eq_env$prev_n_eq,
          n_eq_after  = n_eq,
          stringsAsFactors = FALSE
        ))
      }

      # Index sign flip among proximity-matched equilibria.  Equilibria can
      # be returned in different list orders across iterations, so comparing
      # raw positions can create false flip events.
      matched <- .geio_match_equilibria(eq_env$prev_equilibria, unique_eqs, eq_tol)
      if (nrow(matched) > 0L && length(eq_env$prev_indices) >= max(matched$old)) {
        old_idx <- eq_env$prev_indices[matched$old]
        new_idx <- indices[matched$new]
        flipped <- matched$new[which(old_idx != new_idx)]
        if (length(flipped) > 0L) {
          sigma_vals <- vapply(
            fragilities[flipped], function(f) f$sigma_min, numeric(1))
          sigma_vals <- sigma_vals[is.finite(sigma_vals)]
          min_sigma <- if (length(sigma_vals) > 0L) min(sigma_vals) else NA_real_
          eq_env$bifurcations <- rbind(eq_env$bifurcations, data.frame(
            iter        = k,
            type        = "index_flip",
            sigma_min   = min_sigma,
            n_eq_before = n_eq,
            n_eq_after  = n_eq,
            stringsAsFactors = FALSE
          ))
        }
      }

      # Near-singularity fragility warning based on GEIO-specific threshold.
      sigma_vals <- if (length(fragilities) > 0L)
        vapply(fragilities, function(f) f$sigma_min, numeric(1))
      else numeric(0)
      sigma_vals <- sigma_vals[is.finite(sigma_vals)]
      if (length(sigma_vals) > 0L && min(sigma_vals) < bifurcation_tol) {
        eq_env$bifurcations <- rbind(eq_env$bifurcations, data.frame(
          iter        = k,
          type        = "near_singular",
          sigma_min   = min(sigma_vals),
          n_eq_before = n_eq,
          n_eq_after  = n_eq,
          stringsAsFactors = FALSE
        ))
      }
    }

    eq_env$prev_n_eq    <- n_eq
    eq_env$prev_indices <- indices
  eq_env$prev_equilibria <- unique_eqs

    # --- Store history record ------------------------------------------------
    record <- list(
      iter        = k,
      theta       = theta,
      n_eq        = n_eq,
      equilibria  = lapply(unique_eqs, function(e) e$x),
      indices     = indices,
      index_sum   = idx_sum,
      fragilities = fragilities,
      sigma_mins  = vapply(fragilities, function(f) f$sigma_min, numeric(1))
    )
    eq_env$history[[length(eq_env$history) + 1L]] <- record

    # --- User callback -------------------------------------------------------
    if (!is.null(callback)) {
      tryCatch(
        callback(list(
          iter        = k,
          theta       = theta,
          value       = value,
          equilibria  = record$equilibria,
          indices     = indices,
          fragilities = fragilities,
          bifurcation = nrow(eq_env$bifurcations) > 0 &&
                        eq_env$bifurcations$iter[nrow(eq_env$bifurcations)] == k
        )),
        error = function(e) invisible(NULL)
      )
    }

    invisible(NULL)
  }

  # --- Run GALAHAD -----------------------------------------------------------
  gfit <- .geio_run_galahad(V, gradV, theta0, parts, galahad_ctrl, geio_callback)

  # --- Final equilibrium classification at optimized theta -------------------
  theta_final <- gfit$theta
  starts_final <- list()
  if (!is.null(game$all_pure)) starts_final <- game$all_pure(theta_final)
  if (!is.null(x0))            starts_final <- c(starts_final, list(x0))
  if (!is.null(game$x0_default))
    starts_final <- c(starts_final, list(game$x0_default(theta_final)))
  if (length(starts_final) == 0L && !is.null(n_x))
    starts_final <- list(rep(0.5, n_x))

  raw_final <- list()
  for (s in starts_final) {
    sol <- geio_newton(
      F        = function(x) F_theta(theta_final, x),
      x0       = s,
      DF       = function(x) DF_theta(theta_final, x),
      max_iter = newton_max,
      tol      = newton_tol
    )
    if (sol$converged && is_valid(theta_final, sol$x, eq_tol)) {
      raw_final <- c(raw_final, list(sol))
    }
  }
  unique_final <- .geio_deduplicate(raw_final, eq_tol)

  final_indices     <- integer(length(unique_final))
  final_fragilities <- vector("list", length(unique_final))
  for (i in seq_along(unique_final)) {
    J <- unique_final[[i]]$jacobian
    final_indices[i]     <- geio_index(J = J)
    final_fragilities[[i]] <- .geio_safe_fragility(J, n_x)
  }

  # Merge verified pure equilibria at final theta
  has_pure <- !is.null(game$pure_equilibria)
  if (has_pure) {
    pure_final <- game$pure_equilibria(theta_final)
    for (pe in pure_final) {
      is_dup <- FALSE
      for (ue in unique_final) {
        if (max(abs(pe$x - ue$x)) < eq_tol) { is_dup <- TRUE; break }
      }
      if (!is_dup) {
        unique_final <- c(unique_final, list(list(
          x = pe$x, converged = TRUE, residual = 0,
          iterations = 0L, jacobian = matrix(NA_real_, n_x, n_x)
        )))
        final_indices     <- c(final_indices, pe$index)
        final_fragilities <- c(final_fragilities, list(list(
          sigma_min = NA_real_, sigma_max = NA_real_,
          condition = NA_real_, singular_direction = rep(NA_real_, n_x)
        )))
      }
    }
  }

  # --- Assemble result -------------------------------------------------------
  idx_sum_final <- sum(final_indices)
  expected_n <- if (game$type %in% c("stag_hunt", "bos")) 3L else NA_integer_
  index_sum_applicable <- isTRUE(has_pure) &&
    game$type %in% c("stag_hunt", "bos") &&
    .geio_generic_assumptions_hold(game, theta_final) &&
    all(final_indices != 0L) &&
    length(unique_final) > 0L &&
    (is.na(expected_n) || length(unique_final) == expected_n)

  result <- list(
    theta        = theta_final,
    value        = gfit$value,
    converged    = gfit$converged,
    reason       = gfit$reason,
    iterations   = gfit$iterations,
    galahad      = gfit,
    equilibria   = lapply(unique_final, function(e) e$x),
    indices      = final_indices,
    fragilities  = final_fragilities,
    index_sum    = idx_sum_final,
    index_sum_applicable = index_sum_applicable,
    bifurcations = eq_env$bifurcations,
    eq_history   = eq_env$history,
    certificate  = list(
      converged      = gfit$converged,
      reason         = gfit$reason,
      index_sum      = idx_sum_final,
      index_sum_ok   = if (index_sum_applicable) idx_sum_final == 1L else NA,
      index_sum_applicable = index_sum_applicable,
      has_pure       = has_pure,
      n_equilibria   = length(unique_final),
      n_bifurcations = nrow(eq_env$bifurcations)
    )
  )

  class(result) <- "geio_result"
  result
}


# ============================================================================
# Print method
# ============================================================================

# Item 8.9: Check generic-game assumptions
.geio_generic_assumptions_hold <- function(game, theta) {
  if (is.null(game$type)) return(FALSE)
  if (game$type == "stag_hunt")
    return(length(theta) == 2L && all(is.finite(theta)) && theta[1] > theta[2] && theta[2] > 0)
  if (game$type == "bos")
    return(length(theta) == 2L && all(is.finite(theta)) && theta[1] > 0 && theta[2] > 0)
  FALSE
}

#' @export
print.geio_result <- function(x, ...) {
  cat("GEIO: Game-Equilibrium Index Optimizer\n")
  cat("--------------------------------------\n")
  cat(sprintf("  Converged:       %s (%s)\n", x$converged, x$reason))
  cat(sprintf("  Iterations:      %d\n", x$iterations))
  cat(sprintf("  Objective:       %.6e\n", x$value))
  cat(sprintf("  Theta:           %s\n",
              paste(sprintf("%.4f", x$theta), collapse = ", ")))
  cat(sprintf("  Equilibria:      %d\n", length(x$equilibria)))
  if (isTRUE(x$index_sum_applicable)) {
    cat(sprintf("  Index sum:       %d", x$index_sum))
    if (identical(x$index_sum, 1L)) cat("  [OK: applicable built-in index check]\n")
    else cat("  [CHECK: expected +1 only when all regular equilibria were recovered]\n")
  } else {
    cat("  Index sum:       not evaluated for this custom/non-enumerated game\n")
  }
  cat(sprintf("  Bifurcation diagnostics:    %d\n", nrow(x$bifurcations)))

  if (length(x$equilibria) > 0) {
    cat("\n  Equilibria at final theta:\n")
    for (i in seq_along(x$equilibria)) {
      xi  <- x$equilibria[[i]]
      chi <- x$indices[i]
      sig <- x$fragilities[[i]]$sigma_min
      is_pure <- is.na(sig)
      stab <- if (chi == 1L) "INDEX +1" else if (chi == -1L) "INDEX -1" else "DEGENERATE"
      if (is_pure) {
        cat(sprintf("    [%d] x* = (%s)  chi = %+d  [%s, PURE]\n",
                    i, paste(sprintf("%.4f", xi), collapse = ", "),
                    chi, stab))
      } else {
        cat(sprintf("    [%d] x* = (%s)  chi = %+d  sigma_min = %.2e  [%s]\n",
                    i, paste(sprintf("%.4f", xi), collapse = ", "),
                    chi, sig, stab))
      }
    }
  }

  invisible(x)
}


# ============================================================================
# Internal: deduplication
# ============================================================================

#' @keywords internal
.geio_deduplicate <- function(eq_list, tol) {
  if (length(eq_list) == 0L) return(list())

  unique_eqs <- list(eq_list[[1]])

  for (i in seq_along(eq_list)[-1]) {
    is_dup <- FALSE
    for (j in seq_along(unique_eqs)) {
      if (max(abs(eq_list[[i]]$x - unique_eqs[[j]]$x)) < tol) {
        # Keep the one with smaller residual
        if (eq_list[[i]]$residual < unique_eqs[[j]]$residual)
          unique_eqs[[j]] <- eq_list[[i]]
        is_dup <- TRUE
        break
      }
    }
    if (!is_dup) unique_eqs <- c(unique_eqs, list(eq_list[[i]]))
  }

  unique_eqs
}

#' @keywords internal
.geio_match_equilibria <- function(old_eqs, new_eqs, tol) {
  if (length(old_eqs) == 0L || length(new_eqs) == 0L)
    return(data.frame(old = integer(0), new = integer(0), distance = numeric(0)))

  old_x <- lapply(old_eqs, function(e) e$x)
  new_x <- lapply(new_eqs, function(e) e$x)
  pairs <- data.frame(old = integer(0), new = integer(0), distance = numeric(0))

  for (i in seq_along(old_x)) {
    d <- vapply(new_x, function(xn) max(abs(old_x[[i]] - xn)), numeric(1))
    j <- which.min(d)
    if (length(j) == 1L && is.finite(d[j]) && d[j] <= tol) {
      pairs <- rbind(pairs, data.frame(old = i, new = j, distance = d[j]))
    }
  }

  if (nrow(pairs) == 0L) return(pairs)
  pairs <- pairs[order(pairs$distance), , drop = FALSE]
  pairs <- pairs[!duplicated(pairs$old) & !duplicated(pairs$new), , drop = FALSE]
  rownames(pairs) <- NULL
  pairs
}
