.geio_assert_theta <- function(theta, n, label) {
  if (!is.numeric(theta) || length(theta) != n || !all(is.finite(theta)))
    stop(label, ": theta must be a finite numeric vector of length ", n)
  as.numeric(theta)
}
.geio_assert_x <- function(x, n, label) {
  if (!is.numeric(x) || length(x) != n || !all(is.finite(x)))
    stop(label, ": x must be a finite numeric vector of length ", n)
  as.numeric(x)
}

# ============================================================================
# GEIO -- Built-in Game Definitions
# ============================================================================
# Author : Richard A. Feiss IV, Ph.D.
# Contact: feiss026@umn.edu
# Org    : Minnesota Center for Prion Research and Outreach (MNPRO)
#          University of Minnesota
# License: MIT
# ============================================================================
#
# Provides parameterized game definitions that return the equilibrium
# system F_theta(x) = 0 and its Jacobian D_x F_theta(x) for common
# game classes.
#
# For a 2x2 bimatrix game with payoff matrices A (row player) and
# B (column player), the mixed-strategy Nash equilibrium system is:
#
#   F(p, q) = ( (A - A') q,  (B - B')^T p )
#
# where p and q are mixed strategies (probabilities of playing action 1).
# In 2x2, this reduces to a 2D system in (p, q) in [0,1]^2.
#
# The parameterization theta encodes payoff entries, allowing GEIO
# to optimize over payoff perturbations while tracking equilibrium
# structure.
# ============================================================================


#' Construct a parameterized game for GEIO
#'
#' Creates a game specification object that GEIO can optimize over.
#' The game is defined by its equilibrium system \eqn{F_\theta(x) = 0}
#' and Jacobian \eqn{D_x F_\theta(x)}, parameterized by \eqn{\theta}.
#'
#' Built-in game types:
#' \describe{
#'   \item{\code{"bimatrix"}}{General 2x2 bimatrix game.  \code{theta}
#'     contains the 8 payoff entries
#'     \code{c(a11, a12, a21, a22, b11, b12, b21, b22)}.}
#'   \item{\code{"stag_hunt"}}{Parameterized Stag Hunt.
#'     \code{theta = c(stag_reward, hare_reward)} with
#'     \eqn{\text{stag\_reward} > \text{hare\_reward} > 0}.}
#'   \item{\code{"bos"}}{Battle of the Sexes.
#'     \code{theta = c(alpha, beta)} where \eqn{\alpha} is the
#'     coordination payoff for player 1's preferred outcome and
#'     \eqn{\beta} for player 2's.}
#'   \item{\code{"custom"}}{User-supplied \code{F_theta}, \code{DF_theta},
#'     and equilibrium dimension \code{n_x}.}
#' }
#'
#' @param type Character; one of \code{"bimatrix"}, \code{"stag_hunt"},
#'   \code{"bos"}, or \code{"custom"}.
#' @param F_theta For \code{"custom"}: function(theta, x) -> numeric
#'   vector.
#' @param DF_theta For \code{"custom"}: function(theta, x) -> matrix
#'   (Jacobian w.r.t. x).  If \code{NULL}, numerical Jacobian used.
#' @param n_x For \code{"custom"}: integer dimension of the equilibrium
#'   variable x.
#' @param n_theta For \code{"custom"}: integer dimension of the parameter
#'   vector theta.
#'
#' @return A named list with class \code{"geio_game"}:
#' \describe{
#'   \item{\code{type}}{Character string identifying the game type.}
#'   \item{\code{F_theta}}{Function(theta, x) -> numeric vector.}
#'   \item{\code{DF_theta}}{Function(theta, x) -> Jacobian matrix D_x F.}
#'   \item{\code{n_x}}{Integer; dimension of equilibrium space.}
#'   \item{\code{n_theta}}{Integer; dimension of parameter space.}
#'   \item{\code{x0_default}}{Function(theta) -> reasonable starting
#'     guess for Newton.}
#'   \item{\code{all_pure}}{Function(theta) -> list of pure-strategy
#'     equilibria (for multi-start initialization).}
#'   \item{\code{pure_equilibria}}{Function(theta) -> list of verified
#'     pure Nash equilibria, each a list with elements \code{x} (strategy
#'     profile) and \code{index} (+1 for strict, 0 for degenerate).
#'     \code{NULL} for custom games (user must supply separately).
#'     Required for correct sum-of-indices computation, since the
#'     indifference system \eqn{F_\theta(x) = 0} only captures interior
#'     (fully mixed) equilibria.}
#' }
#'
#' @examples
#' # Stag Hunt with default parameterization
#' game <- geio_game("stag_hunt")
#' # F_theta at theta = (4, 2), mixed strategy (0.5, 0.5)
#' game$F_theta(c(4, 2), c(0.5, 0.5))
#'
#' @export
geio_game <- function(type = c("bimatrix", "stag_hunt", "bos", "custom"),
                       F_theta = NULL, DF_theta = NULL,
                       n_x = NULL, n_theta = NULL) {

  type <- match.arg(type)

  game <- switch(type,
    bimatrix   = .geio_bimatrix(),
    stag_hunt  = .geio_stag_hunt(),
    bos        = .geio_bos(),
    custom     = .geio_custom(F_theta, DF_theta, n_x, n_theta)
  )

  class(game) <- "geio_game"
  game
}


# ============================================================================
# 2x2 Bimatrix game
# ============================================================================
#
# Payoff matrices:
#   A = [[a11, a12], [a21, a22]]   (row player)
#   B = [[b11, b12], [b21, b22]]   (column player)
#
# theta = c(a11, a12, a21, a22, b11, b12, b21, b22)
#
# Mixed strategy: p = P(row plays 1), q = P(col plays 1)
#
# Indifference conditions (equilibrium system):
# Player 1 (row) chooses p = prob of action 1.
# Player 2 (col) chooses q = prob of action 1.
#
# Player 1's expected payoff from action 1: q*a11 + (1-q)*a12
# Player 1's expected payoff from action 2: q*a21 + (1-q)*a22
# Player 1 is indifferent when:
#   q*(a11 - a21) + (1-q)*(a12 - a22) = 0
#   => q*(a11 - a21 - a12 + a22) + (a12 - a22) = 0
#
# Player 2's expected payoff from action 1: p*b11 + (1-p)*b21
# Player 2's expected payoff from action 2: p*b12 + (1-p)*b22
# Player 2 is indifferent when:
#   p*(b11 - b12) + (1-p)*(b21 - b22) = 0
#   => p*(b11 - b12 - b21 + b22) + (b21 - b22) = 0
#
# So F(p,q) maps R^2 -> R^2.
# ============================================================================

#' @keywords internal
.geio_bimatrix <- function() {
  list(
    type    = "bimatrix",
    n_x     = 2L,
    n_theta = 8L,

    F_theta = function(theta, x) {
      .geio_assert_theta(theta, 8L, "bimatrix")
      .geio_assert_x(x, 2L, "bimatrix")
      a11 <- theta[1]; a12 <- theta[2]; a21 <- theta[3]; a22 <- theta[4]
      b11 <- theta[5]; b12 <- theta[6]; b21 <- theta[7]; b22 <- theta[8]
      p <- x[1]; q <- x[2]

      # Player 1 indifference (determines q at interior eq)
      f1 <- q * (a11 - a21 - a12 + a22) + (a12 - a22)
      # Player 2 indifference (determines p at interior eq)
      f2 <- p * (b11 - b12 - b21 + b22) + (b21 - b22)

      c(f1, f2)
    },

    DF_theta = function(theta, x) {
      .geio_assert_theta(theta, 8L, "bimatrix")
      .geio_assert_x(x, 2L, "bimatrix")
      a11 <- theta[1]; a12 <- theta[2]; a21 <- theta[3]; a22 <- theta[4]
      b11 <- theta[5]; b12 <- theta[6]; b21 <- theta[7]; b22 <- theta[8]
      p <- x[1]; q <- x[2]

      # df1/dp = 0  (f1 doesn't depend on p)
      # df1/dq = a11 - a21 - a12 + a22
      # df2/dp = b11 - b12 - b21 + b22
      # df2/dq = 0  (f2 doesn't depend on q)
      matrix(c(0, b11 - b12 - b21 + b22,
               a11 - a21 - a12 + a22, 0), nrow = 2, ncol = 2)
    },

    x0_default = function(theta) c(0.5, 0.5),

    all_pure = function(theta) {
      list(c(0, 0), c(0, 1), c(1, 0), c(1, 1))
    },

    # Verified pure Nash equilibria with determinant-sign local indices.
    # A strict pure NE has index +1; a degenerate one (player indifferent)
    # has index 0.  Convention 1: B[i][j] = col payoff when row=i, col=j.
    is_valid = function(theta, x, tol = 1e-8) {
      length(x) == 2L && all(is.finite(x)) && all(x >= -tol) && all(x <= 1 + tol)
    },

    pure_equilibria = function(theta) {
      a11 <- theta[1]; a12 <- theta[2]; a21 <- theta[3]; a22 <- theta[4]
      b11 <- theta[5]; b12 <- theta[6]; b21 <- theta[7]; b22 <- theta[8]
      eqs <- list()

      # (1,1): row best-responds iff a11 >= a21; col iff b11 >= b12
      if (a11 >= a21 && b11 >= b12)
        eqs <- c(eqs, list(list(
          x = c(1, 1),
          index = if (a11 > a21 && b11 > b12) 1L else 0L)))

      # (1,0): row iff a12 >= a22; col iff b12 >= b11
      if (a12 >= a22 && b12 >= b11)
        eqs <- c(eqs, list(list(
          x = c(1, 0),
          index = if (a12 > a22 && b12 > b11) 1L else 0L)))

      # (0,1): row iff a21 >= a11; col iff b21 >= b22
      if (a21 >= a11 && b21 >= b22)
        eqs <- c(eqs, list(list(
          x = c(0, 1),
          index = if (a21 > a11 && b21 > b22) 1L else 0L)))

      # (0,0): row iff a22 >= a12; col iff b22 >= b21
      if (a22 >= a12 && b22 >= b21)
        eqs <- c(eqs, list(list(
          x = c(0, 0),
          index = if (a22 > a12 && b22 > b21) 1L else 0L)))

      eqs
    }
  )
}


# ============================================================================
# Stag Hunt
# ============================================================================
#
# Parameterization: theta = c(S, H)
#   S = stag_reward (both cooperate)
#   H = hare_reward (defect alone)
#   Convention: S > H > 0
#
# Payoff matrices (Convention 1: M[i][j] = player's payoff when
#   row plays action i, col plays action j):
#   A (row) = [[S, 0], [H, H]]
#   B (col) = [[S, H], [0, H]]  (= A^T by symmetry)
#
# The indifference system uses the symmetric structure directly:
#   Player 1 indiff: q*S - H = 0  =>  q* = H/S
#   Player 2 indiff: p*S - H = 0  =>  p* = H/S
# ============================================================================

#' @keywords internal
.geio_stag_hunt <- function() {
  list(
    type    = "stag_hunt",
    n_x     = 2L,
    n_theta = 2L,

    F_theta = function(theta, x) {
      .geio_assert_theta(theta, 2L, "stag_hunt")
      .geio_assert_x(x, 2L, "stag_hunt")
      S <- theta[1]; H <- theta[2]
      p <- x[1]; q <- x[2]

      # A = B = [[S, 0], [H, H]]
      # Player 1 indiff: q*(S - H - 0 + H) + (0 - H) = q*S - H
      f1 <- q * S - H
      # Player 2 indiff (symmetric): p*S - H
      f2 <- p * S - H

      c(f1, f2)
    },

    DF_theta = function(theta, x) {
      .geio_assert_theta(theta, 2L, "stag_hunt")
      S <- theta[1]
      # df1/dp = 0, df1/dq = S
      # df2/dp = S, df2/dq = 0
      matrix(c(0, S, S, 0), nrow = 2, ncol = 2)
    },

    x0_default = function(theta) c(0.5, 0.5),

    all_pure = function(theta) {
      # (1,1) = both stag, (0,0) = both hare
      list(c(1, 1), c(0, 0))
    },

    # Pure NE for Stag Hunt: (1,1) and (0,0) when S > H > 0.
    # At (1,1): both get S; deviate to hare gets H < S. Strict, index +1.
    # At (0,0): both get H; deviate to stag gets 0 < H. Strict, index +1.
    is_valid = function(theta, x, tol = 1e-8) {
      length(x) == 2L && all(is.finite(x)) && all(x >= -tol) && all(x <= 1 + tol)
    },

    pure_equilibria = function(theta) {
      S <- theta[1]; H <- theta[2]
      eqs <- list()
      # (1,1): both stag. Each gets S, deviate gives H.
      if (S >= H)
        eqs <- c(eqs, list(list(
          x = c(1, 1), index = if (S > H) 1L else 0L)))
      # (0,0): both hare. Each gets H, deviate gives 0.
      if (H >= 0)
        eqs <- c(eqs, list(list(
          x = c(0, 0), index = if (H > 0) 1L else 0L)))
      eqs
    }
  )
}


# ============================================================================
# Battle of the Sexes
# ============================================================================
#
# Parameterization: theta = c(alpha, beta)
#   alpha = coordination payoff at player 1's preferred outcome
#   beta  = coordination payoff at player 2's preferred outcome
#   Convention: alpha > 0, beta > 0
#
# Payoff matrices:
#   A = [[alpha, 0], [0, beta]]
#   B = [[beta, 0],  [0, alpha]]
#
# Both go to opera:    P1 gets alpha, P2 gets beta
# Both go to football: P1 gets beta,  P2 gets alpha
# Mismatch:            both get 0
# ============================================================================

#' @keywords internal
.geio_bos <- function() {
  list(
    type    = "bos",
    n_x     = 2L,
    n_theta = 2L,

    F_theta = function(theta, x) {
      .geio_assert_theta(theta, 2L, "bos")
      .geio_assert_x(x, 2L, "bos")
      alpha <- theta[1]; beta <- theta[2]
      p <- x[1]; q <- x[2]

      # A = [[alpha, 0], [0, beta]]
      # a11=alpha, a12=0, a21=0, a22=beta
      # f1 = q*(alpha - 0 - 0 + beta) + (0 - beta) = q*(alpha+beta) - beta
      f1 <- q * (alpha + beta) - beta

      # B = [[beta, 0], [0, alpha]]
      # b11=beta, b12=0, b21=0, b22=alpha
      # f2 = p*(beta - 0 - 0 + alpha) + (0 - alpha) = p*(alpha+beta) - alpha
      f2 <- p * (alpha + beta) - alpha

      c(f1, f2)
    },

    DF_theta = function(theta, x) {
      .geio_assert_theta(theta, 2L, "bos")
      alpha <- theta[1]; beta <- theta[2]
      # df1/dp = 0, df1/dq = alpha + beta
      # df2/dp = alpha + beta, df2/dq = 0
      ab <- alpha + beta
      matrix(c(0, ab, ab, 0), nrow = 2, ncol = 2)
    },

    x0_default = function(theta) c(0.5, 0.5),

    all_pure = function(theta) {
      # (1,1) = both opera, (0,0) = both football
      list(c(1, 1), c(0, 0))
    },

    # Pure NE for BoS: (1,1) and (0,0) when alpha, beta > 0.
    # A = [[alpha,0],[0,beta]], B = [[beta,0],[0,alpha]]
    # (1,1): row gets alpha vs 0, col gets beta vs 0. Both strict, index +1.
    # (0,0): row gets beta vs 0, col gets alpha vs 0. Both strict, index +1.
    is_valid = function(theta, x, tol = 1e-8) {
      length(x) == 2L && all(is.finite(x)) && all(x >= -tol) && all(x <= 1 + tol)
    },

    pure_equilibria = function(theta) {
      alpha <- theta[1]; beta <- theta[2]
      eqs <- list()
      # (1,1): a11=alpha >= a21=0, b11=beta >= b12=0
      if (alpha >= 0 && beta >= 0)
        eqs <- c(eqs, list(list(
          x = c(1, 1), index = if (alpha > 0 && beta > 0) 1L else 0L)))
      # (0,0): a22=beta >= a12=0, b22=alpha >= b21=0
      if (beta >= 0 && alpha >= 0)
        eqs <- c(eqs, list(list(
          x = c(0, 0), index = if (beta > 0 && alpha > 0) 1L else 0L)))
      eqs
    }
  )
}


# ============================================================================
# Custom game
# ============================================================================

#' @keywords internal
.geio_custom <- function(F_theta, DF_theta, n_x, n_theta) {
  if (is.null(F_theta) || !is.function(F_theta))
    stop("geio_game('custom'): F_theta must be a function")
  if (!is.null(DF_theta) && !is.function(DF_theta))
    stop("geio_game('custom'): DF_theta must be NULL or a function")
  if (is.null(n_x) || !is.numeric(n_x) || length(n_x) != 1L || !is.finite(n_x) || n_x < 1)
    stop("geio_game('custom'): n_x must be a positive integer-like scalar")
  if (is.null(n_theta) || !is.numeric(n_theta) || length(n_theta) != 1L || !is.finite(n_theta) || n_theta < 1)
    stop("geio_game('custom'): n_theta must be a positive integer-like scalar")
  n_x <- as.integer(n_x)
  n_theta <- as.integer(n_theta)

  if (is.null(DF_theta))
    DF_theta <- function(theta, x) geio_numjac(function(xx) F_theta(theta, xx), x)

  list(
    type       = "custom",
    n_x        = as.integer(n_x),
    n_theta    = as.integer(n_theta),
    F_theta    = F_theta,
    DF_theta   = DF_theta,
    x0_default = function(theta) rep(0.5, n_x),
    all_pure   = NULL,
    is_valid   = function(theta, x, tol = 1e-8) length(x) == n_x && all(is.finite(x)),
    # Custom games: user must supply pure_equilibria if boundary
    # equilibria exist.  Without it, only interior equilibria from
    # Newton are reported and the index-sum check may be incomplete.
    pure_equilibria = NULL
  )
}


# ============================================================================
# Belloc et al. (2019) experimental Stag Hunt parameterizations
# ============================================================================
#
# Four Stag Hunt games from Belloc, Bilancini, Boncinelli & D'Alessandro
# (2019). "Intuition and Deliberation in the Stag Hunt Game."
# Scientific Reports 9:14833. doi:10.1038/s41598-019-50556-8
#
# Games differ by basin of attraction of the stag equilibrium:
#   Game 1: basin = 1/4, S=4, H=3 (A=[3,0;3,4])  -- baseline
#   Game 2: basin = 1/4, bimatrix (A=[4,1;4,5])    -- shifted +1
#   Game 3: basin = 3/8, S=4, H=2.5                -- widest basin
#   Game 4: basin = 1/8, S=4, H=3.5                -- narrowest basin
#
# Experimental result: 62.85% play stag under time pressure vs 52.32%
# under deliberation.  Basin of attraction positively predicts stag.
# ============================================================================

#' Belloc et al. (2019) experimental Stag Hunt games
#'
#' Returns the four parameterized Stag Hunt games from the Belloc et al.
#' (2019) laboratory experiment.  Each game is returned as a
#' \code{geio_game} object with the appropriate \code{theta}.  Games
#' are ordered by basin of attraction of the stag equilibrium (widest
#' to narrowest).
#'
#' @return A named list with elements \code{game3}, \code{game1},
#'   \code{game2}, \code{game4}, each a list containing:
#'   \describe{
#'     \item{\code{game}}{A \code{geio_game} object.}
#'     \item{\code{theta}}{Numeric parameter vector for this game.}
#'     \item{\code{basin}}{Basin of attraction of the stag equilibrium.}
#'     \item{\code{p_star}}{Mixed-strategy equilibrium probability.}
#'     \item{\code{label}}{Human-readable label.}
#'   }
#'
#'   Games 1, 3, and 4 use the symmetric \code{"stag_hunt"} template
#'   (stag-alone payoff 0).  Game 2 adds one point to every outcome, so
#'   its stag-alone payoff is nonzero; it is therefore built with the
#'   \code{"bimatrix"} template.  As a result, Game 2's equilibrium
#'   enumeration follows the bimatrix path rather than the stag-hunt
#'   pure-equilibrium path, though all four games share the same stag
#'   basin structure described above.
#'
#' @examples
#' bg <- geio_belloc_games()
#' # Game 3 has widest basin (3/8) => easiest to sustain stag
#' bg$game3$basin
#' bg$game3$game$F_theta(bg$game3$theta, c(bg$game3$p_star, bg$game3$p_star))
#'
#' @references
#' Belloc, M., Bilancini, E., Boncinelli, L., & D'Alessandro, S. (2019).
#' Intuition and Deliberation in the Stag Hunt Game.
#' \emph{Scientific Reports}, 9, 14833.
#' \doi{10.1038/s41598-019-50556-8}
#'
#' @export
geio_belloc_games <- function() {

  # Game 1: A(hare)=3, A(stag,stag)=4, A(stag,hare)=0
  # Symmetric: A = [[4,0],[3,3]], theta = (S=4, H=3)
  # p* = H/S = 3/4, basin = 1 - 3/4 = 1/4
  g1 <- list(
    game   = geio_game("stag_hunt"),
    theta  = c(S = 4, H = 3),
    basin  = 1/4,
    p_star = 3/4,
    label  = "Belloc Game 1 (baseline, basin=1/4)"
  )

  # Game 2: A(hare)=4, A(stag,stag)=5, A(stag,hare)=1
  # NOT symmetric in our stag_hunt template (stag-alone != 0)
  # Use bimatrix: A=[[5,1],[4,4]], B=A^T in row/column payoff convention
  # a11=5,a12=1,a21=4,a22=4, b11=5,b12=4,b21=1,b22=4
  # Player 1 indiff: q*(5-4-1+4)+(1-4) = 4q-3 = 0 => q=3/4
  # Player 2 indiff: p*(5-4-1+4)+(1-4) = 4p-3 = 0 => p=3/4
  g2 <- list(
    game   = geio_game("bimatrix"),
    theta  = c(5, 1, 4, 4, 5, 4, 1, 4),
    basin  = 1/4,
    p_star = 3/4,
    label  = "Belloc Game 2 (shifted +1, basin=1/4)"
  )

  # Game 3: A(hare)=2.5, A(stag,stag)=4, A(stag,hare)=0
  # theta = (S=4, H=2.5), p* = 2.5/4 = 5/8, basin = 3/8
  g3 <- list(
    game   = geio_game("stag_hunt"),
    theta  = c(S = 4, H = 2.5),
    basin  = 3/8,
    p_star = 5/8,
    label  = "Belloc Game 3 (widest basin=3/8)"
  )

  # Game 4: A(hare)=3.5, A(stag,stag)=4, A(stag,hare)=0
  # theta = (S=4, H=3.5), p* = 3.5/4 = 7/8, basin = 1/8
  g4 <- list(
    game   = geio_game("stag_hunt"),
    theta  = c(S = 4, H = 3.5),
    basin  = 1/8,
    p_star = 7/8,
    label  = "Belloc Game 4 (narrowest basin=1/8)"
  )

  list(game3 = g3, game1 = g1, game2 = g2, game4 = g4)
}
