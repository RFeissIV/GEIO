## GEIO 0.1.0 — deduplicated, claim-disciplined test suite

test_that("package exports exist", {
  for (fn in c("GEIO","geio_game","geio_numjac","geio_newton","geio_index","geio_fragility","geio_belloc_games"))
    expect_true(is.function(get(fn)), label = fn)
})

# --- numjac ---
test_that("geio_numjac approximates analytical Jacobian", {
  J <- geio_numjac(function(x) c(x[1]^2+x[2], x[1]*x[2]-1), c(2,3))
  expect_lt(max(abs(J - matrix(c(4,3,1,2),2,2))), 1e-5)
})
test_that("geio_numjac rejects non-finite x", {
  expect_error(geio_numjac(identity, c(1, NA)), "finite")
})

# --- newton ---
test_that("geio_newton solves linear system", {
  sol <- geio_newton(function(x) c(2*x[1]+x[2]-5, x[1]+3*x[2]-7), c(0,0))
  expect_true(sol$converged)
})
test_that("geio_newton reports non-convergence", {
  sol <- geio_newton(function(x) c(x[1]^2+x[2]^2+1, x[1]+x[2]), c(0,0), max_iter=10L)
  expect_false(sol$converged)
})
test_that("geio_newton rejects bad damping", {
  expect_error(geio_newton(function(x) x^2, 0.5, damping=-1), "positive")
})

# --- index ---
test_that("geio_index signs correct", {
  expect_equal(geio_index(matrix(2)), 1L)
  expect_equal(geio_index(matrix(-2)), -1L)
  expect_equal(geio_index(matrix(0)), 0L)
})
test_that("geio_index non-finite J returns 0", {
  expect_equal(geio_index(matrix(NA_real_,1,1)), 0L)
})
test_that("geio_index rejects non-numeric J", {
  expect_error(geio_index(matrix("a",1,1)), "numeric")
})
test_that("geio_index rejects bad tol", {
  expect_error(geio_index(J=diag(2), tol=-1), "non-negative")
})

# --- fragility ---
test_that("geio_fragility correct values", {
  f <- geio_fragility(J=diag(c(5,0.01)))
  expect_equal(f$sigma_min, 0.01); expect_equal(f$sigma_max, 5)
})
test_that("geio_fragility rejects non-square J", {
  expect_error(geio_fragility(J=matrix(1:6,2,3)), "square")
})
test_that("geio_fragility rejects non-finite J", {
  expect_error(geio_fragility(J=matrix(NA_real_,1,1)), "finite")
})
test_that("geio_fragility rejects non-numeric J", {
  expect_error(geio_fragility(J=matrix("a",1,1)), "numeric")
})

# --- games ---
test_that("stag_hunt creates valid game", {
  game <- geio_game("stag_hunt")
  expect_s3_class(game, "geio_game")
  expect_lt(max(abs(game$F_theta(c(4,2), c(0.5,0.5)))), 1e-10)
})
test_that("bos creates valid game", {
  expect_lt(max(abs(geio_game("bos")$F_theta(c(3,1), c(0.75,0.25)))), 1e-10)
})
test_that("custom game works", {
  g <- geio_game("custom", F_theta=function(t,x) c(t[1]*x[1]-1), n_x=1L, n_theta=1L)
  expect_equal(g$type, "custom")
})
test_that("custom rejects bad dimensions", {
  expect_error(geio_game("custom", F_theta=identity, n_x=-1, n_theta=1), "positive")
})
test_that("built-in games reject wrong theta length", {
  expect_error(geio_game("stag_hunt")$F_theta(c(1), c(0.5,0.5)), "theta")
})
test_that("built-in games reject wrong x length", {
  expect_error(geio_game("stag_hunt")$F_theta(c(4,2), c(0.5)), "x")
})

# --- mathematical sanity ---
test_that("stag hunt mixed eq has negative local index", {
  expect_equal(geio_index(geio_game("stag_hunt")$DF_theta(c(4,2), c(0.5,0.5))), -1L)
})
test_that("bos mixed eq has negative local index", {
  expect_equal(geio_index(geio_game("bos")$DF_theta(c(3,1), c(0.75,0.25))), -1L)
})
test_that("stag hunt index sum +1 under S>H>0", {
  game <- geio_game("stag_hunt"); theta <- c(4,1)
  mixed <- geio_index(game$DF_theta(theta, c(0.25,0.25)))
  pure <- vapply(game$pure_equilibria(theta), function(z) z$index, integer(1))
  expect_equal(sum(c(mixed, pure)), 1L)
})
test_that("bos index sum +1 under alpha,beta>0", {
  game <- geio_game("bos"); theta <- c(3,2)
  mixed <- geio_index(game$DF_theta(theta, c(0.6,0.4)))
  pure <- vapply(game$pure_equilibria(theta), function(z) z$index, integer(1))
  expect_equal(sum(c(mixed, pure)), 1L)
})

# --- Belloc ---
test_that("Belloc games return 4 entries in documented order", {
  bg <- geio_belloc_games()
  expect_equal(length(bg), 4)
  expect_equal(names(bg), c("game3","game1","game2","game4"))
})
test_that("Belloc game1 mixed eq at p*=3/4", {
  g <- geio_belloc_games()$game1
  expect_lt(max(abs(g$game$F_theta(g$theta, c(g$p_star,g$p_star)))), 1e-10)
})
test_that("all Belloc games have finite diagnostics", {
  for (g in geio_belloc_games()) {
    J <- g$game$DF_theta(g$theta, c(g$p_star,g$p_star))
    expect_true(all(is.finite(J)))
    expect_true(is.finite(geio_fragility(J=J)$sigma_min))
  }
})

# --- GALAHAD integration ---
test_that("GEIO optimizer runs", {
  skip_on_cran(); skip_if_not_installed("GALAHAD")
  fit <- GEIO(geio_game("stag_hunt"), function(t) sum((t-c(3,1))^2),
              function(t) 2*(t-c(3,1)), c(4,2),
              list(positive=c(1L,2L), euclidean=integer(0)),
              control=list(max_iter=10))
  expect_s3_class(fit, "geio_result")
})
test_that("print.geio_result works", {
  skip_on_cran(); skip_if_not_installed("GALAHAD")
  fit <- GEIO(geio_game("stag_hunt"), function(t) sum((t-c(3,1))^2),
              function(t) 2*(t-c(3,1)), c(4,2),
              list(positive=c(1L,2L), euclidean=integer(0)),
              control=list(max_iter=5))
  expect_output(print(fit), "GEIO")
})
