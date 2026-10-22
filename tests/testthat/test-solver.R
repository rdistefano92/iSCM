test_that("solve_iscm recovers known effects to machine precision", {
  set.seed(1)
  m <- 4L
  Omega <- diag(1, m)
  for (i in seq_len(m)) Omega[i, -i] <- -runif(m - 1L, 0, 0.25)
  theta_true <- matrix(rnorm(m * 5), nrow = m)
  beta <- Omega %*% theta_true
  theta_hat <- solve_iscm(beta, Omega)
  expect_equal(theta_hat, theta_true, tolerance = 1e-10)
})

test_that("solve_iscm accepts a vector and returns a vector", {
  Omega <- matrix(c(1, -0.42, -0.33, 1), 2, 2, byrow = TRUE)
  beta <- c(1.5, -0.2)
  out <- solve_iscm(beta, Omega)
  expect_null(dim(out))
  expect_length(out, 2L)
  # closed form (eq. 4.1)
  den <- 1 - 0.42 * 0.33
  expect_equal(out[1], (beta[1] + 0.42 * beta[2]) / den, tolerance = 1e-12)
})

test_that("solve_iscm errors on a numerically singular Omega", {
  Omega <- matrix(c(1, -1, -1, 1), 2, 2)  # det = 0
  expect_error(solve_iscm(c(1, 1), Omega), "singular")
})

test_that("build_omega places cross-weights correctly", {
  cross <- list(c(u2 = 0.4, u3 = 0.1),
                c(u1 = 0.3),
                c(u1 = 0.2, u2 = 0.05))
  Om <- build_omega(cross, c("u1", "u2", "u3"))
  expect_equal(diag(Om), c(u1 = 1, u2 = 1, u3 = 1))
  expect_equal(Om["u1", "u2"], -0.4)
  expect_equal(Om["u2", "u1"], -0.3)
  expect_equal(Om["u3", "u2"], -0.05)
  expect_equal(Om["u1", "u3"], -0.1)
})
