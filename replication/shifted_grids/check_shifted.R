args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, dir.exists(args[1]))
out <- args[1]
writeLines(capture.output(sessionInfo()), file.path(out, "r_session_info.txt"))

t0 <- 1 / 2
tau <- 1 / 2
envelope <- 1 / 3
alpha <- 0.05
tolerance <- 1e-10
clip <- function(x) pmax(0, pmin(1, x))

lower_square_cdf <- function(u, v, gap) {
  small <- pmin(u, v)
  large <- pmax(u, v)
  end_cone <- pmin(large, 2 * gap, 2 * small)
  cone_width <- pmax(0, end_cone - small)
  cone <- 2 * small * log1p(cone_width / pmax(small, 1e-100)) - cone_width
  start_strip <- pmax(small, 2 * gap)
  end_strip <- pmin(large, small + gap)
  strip_width <- pmax(0, end_strip - start_strip)
  strip <- strip_width * (1 + (small - (start_strip + end_strip) / 2) / gap)
  tau * (small + (cone + strip) / 2)
}

baseline_cdf <- function(u, v, gap) {
  ul <- pmin(u, t0)
  vl <- pmin(v, t0)
  uh <- pmax(0, u - t0)
  vh <- pmax(0, v - t0)
  deficit <- function(a) a - lower_square_cdf(a, t0, gap)
  lower_square_cdf(ul, vl, gap) + 2 * uh * deficit(vl) +
    2 * vh * deficit(ul) + 2 * tau * uh * vh
}

rectangles <- function(cdf, u, v) {
  probabilities <- outer(u, v, cdf)
  m <- length(u)
  n <- length(v)
  probabilities[-1L, -1L] - probabilities[-m, -1L] -
    probabilities[-1L, -n] + probabilities[-m, -n]
}

lower_square_quadrature <- function(u, v, gap) {
  small <- min(u, v)
  large <- max(u, v)
  if (small <= 0 || large == small) return(tau * small)
  breakpoints <- sort(unique(c(small, large, 2 * gap, 2 * small, small + gap)))
  breakpoints <- breakpoints[breakpoints >= small & breakpoints <= large]
  tail_fraction <- function(s) {
    width <- pmin(s / 2, gap)
    clip((small - s + width) / width)
  }
  integrals <- vapply(seq_len(length(breakpoints) - 1L), function(i) {
    integrate(tail_fraction, breakpoints[i], breakpoints[i + 1L],
              rel.tol = 1e-10, abs.tol = 1e-13)$value
  }, numeric(1))
  tau * small + tau * sum(integrals) / 2
}

psi <- function(x) {
  scaled <- clip(2 * x - 1)
  1 - 3 * scaled^2 + 2 * scaled^3
}

phi <- function(x) pmax(0, 1 - x^2)^3

perturbation <- function(u, v, gap, h) {
  total <- u + v
  total / 2 * phi(2 * (u - v) / gap) * phi(4 * (u - v) / pmax(total, 1e-100)) *
    psi(u / h) * psi(v / h)
}

check_construction <- function(cells, fraction, rho, shift_name) {
  delta <- 1 / cells
  gap <- delta * min(fraction, 1 - fraction)
  h <- min(t0, gap^(1 / (rho + 1)))
  eta <- h^rho / 30000
  gu <- (0:cells) * delta
  gv <- sort(unique(c(0, 1, ((0:(cells - 1)) + fraction) * delta)))
  c0 <- function(u, v) baseline_cdf(u, v, gap)
  c1 <- function(u, v) c0(u, v) + eta * perturbation(u, v, gap, h)
  label0 <- rectangles(c0, gu, gv)
  label1 <- rectangles(c1, gu, gv)
  mass_change <- max(abs(label0 - label1))

  inside <- h * c(1 / 64, 1 / 8, 1 / 4, 3 / 8, 1 / 2, 3 / 4, 1)
  offsets <- pmin(t0, inside + gap / 3)
  ref_u <- sort(unique(c(gu, inside, offsets)))
  ref_v <- sort(unique(c(gv, inside, offsets)))
  refined0 <- rectangles(c0, ref_u, ref_v)
  refined1 <- rectangles(c1, ref_u, ref_v)
  row_error <- max(abs(rowSums(refined1) - diff(ref_u)))
  col_error <- max(abs(colSums(refined1) - diff(ref_v)))
  min_refined_mass <- min(refined0, refined1)
  max_refined_change <- max(abs(refined1 - refined0))

  p <- sort(unique(c(inside, t0, exp(seq(log(h / 100), log(t0), length.out = 63)))))
  limit <- tau + eta
  diagonal <- c1(p, p) / p
  max_envelope_excess <- max(abs(diagonal - limit) - envelope * p^rho)
  departure_in_core <- max(abs(diagonal[p <= h / 2] - limit))

  quad_u <- c(h / 8, gap / 2, 2 * gap / 3, h / 2, t0 / 2, t0 * 7 / 8)
  quad_v <- pmin(t0, quad_u + c(gap / 3, h / 3, gap, gap / 2, h / 4, gap / 5))
  quadrature <- mapply(lower_square_quadrature, quad_u, quad_v, MoreArgs = list(gap = gap))
  integration_error <- max(abs(quadrature - lower_square_cdf(quad_u, quad_v, gap)))

  stopifnot(mass_change <= tolerance, min_refined_mass >= -tolerance,
            max_refined_change > 0, row_error <= tolerance, col_error <= tolerance,
            max_envelope_excess <= tolerance, departure_in_core <= tolerance,
            integration_error <= tolerance, eta * 15000 * h / gap <= 1 / 2 + tolerance)
  data.frame(
    cells, shift_name, fraction, delta, gap, rho, h, lambda_separation = eta,
    observed_cells = length(label0), max_observed_mass_change = mass_change,
    refined_cells = length(refined1), max_refined_mass_change = max_refined_change,
    min_refined_mass, row_error, col_error, departure_in_core,
    max_envelope_excess, integration_error
  )
}

clopper_pearson <- function(s, n) {
  lower <- rep(0, length(s))
  upper <- rep(1, length(s))
  lower[s > 0] <- qbeta(alpha / 2, s[s > 0], n - s[s > 0] + 1)
  upper[s < n] <- qbeta(1 - alpha / 2, s[s < n] + 1, n - s[s < n])
  list(lower, upper)
}

check_interval <- function(cells, fraction, n, shift_name) {
  delta <- 1 / cells
  gap <- delta * min(fraction, 1 - fraction)
  p <- max(n^(-1 / 3), delta, sqrt(gap))
  # For rho = 1 the eligibility inequalities use the dyadic cutpoints
  # directly and preserve equality at an exact boundary.
  eligible <- function(t) t >= delta & t * t >= gap & n * t * t * t >= 1
  if (!eligible(1 / 4)) {
    return(data.frame(cells, shift_name, n, delta, gap, p, x = NA_real_,
                      y = NA_real_, q = NA_real_, coverage = 1, expected_length = 1))
  }
  candidate_x <- (1:(cells - 1)) * delta
  offset <- if (fraction <= 1 / 2) gap else -gap
  candidate_y <- candidate_x + offset
  first <- match(TRUE, eligible(pmin(candidate_x, candidate_y)))
  stopifnot(!is.na(first))
  x <- candidate_x[first]
  y <- candidate_y[first]
  t <- min(x, y)
  b <- min(abs(x - y), t)
  gu <- (0:cells) * delta
  gv <- sort(unique(c(0, 1, ((0:(cells - 1)) + fraction) * delta)))
  stopifnot(min(abs(x - gu)) <= tolerance, min(abs(y - gv)) <= tolerance,
            t >= p - tolerance, t <= p + delta + tolerance,
            abs(abs(x - y) - gap) <= tolerance)
  # C(x,y) and its diagonal limit for the Clayton copula with parameter one.
  q <- x * y / (x + y - x * y)
  true_lambda <- 1 / 2
  counts <- 0:n
  bounds <- clopper_pearson(counts, n)
  lower <- clip(pmin(t, bounds[[1]]) / t - b / t - envelope * t)
  upper <- clip(pmin(t, bounds[[2]]) / t + envelope * t)
  weights <- dbinom(counts, n, q)
  covered <- lower <= true_lambda & true_lambda <= upper
  coverage <- sum(weights[covered])
  expected_length <- sum(weights * (upper - lower))
  stopifnot(coverage >= 1 - alpha - tolerance, coverage <= 1 + tolerance,
            abs(sum(weights) - 1) <= tolerance, all(lower <= upper))
  data.frame(cells, shift_name, n, delta, gap, p, x, y, q, coverage, expected_length)
}

construction <- list()
for (cells in c(2L, 8L, 32L, 64L)) {
  shifts <- c(quarter = 1 / 4, half = 1 / 2, three_quarters = 3 / 4,
              near_one = 1 - 1 / cells, near_zero = 1 / cells^2)
  for (rho in c(1 / 2, 1, 2)) {
    for (shift in names(shifts)) {
      construction[[length(construction) + 1L]] <- check_construction(
        cells, shifts[[shift]], rho, shift
      )
    }
  }
}
interval <- list()
for (cells in c(8L, 32L, 64L)) {
  shifts <- c(aligned = 0, half = 1 / 2,
              near_one = 1 - 1 / cells, near_zero = 1 / cells^2)
  for (n in c(2048L, 32768L)) {
    for (shift in names(shifts)) {
      interval[[length(interval) + 1L]] <- check_interval(cells, shifts[[shift]], n, shift)
    }
  }
}
construction <- do.call(rbind, construction)
interval <- do.call(rbind, interval)
write.table(construction, file.path(out, "construction_checks.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(interval, file.path(out, "exact_binomial_checks.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
cat(sprintf("Deterministic: %d construction cases, %d binomial cases. No random seed.\n",
            nrow(construction), nrow(interval)))
cat(sprintf("Maximum observed cell change: %.12g. Maximum quadrature error: %.12g.\n",
            max(construction$max_observed_mass_change), max(construction$integration_error)))
cat(sprintf("Lowest exact binomial coverage in the checked synthetic model: %.12f.\n",
            min(interval$coverage)))
