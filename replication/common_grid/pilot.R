args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, dir.exists(args[1]))
out <- args[1]
writeLines(capture.output(sessionInfo()), file.path(out, "r_session_info.txt"))

alpha <- 0.05
rho <- 1
envelope <- 1 / 3
t0 <- 1 / 2
seed <- 20261124L
replicates <- 1000L
sample_sizes <- c(2048L, 32768L)
grids <- list(
  fine = (1:64) / 128,
  coarse = c(1 / 8, 1 / 4, 3 / 8, 1 / 2),
  gapped = c(1 / 128, 1 / 4, 3 / 8, 1 / 2)
)
laws <- list(
  linear = list(diagonal = function(p) p / 2, lambda = 1 / 2),
  clayton_1 = list(diagonal = function(p) p / (2 - p), lambda = 1 / 2),
  zero_tail = list(diagonal = function(p) p^2 / 3, lambda = 0)
)
method_names <- c("aligned", "weighted_intersection", "bracket", "latent_reference")
scenario_count <- length(sample_sizes) * length(grids) * length(laws)
summary_family_size <- scenario_count * length(method_names)

clip <- function(x) pmin(1, pmax(0, x))

binomial_limits <- function(count, n, error) {
  error <- rep_len(error, length(count))
  lower <- numeric(length(count))
  upper <- rep(1, length(count))
  nonzero <- count > 0
  below_n <- count < n
  lower[nonzero] <- qbeta(error[nonzero] / 2, count[nonzero], n - count[nonzero] + 1)
  upper[below_n] <- qbeta(1 - error[below_n] / 2, count[below_n] + 1, n - count[below_n])
  list(lower = lower, upper = upper)
}

expand <- function(lower, upper, p) {
  list(
    lower = clip(clip(lower / p) - envelope * p^rho),
    upper = clip(clip(upper / p) + envelope * p^rho)
  )
}

fixed_threshold <- function(count, n, p, error = alpha) {
  limits <- binomial_limits(count, n, error)
  expand(limits$lower, limits$upper, p)
}

max_count_matrix <- function(cutoffs, n, law) {
  atom_probs <- diff(c(0, law$diagonal(cutoffs), 1))
  stopifnot(all(atom_probs >= 0), abs(sum(atom_probs) - 1) < 1e-12)
  atoms <- rmultinom(replicates, n, atom_probs)
  t(vapply(seq_len(replicates), function(index) cumsum(atoms[-nrow(atoms), index]),
           numeric(length(cutoffs))))
}

exact_metrics <- function(n, p, law) {
  counts <- 0:n
  intervals <- fixed_threshold(counts, n, p)
  weights <- dbinom(counts, n, law$diagonal(p))
  c(
    mean_length = sum(weights * (intervals$upper - intervals$lower)),
    coverage = sum(weights * (intervals$lower <= law$lambda & law$lambda <= intervals$upper))
  )
}

cutoff_score <- function(p, n) envelope * p^rho + 1 / sqrt(n * p)

run_scenario <- function(law, grid, n, law_name, grid_name) {
  latent_threshold <- min(t0, (1 / (2 * rho * envelope * sqrt(n)))^(1 / (rho + 0.5)))
  cutoffs <- sort(unique(c(grid, latent_threshold)))
  sample_counts <- max_count_matrix(cutoffs, n, law)
  grid_counts <- sample_counts[, match(grid, cutoffs), drop = FALSE]
  cutoff_rank <- rank(cutoff_score(grid, n), ties.method = "first")
  chosen_index <- which(cutoff_rank == 1)
  chosen <- grid[chosen_index]
  chosen_interval <- fixed_threshold(grid_counts[, chosen_index], n, chosen)
  latent_interval <- fixed_threshold(
    sample_counts[, match(latent_threshold, cutoffs)], n, latent_threshold
  )

  allocation <- alpha * cutoff_rank^(-2) / sum(cutoff_rank^(-2))
  simultaneous <- fixed_threshold(as.vector(grid_counts), n,
                                  rep(grid, each = replicates),
                                  rep(allocation, each = replicates))
  lower_matrix <- matrix(simultaneous$lower, nrow = replicates)
  upper_matrix <- matrix(simultaneous$upper, nrow = replicates)
  lower <- apply(lower_matrix, 1, max)
  upper <- apply(upper_matrix, 1, min)
  empty <- lower > upper
  lower[empty] <- 0
  upper[empty] <- 1

  preceding <- max(c(0, grid[grid <= latent_threshold]))
  following <- min(grid[grid >= latent_threshold])
  preceding_counts <- if (preceding == 0) rep(0L, replicates) else {
    grid_counts[, match(preceding, grid)]
  }
  lower_bracket <- binomial_limits(preceding_counts, n, alpha)$lower
  upper_bracket <- binomial_limits(grid_counts[, match(following, grid)], n, alpha)$upper
  bracket_interval <- expand(lower_bracket, upper_bracket, latent_threshold)

  intervals <- list(chosen_interval, list(lower = lower, upper = upper),
                    bracket_interval, latent_interval)
  finite_sum <- rbind(
    exact_metrics(n, chosen, law), c(NA_real_, NA_real_),
    c(NA_real_, NA_real_), exact_metrics(n, latent_threshold, law)
  )
  rownames(finite_sum) <- method_names
  metadata <- data.frame(
    law = law_name, grid = grid_name, n = n, true_lambda = law$lambda,
    aligned_threshold = chosen, continuous_threshold = latent_threshold,
    preceding = preceding, following = following
  )
  summaries <- vector("list", length(method_names))
  trials <- vector("list", length(method_names))
  for (index in seq_along(intervals)) {
    interval <- intervals[[index]]
    width <- interval$upper - interval$lower
    covered <- interval$lower <= law$lambda & law$lambda <= interval$upper
    stopifnot(all(width >= 0 & width <= 1), all(is.finite(width)))
    rate <- mean(covered)
    simulation_ci <- binomial_limits(sum(covered), replicates, alpha / summary_family_size)
    width_se <- sd(width) / sqrt(replicates)
    width_delta <- qt(1 - alpha / (2 * summary_family_size), replicates - 1L) * width_se
    summaries[[index]] <- cbind(
      metadata, method = method_names[index], replicates = replicates,
      covered = sum(covered), coverage = rate, mcse_coverage = sqrt(rate * (1 - rate) / replicates),
      mc_coverage_lower = simulation_ci$lower, mc_coverage_upper = simulation_ci$upper,
      mean_length = mean(width), mcse_length = width_se,
      mc_length_lower = clip(mean(width) - width_delta),
      mc_length_upper = clip(mean(width) + width_delta),
      fallback_count = if (index == 2) sum(empty) else 0L,
      exact_coverage = finite_sum[index, 2L], exact_mean_length = finite_sum[index, 1L]
    )
    trials[[index]] <- cbind(
      metadata, method = method_names[index], replicate = seq_len(replicates),
      lower = interval$lower, upper = interval$upper, covered = covered, interval_length = width,
      fallback = if (index == 2) empty else FALSE
    )
  }
  difference <- trials[[1]]$interval_length - trials[[2]]$interval_length
  difference_se <- sd(difference) / sqrt(replicates)
  difference_delta <- qt(1 - alpha / (2 * scenario_count), replicates - 1L) * difference_se
  paired <- cbind(metadata, mean_length_difference = mean(difference), mcse = difference_se,
                  mc_lower = mean(difference) - difference_delta,
                  mc_upper = mean(difference) + difference_delta)
  list(summary = do.call(rbind, summaries), trials = do.call(rbind, trials), paired = paired)
}

lower_square_cdf <- function(u, v) {
  small <- pmin(u, v)
  large <- pmax(u, v, .Machine$double.xmin)
  (1 / 2) * (1.5 * small - small^2 / (2 * large))
}

baseline_cdf <- function(u, v) {
  u_low <- pmin(u, t0)
  v_low <- pmin(v, t0)
  u_high <- pmax(0, u - t0)
  v_high <- pmax(0, v - t0)
  deficit <- function(x) x - lower_square_cdf(x, t0)
  lower_square_cdf(u_low, v_low) + 2 * u_high * deficit(v_low) +
    2 * v_high * deficit(u_low) + u_high * v_high
}

cutoff_function <- function(x) {
  z <- clip(2 * x - 1)
  1 - 3 * z^2 + 2 * z^3
}

perturbation_cdf <- function(u, v, h) {
  ratio <- ifelse(u + v > 0, (u * v)^2 / pmax(u + v, 1e-100)^3, 0)
  ratio * cutoff_function(u / h) * cutoff_function(v / h)
}

cell_probabilities <- function(cdf, u, v) {
  cdf_grid <- outer(u, v, cdf)
  nx <- length(u)
  ny <- length(v)
  cdf_grid[-1L, -1L] - cdf_grid[-nx, -1L] -
    cdf_grid[-1L, -ny] + cdf_grid[-nx, -ny]
}

check_perturbation <- function(h, a) {
  u <- sort(unique(c(0, a / 2, a, h, 0.5, 1)))
  v <- sort(unique(c(0, h * (1:47) / 48, h, 0.5, 1)))
  eta <- min(1 / (2 * 134), 8 * envelope / 2) * h
  q0 <- cell_probabilities(baseline_cdf, u, v)
  q1 <- cell_probabilities(
    function(x, y) baseline_cdf(x, y) + eta * perturbation_cdf(x, y, h), u, v
  )
  change <- max(abs(q1 - q0))
  chi_squared <- sum((q1[q0 > 0] - q0[q0 > 0])^2 / q0[q0 > 0])
  bound <- 2 * 134^2 * eta^2 * a
  row_error <- max(abs(rowSums(q1) - diff(u)))
  col_error <- max(abs(colSums(q1) - diff(v)))
  diagonal <- sort(unique(c(h / 2, h, 1 / 2, h * (1:101) / 101, h / 2^14)))
  lambda1 <- 1 / 2 + eta / 8
  ratio <- (baseline_cdf(diagonal, diagonal) + eta * perturbation_cdf(diagonal, diagonal, h)) /
    diagonal
  envelope_excess <- max(abs(ratio - lambda1) - envelope * diagonal^rho)
  stopifnot(min(q1) >= -1e-12, row_error < 1e-12, col_error < 1e-12,
            chi_squared <= bound + 1e-12, envelope_excess <= 1e-12,
            abs(ratio[1L] - lambda1) < 1e-12)
  if (a == 0) stopifnot(change < 1e-12)
  data.frame(h = h, a = a, min_cell_mass = min(q1), max_cell_change = change,
             max_row_error = row_error, max_column_error = col_error, chi_squared = chi_squared,
             chi_squared_bound = bound, envelope_excess = envelope_excess,
             lambda_separation = eta / 8)
}

math_check <- do.call(rbind, lapply(c(1 / 16, 1 / 4), function(h) {
  rbind(check_perturbation(h, 0), check_perturbation(h, h / 8), check_perturbation(h, h / 2))
}))
write.csv(math_check, file.path(out, "deterministic_checks.csv"), row.names = FALSE)

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(seed)
results <- list()
index <- 0L
for (n in sample_sizes) {
  for (grid_name in names(grids)) {
    for (law_name in names(laws)) {
      index <- index + 1L
      results[[index]] <- run_scenario(laws[[law_name]], grids[[grid_name]], n, law_name, grid_name)
    }
  }
}
summary <- do.call(rbind, lapply(results, function(result) result$summary))
trials <- do.call(rbind, lapply(results, function(result) result$trials))
paired <- do.call(rbind, lapply(results, function(result) result$paired))
stopifnot(nrow(summary) == summary_family_size, nrow(trials) == summary_family_size * replicates,
          nrow(paired) == scenario_count,
          all(summary$exact_coverage[!is.na(summary$exact_coverage)] >= 1 - alpha - 1e-12))
write.csv(summary, file.path(out, "pilot_summary.csv"), row.names = FALSE)
write.csv(trials, file.path(out, "pilot_trials.csv"), row.names = FALSE)
write.csv(paired, file.path(out, "pilot_paired_length.csv"), row.names = FALSE)
cat(sprintf("Synthetic pilot. R seed %d. %d scenarios, %d replicates each. Six deterministic checks.\n",
            seed, scenario_count, replicates))
cat(sprintf("Saved %d interval outcomes per method, %d in total.\n",
            scenario_count * replicates, nrow(trials)))
