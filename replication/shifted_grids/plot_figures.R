args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, file.exists(args[1]), dir.exists(args[2]))
exact <- read.delim(args[1], stringsAsFactors = FALSE)
out <- args[2]
writeLines(capture.output(sessionInfo()), file.path(out, "r_session_info.txt"))

geometry <- function(file) {
  cairo_pdf(file, width = 6.35, height = 3.55, pointsize = 10, family = "sans")
  on.exit(dev.off())
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, 0.14))
  delta <- 1 / 8
  threshold <- 2 * delta
  axis_positions <- (0:4) * delta
  axis_labels <- c("0", "0.125", "0.250", "0.375", "0.500")
  square_fill <- "#C7DFEB"
  strip_fill <- "#EFB989"
  grid <- "#B8BFC5"
  outline <- "#1D4D60"

  for (panel in 1:2) {
    gap <- c(0, delta / 2)[panel]
    label <- c("(a) Aligned cutpoints", "(b) Half-cell shift")[panel]
    par(mar = c(3.2, 4.05, 2, 0.65), mgp = c(1.9, 0.48, 0), tcl = -0.25,
        xaxs = "i", yaxs = "i", las = 1)
    plot.new()
    plot.window(xlim = c(0, 4 * delta), ylim = c(0, 4 * delta), asp = 1)
    rect(0, 0, threshold, threshold, col = square_fill, border = NA)
    if (gap > 0) {
      rect(0, threshold, threshold, threshold + gap, col = strip_fill, border = NA)
      rect(0, threshold, threshold, threshold + gap,
           density = 15, angle = 45, col = "#955D35", border = NA)
    }
    abline(v = axis_positions, h = axis_positions + gap, col = grid, lwd = 0.8)
    abline(a = 0, b = 1, col = "#52626B", lty = 3, lwd = 1.15)
    rect(0, 0, threshold, threshold + gap, border = outline, lwd = 2)
    if (gap > 0) {
      segments(0, threshold, threshold, threshold, col = outline, lty = 2, lwd = 1.2)
      arrows(threshold + 0.025, threshold, threshold + 0.025, threshold + gap,
             angle = 90, code = 3, length = 0.045, lwd = 1, col = "#8B421F")
      text(threshold + 0.053, threshold + gap / 2, "d", col = "#6E321A")
    }
    axis(1, at = axis_positions, labels = axis_labels)
    axis(2, at = axis_positions, labels = axis_labels)
    mtext("u", side = 1, line = 1.65)
    mtext("v", side = 2, line = 3.1, las = 0)
    title(main = label, cex.main = 0.98, font.main = 2)
  }

  par(mar = rep(0, 4))
  plot.new()
  legend("center", c("Diagonal event", "Additional strip"),
         fill = c(square_fill, strip_fill), border = NA, ncol = 2, bty = "n",
         cex = 0.94, x.intersp = 0.75)
}

lengths <- function(data, file) {
  sizes <- c(2048L, 32768L)
  cells <- c(8L, 32L, 64L)
  shifts <- c("aligned", "half", "near_one", "near_zero")
  stopifnot(nrow(data) == 24L, all(table(data$n, data$cells, data$shift_name) == 1L),
            setequal(data$n, sizes), setequal(data$cells, cells),
            setequal(data$shift_name, shifts), all(is.finite(data$expected_length)))
  palette <- c("#166083", "#A64819", "#317563", "#555F67")
  line_types <- c(1, 2, 3, 4)
  marks <- c(16, 17, 15, 1)
  ceiling <- max(0.6, max(data$expected_length))
  cairo_pdf(file, width = 6.65, height = 3.45, pointsize = 10, family = "sans")
  on.exit(dev.off())
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, 0.15))
  for (panel in seq_along(sizes)) {
    n <- sizes[panel]
    par(mar = c(3.6, 3.5, 2, 0.65), mgp = c(1.95, 0.48, 0), tcl = -0.28,
        xaxs = "i", yaxs = "i", las = 1)
    plot.new()
    plot.window(xlim = range(log2(cells)) + c(-0.14, 0.14), ylim = c(0, ceiling))
    y_ticks <- seq(0, 0.6, by = 0.1)
    abline(h = y_ticks, col = "#D6DCE0", lty = 3)
    for (j in seq_along(shifts)) {
      series <- data[data$n == n & data$shift_name == shifts[j], ]
      series <- series[order(series$cells), ]
      lines(log2(series$cells), series$expected_length, type = "o", col = palette[j],
            lty = line_types[j], lwd = 1.55, pch = marks[j], cex = 0.95 + 0.25 * (j == 4))
    }
    axis(1, at = log2(cells), labels = cells)
    axis(2, at = y_ticks, labels = formatC(y_ticks, digits = 1, format = "f"))
    box(bty = "l")
    mtext("Cells in first margin, M", side = 1, line = 1.95)
    if (panel == 1L) {
      mtext("Expected interval length", side = 2, line = 2.19, las = 0)
    }
    label <- paste0("(", letters[panel], ") n = ", format(n, big.mark = ","))
    title(main = label, cex.main = 1, font.main = 2)
  }
  par(mar = rep(0, 4))
  plot.new()
  keys <- c("Aligned", "Half-cell", "Near 1", "Near 0")
  legend("center", keys, col = palette, lty = line_types, pch = marks,
         lwd = 1.55, cex = 0.92, ncol = 4, bty = "n")
}

geometry(file.path(out, "Fig1.pdf"))
lengths(exact, file.path(out, "Fig2.pdf"))
cat("Two vector PDFs. The interval length panels contain all 24 fixed configurations.\n")
