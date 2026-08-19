#!/usr/bin/env Rscript
# Manuscript version of the PTC panel (figure 1G): the same three groups as fig04 of
# analysis/make_figures.py, drawn to the conventions Nature-family journals ask for.
#
# What "Nature style" means here, concretely:
#   * one sans-serif face at one size -- 6 pt for every label, sentence case, no bold anywhere
#     except a panel letter (which is set by the person assembling the figure, not here);
#   * no grid, no frame, no background tint: a single hairline x axis with outward ticks, and
#     nothing else that is not data;
#   * hairlines at 0.25 pt, the weight a printer can still hold at final size;
#   * the figure laid out at its printed width (183 mm, the double-column measure), so 6 pt is
#     6 pt on the page and nothing is rescaled afterwards;
#   * sample sizes on the figure itself rather than only in the caption.
#
# Why R rather than the matplotlib script that draws every other figure: the cluster has no Arial
# and no Helvetica, and matplotlib can only embed a font it can find, so its vector output is set
# in DejaVu Sans. R's pdf() device uses the base-14 PDF fonts, so asking it for Helvetica gives a
# PDF that genuinely references Helvetica and opens as Helvetica in Illustrator.
#
#   module load matplotlib/3.7.0-gfbf-2022b && python3 analysis/export_panel_data.py
#   module load R/4.3.1-foss-2022b && Rscript analysis/fig04_nature.R
#
# in   analysis/data/fig04_ptc_yesno.tsv
# out  analysis/pics/fig1G_ptc_yesno.png             600 dpi
#
# The reviewer asked for vector art, and the vector path is written and works -- set vector_out to
# TRUE below and the script also writes the PDF (real Helvetica, see the note at the device) and an
# SVG. It is off because the PNG is what is being circulated at the moment.

vector_out <- FALSE     # TRUE also writes the PDF and SVG of the same figure

suppressPackageStartupMessages({
  library(ggplot2)
  if (vector_out) library(svglite)
})

here  <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
if (is.na(here) || here == "") here <- "analysis"
pics  <- file.path(here, "pics")
src   <- file.path(here, "data", "fig04_ptc_yesno.tsv")
if (!file.exists(src)) stop("run analysis/export_panel_data.py first -- no ", src)

d <- read.delim(src, stringsAsFactors = FALSE)

# Bottom-to-top order of the rows, and left-to-right order of the panels.
grp_levels <- c("outside_gene", "no_ptc", "ptc")
grp_labels <- c(outside_gene = "Outside a gene", no_ptc = "No PTC created", ptc = "PTC created")
d$group   <- factor(d$group, levels = grp_levels, labels = grp_labels[grp_levels])
d$library_label <- factor(d$library_label, levels = c("Regular", "PCPA"))

# Fill is a pale tint of the outline colour: at 6 pt on a printed page a saturated fill fights the
# median line, and the tint still separates the three groups in greyscale.
stroke <- c("Outside a gene" = "#6f6f6a", "No PTC created" = "#2f6b2a", "PTC created" = "#3a2d85")
fill   <- c("Outside a gene" = "#d8d8d5", "No PTC created" = "#b9dcb3", "PTC created" = "#c3bce2")

# Whiskers are the 5th and 95th percentiles, not 1.5 IQR: a few barcodes sit near log2FC 12 and the
# full range would squash every box. geom_boxplot(stat = "identity") draws exactly these five
# numbers, so nothing is recomputed behind our back and no outlier layer is added.
q <- function(v, p) as.numeric(stats::quantile(v, p, na.rm = TRUE, type = 7))
stats <- do.call(rbind, lapply(split(d, list(d$library_label, d$group), drop = TRUE), function(s) {
  data.frame(library_label = s$library_label[1], group = s$group[1],
             ymin = q(s$lfc, 0.05), lower = q(s$lfc, 0.25), middle = q(s$lfc, 0.50),
             upper = q(s$lfc, 0.75), ymax = q(s$lfc, 0.95), n = nrow(s))
}))

# Axis range from the whiskers themselves, so no tail is clipped, plus a narrow strip on the right
# where the sample sizes line up as a column. Breaks stop at the data, so the added strip does not
# grow a run of empty ticks.
lo <- min(stats$ymin); hi <- max(stats$ymax)
pad <- 0.04 * (hi - lo)
xlim <- c(lo - pad, hi + pad + 0.13 * (hi - lo))
n_x <- xlim[2]

pt <- 0.3527777778   # mm per pt: ggplot linewidth is mm, font size is mm of cap height
hair <- 0.25 * pt    # 0.25 pt hairline

p <- ggplot(stats, aes(y = group)) +
  # No change at zero is the comparison the whole panel is about, so the reference line is drawn
  # at about three times the hairline weight and dark enough to read against the boxes.
  geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.75 * pt, linetype = "dotted") +
  geom_boxplot(aes(xmin = ymin, xlower = lower, xmiddle = middle, xupper = upper, xmax = ymax,
                   fill = group, colour = group),
               stat = "identity", width = 0.55, linewidth = hair) +
  # The median is the number the reader takes away, so it is the one line drawn in black.
  geom_segment(aes(x = middle, xend = middle,
                   y = as.numeric(group) - 0.275, yend = as.numeric(group) + 0.275),
               colour = "black", linewidth = 0.6 * pt) +
  geom_text(aes(x = n_x, label = paste0("n = ", format(n, big.mark = ",", trim = TRUE))),
            hjust = 1, vjust = -1.15, size = 6 * pt, colour = "grey35") +
  scale_fill_manual(values = fill) +
  scale_colour_manual(values = stroke) +
  scale_x_continuous(breaks = seq(2 * ceiling((lo - pad) / 2), 2 * floor((hi + pad) / 2), 2),
                     expand = c(0, 0)) +
  coord_cartesian(xlim = xlim, clip = "off") +
  facet_wrap(~ library_label, nrow = 1) +
  labs(x = expression(log[2] * " fold change (CHX/EtOH)"), y = NULL) +
  theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 6, hjust = 0, margin = margin(b = 1.5, unit = "mm")),
    axis.line.x = element_line(colour = "black", linewidth = hair, lineend = "square"),
    axis.line.y = element_blank(),
    axis.ticks.x = element_line(colour = "black", linewidth = hair),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(1, "mm"),
    axis.text = element_text(size = 6, colour = "black"),
    axis.title.x = element_text(size = 6, colour = "black", margin = margin(t = 1, unit = "mm")),
    panel.spacing = unit(4, "mm"),
    plot.margin = margin(1.5, 2, 1, 1, unit = "mm")
  )

# 140 mm rather than the full 183 mm double-column measure: three boxes and three short labels do
# not fill 183 mm, and stretching them there is the empty-panel look Nature figures avoid.
w_mm <- 140; h_mm <- 40
w_in <- w_mm / 25.4; h_in <- h_mm / 25.4
dir.create(pics, showWarnings = FALSE, recursive = TRUE)
out <- file.path(pics, "fig1G_ptc_yesno")

ragg::agg_png(paste0(out, ".png"), width = w_in, height = h_in, units = "in", res = 600,
              background = "white")
print(p); invisible(grDevices::dev.off())
written <- "png"

if (vector_out) {
  # pdf() and not cairo_pdf(): cairo would embed whatever sans the machine has (DejaVu here), while
  # the base device references Helvetica, which is what the journal expects and what opens as real
  # Helvetica in Illustrator.
  grDevices::pdf(paste0(out, ".pdf"), width = w_in, height = h_in,
                 family = "Helvetica", useDingbats = FALSE, bg = "white")
  print(p); invisible(grDevices::dev.off())

  svglite::svglite(paste0(out, ".svg"), width = w_in, height = h_in, bg = "white")
  print(p); invisible(grDevices::dev.off())
  # svglite names the font it actually found on this machine, and this machine only has DejaVu
  # Sans, so the SVG would ask every reader for a font the figure is not meant to be set in.
  # Rewrite the family to the intended stack. Only the declaration changes; the geometry stays as
  # laid out (with DejaVu metrics, which is why the PDF, not the SVG, is the file to judge the
  # spacing on).
  svg_path <- paste0(out, ".svg")
  writeLines(gsub("font-family: [^;]+;", "font-family: Helvetica, Arial, sans-serif;",
                  readLines(svg_path, warn = FALSE)), svg_path)
  written <- "{png,pdf,svg}"
}

cat(sprintf("  pics/fig1G_ptc_yesno.%s   %g x %g mm\n", written, w_mm, h_mm))
