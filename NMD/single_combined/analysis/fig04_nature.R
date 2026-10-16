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

vector_out <- TRUE      # also writes the PDF and SVG of the same figure

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
grp_labels <- c(outside_gene = "Outside a gene", no_ptc = "No PTC created",
                ptc = "PTC expected to trigger NMD")
d$group   <- factor(d$group, levels = grp_levels, labels = grp_labels[grp_levels])
d$library_label <- factor(d$library_label, levels = c("5ss_noP", "CPA_noP"))

# Fill is a pale tint of the outline colour: at 6 pt on a printed page a saturated fill fights the
# median line, and the tint still separates the three groups in greyscale.
# Keyed off grp_labels, not off literal strings: the third label changed once already, and a
# hardcoded key silently becomes NA - which ggplot draws as a colourless box rather than an error.
stroke <- setNames(c("#6f6f6a", "#2f6b2a", "#3a2d85"), grp_labels[grp_levels])
fill   <- setNames(c("#d8d8d5", "#b9dcb3", "#c3bce2"), grp_labels[grp_levels])

# Whiskers are the 5th and 95th percentiles, not 1.5 IQR: a few barcodes sit near log2FC 12 and the
# full range would squash every box. geom_boxplot(stat = "identity") draws exactly these five
# numbers, so nothing is recomputed behind our back and no outlier layer is added.
q <- function(v, p) as.numeric(stats::quantile(v, p, na.rm = TRUE, type = 7))
stats <- do.call(rbind, lapply(split(d, list(d$library_label, d$group), drop = TRUE), function(s) {
  data.frame(library_label = s$library_label[1], group = s$group[1],
             ymin = q(s$lfc, 0.05), lower = q(s$lfc, 0.25), middle = q(s$lfc, 0.50),
             upper = q(s$lfc, 0.75), ymax = q(s$lfc, 0.95), n = nrow(s))
}))

# Wilcoxon rank-sum tests: group 1 against group 2, and group 1 against group 3, where the groups
# are numbered AS THE FIGURE READS - top to bottom. So group 1 is the top row, "PTC expected to
# trigger NMD", and it is the baseline both comparisons are made against.
#
# Note this is the reverse of grp_levels, which is written bottom-to-top because that is the order
# ggplot lays a discrete y axis out in. Reading the numbering off the factor instead of off the
# printed figure put the baseline on "Outside a gene" once already.
#
# Two-sided and unpaired (wilcox.test's defaults), computed per library on the raw per-barcode
# log2FC, not on the five summary numbers the boxes are drawn from.
#
# With n in the thousands R uses the normal approximation with a continuity correction rather than
# the exact distribution. That is the intended behaviour and the reason for suppressWarnings, which
# silences only the "cannot compute exact p-value with ties" note.
n_grp <- length(grp_levels)
ref_g <- n_grp                       # top row of the figure = last factor level
other_g <- c(n_grp - 1L, n_grp - 2L) # the 2nd and 3rd rows from the top
tests <- do.call(rbind, lapply(levels(d$library_label), function(lib) {
  do.call(rbind, lapply(other_g, function(g) {
    a <- d$lfc[d$library_label == lib & d$group == levels(d$group)[ref_g]]
    b <- d$lfc[d$library_label == lib & d$group == levels(d$group)[g]]
    w <- suppressWarnings(stats::wilcox.test(a, b))
    # lo_g / hi_g are the bracket's ends in plotting coordinates, so they are sorted; ref_g and g
    # say which side is the baseline, which is what the console line reports.
    data.frame(library_label = lib, ref_g = ref_g, other_g = g,
               lo_g = min(ref_g, g), hi_g = max(ref_g, g), p = w$p.value)
  }))
}))
tests$library_label <- factor(tests$library_label, levels = levels(d$library_label))

# Holm-Bonferroni WITHIN each library. The family is the set of comparisons against the shared
# baseline of one experiment; ave() splits on library_label so 5ss_noP and CPA_noP are corrected
# independently of each other. They are different constructs on different barcodes - separate
# experiments, not two looks at one question - so a family spanning both would make each penalise
# the other. Holm rather than plain Bonferroni: same guarantee on the family-wise error rate,
# uniformly more powerful, no reason to prefer the weaker one.
tests$p_adj <- ave(tests$p, tests$library_label,
                   FUN = function(v) p.adjust(v, method = "holm"))

# The number itself, not a star. Stars were tried and dropped: the ladder saturates at four, so
# 5.8e-36 and 2.9e-56 print identically, and every comparison on this panel is far past that
# ceiling - the whole figure read "**** **** **** ****" and said nothing about how far past.
#
# The ADJUSTED p is what is drawn, so the panel shows what survives the correction. The raw value
# goes to the console beside it.
#
# Written as plotmath, not as unicode superscripts. Unicode was tried first and is wrong here: the
# base-14 PDF device cannot encode U+00D7 or U+207B, and pdf() replaces each byte with a dot, so
# the PNG read "P = 5.8 x 10^-36" while the PDF - the file that goes to the journal - read
# "P = 5.8 . 10.........". 48 conversion warnings, one per character, and a silently broken figure.
# plotmath renders on every device, and it gives an italic P for free, which is the convention.
#
# Each label is parsed on its own, so exponents of different widths across the facets are fine.
fmt_p <- function(p) {
  if (!is.finite(p)) return("")
  if (p >= 0.01)  return(sprintf("italic(P) == %.2f", p))
  if (p >= 0.001) return(sprintf("italic(P) == %.3f", p))
  e <- floor(log10(p))
  m <- p / 10^e
  sprintf("italic(P) == %.1f %%*%% 10^%d", m, e)
}
tests$label <- vapply(tests$p_adj, fmt_p, character(1))
cat("  Wilcoxon rank-sum, two-sided; Holm-adjusted within each library - for the caption:\n")
for (i in seq_len(nrow(tests))) {
  cat(sprintf("    %-8s %-28s vs %-28s p = %.3g, adjusted %.3g   %s\n",
              as.character(tests$library_label[i]),
              levels(d$group)[tests$ref_g[i]], levels(d$group)[tests$other_g[i]],
              tests$p[i], tests$p_adj[i], tests$label[i]))
}

# Axis range from the whiskers themselves, so no tail is clipped, plus two strips on the right: the
# sample sizes as a column, then the comparison brackets with their p-values. Breaks stop at the
# data, so the added strips do not grow a run of empty ticks.
lo <- min(stats$ymin); hi <- max(stats$ymax)
span <- hi - lo
pad <- 0.04 * span
pt <- 0.3527777778   # mm per pt: ggplot linewidth is mm, font size is mm of cap height
hair <- 0.25 * pt    # 0.25 pt hairline
# Two strips, both INSIDE the panel. The sample sizes go left, where the boxes leave room, because
# each library has its own and a single column outside the first panel would lose the second one's.
# The stars go right; being four characters at most they need a fraction of the room "P = 5.8 x
# 10^-42" did, which is why the panel no longer runs far past its own scale.
n_x <- lo - 0.24 * span
br_x <- c(hi + 0.05 * span, hi + 0.10 * span)     # one x per comparison, the wider one further out
br_tick <- 0.02 * span                            # short arms reaching back towards the boxes
# The stars are set rotated, so the right pad only has to clear the height of a line of 6 pt type
# rather than the width of four characters.
# The right pad now has to hold a whole "P = 5.8 x 10^-36" set horizontally at 6 pt. That is only
# affordable because the x rule is drawn as its own annotate() segment clipped to the breaks - with
# the theme's axis.line.x the rule would follow the strip out into the annotation.
xlim <- c(lo - 0.30 * span, hi + 0.50 * span)

# The taller bracket sits further out, whichever direction it points.
tests$x <- br_x[tests$hi_g - tests$lo_g]
# Midpoint of the bracket.
tests$y <- (tests$lo_g + tests$hi_g) / 2

# Both labels horizontal, one size, left-aligned at a single x clear of the outer bracket. The
# per-row angle and size that the star notation needed are gone with it: "P = 5.8 x 10^-36" has a
# reading direction whatever the comparison, and at 6 pt it is about 17 mm long - too long to set
# turned against a bracket that spans one row, which is only 10 mm of panel height.
#
# They sit half a row apart vertically, roughly 5 mm against a 2 mm line, so they cannot meet.
tests$lab_angle <- 0
tests$lab_size <- 4.6 * pt
tests$lab_x <- max(br_x) + 0.035 * span

# The sample sizes, one row per library and group, placed inside the panel that owns them.
sizes <- stats[, c("library_label", "group", "n")]

# The breaks, computed before the plot so the x rule can be clipped to them.
brk <- seq(2 * ceiling((lo - pad) / 2), 2 * floor((hi + pad) / 2), 2)

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
  geom_text(data = sizes, inherit.aes = FALSE,
            aes(x = n_x, y = group,
                label = paste0("n = ", format(n, big.mark = ",", trim = TRUE))),
            hjust = 0, size = 5.5 * pt, colour = "grey45") +
  # The x axis drawn as data, not as axis.line.x: the theme's rule spans the whole panel and would
  # run out into the two annotation strips. A scale has to stop where its numbers stop.
  annotate("segment", x = min(brk), xend = max(brk), y = -Inf, yend = -Inf,
           colour = "black", linewidth = hair, lineend = "square") +
  # The two comparison brackets: a vertical spine at its own x, two short arms reaching back
  # towards the rows it joins, and the p-value beside it. Drawn from `tests`, so the bracket and
  # the number can never disagree about which pair they refer to.
  geom_segment(data = tests, inherit.aes = FALSE,
               aes(x = x, xend = x, y = lo_g, yend = hi_g),
               colour = "grey35", linewidth = hair) +
  geom_segment(data = tests, inherit.aes = FALSE,
               aes(x = x, xend = x - br_tick, y = lo_g, yend = lo_g),
               colour = "grey35", linewidth = hair) +
  geom_segment(data = tests, inherit.aes = FALSE,
               aes(x = x, xend = x - br_tick, y = hi_g, yend = hi_g),
               colour = "grey35", linewidth = hair) +
  # Rotated to read along its own bracket. It also costs almost no horizontal room, which is why
  # the panel no longer needs a wide strip on the right. Angle, size and x come from the row (see
  # where `tests` gets them), so the two labels can be set differently in one layer.
  geom_text(data = tests, inherit.aes = FALSE,
            aes(x = lab_x, y = y, label = label, angle = lab_angle, size = lab_size),
            parse = TRUE, hjust = 0, vjust = 0.5, colour = "grey35") +
  scale_fill_manual(values = fill) +
  scale_colour_manual(values = stroke) +
  # `size` is already in the mm ggplot wants, so it is passed through rather than scaled.
  scale_size_identity() +
  scale_x_continuous(breaks = brk, expand = c(0, 0)) +
  coord_cartesian(xlim = xlim, clip = "off") +
  facet_wrap(~ library_label, nrow = 1) +
  labs(x = expression(log[2] * " fold change (CHX/EtOH)"), y = NULL) +
  theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 6, hjust = 0, margin = margin(b = 1.5, unit = "mm")),
    axis.line.x = element_blank(),      # drawn as an annotate() segment, clipped to the breaks
    axis.line.y = element_blank(),
    axis.ticks.x = element_line(colour = "black", linewidth = hair),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(1, "mm"),
    axis.text = element_text(size = 6, colour = "black"),
    axis.title.x = element_text(size = 6, colour = "black", margin = margin(t = 1, unit = "mm")),
    panel.spacing = unit(4, "mm"),
    plot.margin = margin(1.5, 2, 1, 1, unit = "mm")
  )

# 183 mm, the full double-column measure. It used to be 140 mm on the argument that three boxes and
# three short labels do not fill 183 mm and stretching them there gives the empty-panel look Nature
# figures avoid. The comparison brackets and their p-values added a strip of real content on the
# right of each panel, so the width is now earned rather than padding.
w_mm <- 183; h_mm <- 40
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
