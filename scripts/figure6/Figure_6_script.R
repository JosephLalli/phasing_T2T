# ============================================================================
# Figure 6 Script: Genomic Ideogram Plots with SVG Export
# ============================================================================
# This script creates ideogram plots showing variant density, switch error
# rates, and genotype error rates across chromosomes for both CHM13 and GRCh38
# reference genomes. It generates separate SVG files for each genome; use
# stitch_svgs.py to combine them into a single vector figure.
# ============================================================================

# Load required libraries
library(karyoploteR)
library(arrow)
library(PlotTools)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(svglite)

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]])))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }

  stop("Could not determine the Figure 6 script path.")
}

script_path <- get_script_path()
script_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(script_dir, "..", ".."))
setwd(project_dir)

rolling_data_candidates <- c(
  "intermediate_data_whole_genome/rolling_stats_500k_window.parquet",
  "intermediate_data/rolling_stats_500k_window.parquet"
)
rolling_data_matches <- rolling_data_candidates[file.exists(rolling_data_candidates)]

if (length(rolling_data_matches) == 0) {
  stop(
    "Could not find rolling_stats_500k_window.parquet in intermediate_data_whole_genome/ or intermediate_data/."
  )
}

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================

# File paths
CONFIG <- list(
  # Input data files
  rolling_data_file = rolling_data_matches[[1]],
  t2t_cytobands_file = 'resources/chm13v2.0_cytobands_allchrs.bed',
  grch38_cytobands_file = 'resources/grch38_cytobands_allchrs.bed',
  t2t_cytobands_header_file = "resources/chm13v2.0_cytobands_allchrs.w_header.bed",
  grch38_cytobands_header_file = "resources/grch38_cytobands_allchrs.w_header.bed",
  cnv_file = 'resources/decipher_syndromes.txt',
  t2t_segdup_file = "resources/CHM13_segdups_gt10kb.bed",
  grch38_segdup_file = "resources/GRCh38_segdups_gt10kb.bed",

  # Output directory
  output_dir = 'figures/figure6',

  # Plot parameters
  plot = list(
    # Error rate limits
    ser_ymax = 10,
    gt_ymax = 10,
    num_ticks = 6,
    
    # Colors
    var_density_palette = c("white", "#F2AF29", "#EB6424"),
    ser_color = "#EB6424",
    gt_color = '#EEB902',
    chm13_ser_color = "#593F62",
    chm13_gt_color = '#7B6D8D',
    grch38_ser_color = "#8499B1",
    grch38_gt_color = '#A5C4D4',
    cnv_highlight_color = "#FFDDDD",
    segdup_highlight_color = "#DDDDFF",
    
    # Font sizes (in points)
    default_pts = 6,
    linewidth_pts = 0.5,
    percent_size_pts = 5,
    chrom_size_pts = 8,
    label_size_pts = 6,
    base_number_size_pts = 5,
    regional_var_label_margin = 0.05,
    regional_data_label_margin = 0.09,
    figure6_var_label_margin = 0.02,
    figure6_data_label_margin = 0.05,
    
    # Region padding
    slop = 1.5e6,

    # Image dimensions (used by svglite)
    width_inches = 3.5,
    height_inches = 3.5,
    figure6_panel_width_inches = 3.0,
    figure6_panel_height_inches = 1.62,
    figure6_colorbar_width_inches = 3.4,
    figure6_colorbar_height_inches = 0.6,
    figure6_colorbar_axis_width_pt = 219.0
  ),
  
  # Genomic regions of interest (CHM13 coordinates)
  chm13_regions = list(
    chr22_22q11 = list(chr = "chr22", start = 19397653, end = 23803198, name = "chm13_22q11"),
    chr15_anglemans = list(chr = "chr15", start = 20807475, end = 25935855, name = "chm13_anglemans"),
    chr16_16p12.1 = list(chr = "chr16", start = 21398763, end = 30473116, name = "chm13_16p12.1"),
    chr1_1q21.1 = list(chr = "chr1", start = 145165377, end = 147746604, name = "chm13_1q21.1"),
    chr8_8p23.1 = list(chr = "chr8", start = 7833176, end = 11501098, name = "chm13_8p23.1")
  ),
  
  # Genomic regions of interest (GRCh38 coordinates)
  grch38_regions = list(
    chr22_22q11 = list(chr = "chr22", start = 19022279, end = 23380258, name = "grch38_22q11"),
    chr15_anglemans = list(chr = "chr15", start = 23123712, end = 28193120, name = "grch38_anglemans"),
    chr1_1q21.1 = list(chr = "chr1", start = 145686995, end = 148411223, name = "grch38_1q21.1"),
    chr16_16p12.1 = list(chr = "chr16", start = 21398763, end = 30188533, name = "grch38_16p12.1"),
    chr8_8p23.1 = list(chr = "chr8", start = 8242533, end = 11907120, name = "grch38_8p23.1")
  )
)

CONFIG$figure6_regions <- c("chr15_anglemans", "chr22_22q11")

# Layout r-values for data bands within data.panel = 1
CONFIG$plot$layout <- list(
  var_band_start = 0,
  var_band_stop = 0.1,
  data_1_start = 0.175,
  data_1_end = 0.55,
  data_2_start = 0.625,
  data_2_end = 1
)

# Convert linewidth from points to R units
CONFIG$plot$linewidth <- CONFIG$plot$linewidth_pts / 0.75

# ============================================================================
# UTILITY FUNCTIONS (from refactored script)
# ============================================================================

# Function for named group splitting
named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))
  
  grouped %>%
    group_split() %>%
    rlang::set_names(names)
}

# Clamp values between min and max
clamp <- function(x, min, max) {
  case_when(
    x < min ~ min,
    x > max ~ max,
    .default = x
  )
}

strip_unit <- function(value) {
  as.numeric(gsub("[^0-9.-]", "", value))
}

expand_figure6_source_canvas <- function(svg_path, left_pad_pt = 12, bottom_pad_pt = 4) {
  svg_lines <- readLines(svg_path, warn = FALSE)
  root_idx <- grep("<svg ", svg_lines, fixed = TRUE)[1]
  if (is.na(root_idx)) {
    stop(paste("Could not find root <svg> tag in", svg_path))
  }

  root_line <- svg_lines[root_idx]
  width_pt <- strip_unit(sub(".*width='([^']+)'.*", "\\1", root_line))
  height_pt <- strip_unit(sub(".*height='([^']+)'.*", "\\1", root_line))
  padded_width_pt <- width_pt + left_pad_pt
  padded_height_pt <- height_pt + bottom_pad_pt

  svg_lines[root_idx] <- sub(
    "width='[^']+'",
    sprintf("width='%.2fpt'", padded_width_pt),
    svg_lines[root_idx]
  )
  svg_lines[root_idx] <- sub(
    "height='[^']+'",
    sprintf("height='%.2fpt'", padded_height_pt),
    svg_lines[root_idx]
  )
  svg_lines[root_idx] <- sub(
    "viewBox='[^']+'",
    sprintf("viewBox='%.2f 0 %.2f %.2f'", -left_pad_pt, padded_width_pt, padded_height_pt),
    svg_lines[root_idx]
  )

  background_idx <- grep("<rect width='100%' height='100%'", svg_lines, fixed = TRUE)
  for (idx in background_idx) {
    svg_lines[idx] <- sub(
      "<rect width='100%' height='100%'",
      sprintf("<rect x='%.2f' width='%.2f' height='100%%'", -left_pad_pt, padded_width_pt),
      svg_lines[idx],
      fixed = TRUE
    )
  }

  full_clip_rect <- sprintf(
    "<rect x='0.00' y='0.00' width='%.2f' height='%.2f' />",
    width_pt,
    height_pt
  )
  full_clip_rect_replacement <- sprintf(
    "<rect x='%.2f' y='0.00' width='%.2f' height='%.2f' />",
    -left_pad_pt,
    padded_width_pt,
    padded_height_pt
  )
  svg_lines <- sub(full_clip_rect, full_clip_rect_replacement, svg_lines, fixed = TRUE)

  inner_clip_rect <- sprintf(
    "<rect x='0.72' y='0.72' width='%.2f' height='%.2f' />",
    width_pt - 1.44,
    height_pt - 1.44
  )
  inner_clip_rect_replacement <- sprintf(
    "<rect x='%.2f' y='0.72' width='%.2f' height='%.2f' />",
    0.72 - left_pad_pt,
    padded_width_pt - 1.44,
    padded_height_pt - 1.44
  )
  svg_lines <- sub(inner_clip_rect, inner_clip_rect_replacement, svg_lines, fixed = TRUE)

  writeLines(svg_lines, svg_path)
}

# Clip overlapping genomic ranges
clip_data <- function(gr1, gr2) {
  o <- findOverlaps(gr1, gr2)
  grl1 <- split(gr1[queryHits(o)], seq_len(length(o)))
  grl2 <- split(gr2[subjectHits(o)], seq_len(length(o)))
  foo <- function(x, y) {
    rv <- x
    start(rv) <- max(start(x), start(y))
    end(rv) <- min(end(x), end(y))
    return(rv)
  }
  unlist(mendoapply(foo, grl1, y = grl2))
}

load_bed_ranges <- function(path) {
  bed_df <- read_tsv(
    path,
    col_names = c("chrom", "start", "end"),
    col_types = cols(
      chrom = col_character(),
      start = col_integer(),
      end = col_integer()
    ),
    comment = "#",
    progress = FALSE,
    show_col_types = FALSE
  )

  makeGRangesFromDataFrame(
    bed_df,
    seqnames.field = "chrom",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE,
    starts.in.df.are.0based = TRUE
  )
}

clip_and_merge_highlight_data <- function(highlight_data, region) {
  clipped_data <- clip_data(highlight_data, region)
  if (length(clipped_data) == 0) {
    return(clipped_data)
  }

  GenomicRanges::reduce(clipped_data, ignore.strand = TRUE)
}

add_highlight_rectangles <- function(kp, highlight_data, fill_color, r0, r1, data_panel, config) {
  if (length(highlight_data) == 0) {
    return(kp)
  }

  kpRect(
    kp,
    data = highlight_data,
    y0 = 0,
    y1 = 1,
    col = fill_color,
    border = NA,
    r0 = r0,
    r1 = r1,
    data.panel = data_panel,
    clipping = FALSE,
    lwd = config$plot$linewidth
  )
}

# Create GRanges objects from region list
create_genomic_ranges <- function(regions_list, slop = 0) {
  ranges_list <- list()
  names_list <- list()
  
  for (i in seq_along(regions_list)) {
    region <- regions_list[[i]]
    ranges_list[[i]] <- toGRanges(data.frame(
      region$chr, 
      region$start - slop, 
      region$end + slop
    ))
    names_list[[i]] <- region$name
  }
  
  return(list(ranges = ranges_list, names = names_list))
}

create_padded_region <- function(region_info, slop = 0) {
  toGRanges(data.frame(
    region_info$chr,
    region_info$start - slop,
    region_info$end + slop
  ))
}

max_variants_in_region <- function(data, region, label) {
  region_data <- clip_data(data, region)

  if (length(region_data) == 0) {
    stop(paste("No variant density data found for", label))
  }

  max(region_data$variants_per_1k, na.rm = TRUE)
}

calculate_figure6_max_var_density <- function(grch38_data, chm13_data, config) {
  figure6_max_values <- c()

  for (region_name in config$figure6_regions) {
    grch38_region <- create_padded_region(config$grch38_regions[[region_name]], config$plot$slop)
    chm13_region <- create_padded_region(config$chm13_regions[[region_name]], config$plot$slop)

    figure6_max_values <- c(
      figure6_max_values,
      max_variants_in_region(grch38_data, grch38_region, paste("GRCh38", region_name)),
      max_variants_in_region(chm13_data, chm13_region, paste("CHM13v2.0", region_name))
    )
  }

  max(figure6_max_values, na.rm = TRUE)
}

format_colorbar_breaks <- function(max_var_density, num_breaks = 6) {
  sprintf("%.2f", seq(0, max_var_density, length.out = num_breaks))
}

# `kpArea` and `kpHeatmap` handle GRanges x-coordinates differently. We reuse
# one constructor so the script can keep the original rolling windows for the
# area plots and a second, derived coordinate set for the density heatmap.
rolling_df_to_granges_map <- function(df, start_field = "start", end_field = "end", drop_cols = character()) {
  df %>%
    select(-any_of(drop_cols)) %>%
    within(method <- paste(method_of_phasing, ground_truth_data_source, sep = "_")) %>%
    named_group_split(method) %>%
    map(
      makeGRangesFromDataFrame,
      seqnames.field = "chrom",
      start.field = start_field,
      end.field = end_field,
      keep.extra.columns = TRUE,
      ignore.strand = TRUE,
      starts.in.df.are.0based = TRUE
    )
}

# The rolling parquet currently stores the window midpoint in `start` and
# computes `end` from that midpoint. Shift both bounds back by half a window so
# the plotted intervals match the genomic span summarized by each row.
recenter_rolling_window_bounds <- function(df) {
  window_widths <- unique(df$end - df$start)
  if (length(window_widths) != 1) {
    stop(
      paste(
        "Expected a single rolling window width in",
        CONFIG$rolling_data_file,
        "but found:",
        paste(sort(window_widths), collapse = ", ")
      )
    )
  }

  window_width <- window_widths[[1]]
  if (window_width %% 2 != 0) {
    stop(paste("Expected an even rolling window width but found", window_width))
  }

  half_window <- as.integer(window_width / 2)
  df %>%
    mutate(
      start = as.integer(start - half_window),
      end = as.integer(end - half_window)
    )
}

# The rolling parquet stores 500 kb summaries every 10 kb. `kpHeatmap` draws
# the full interval width, so plotting the raw windows would smear each sample
# across 500 kb and visually offset density gaps relative to the midpoint-based
# error-rate tracks. This helper condenses each row into a 10 kb tile centered
# on the original rolling-window midpoint.
create_variant_density_tiles <- function(df,
                                         group_cols = c("genome", "chrom", "method_of_phasing", "ground_truth_data_source")) {
  window_widths <- unique(df$end - df$start)
  if (length(window_widths) != 1) {
    stop(
      paste(
        "Expected a single rolling window width in",
        CONFIG$rolling_data_file,
        "but found:",
        paste(sort(window_widths), collapse = ", ")
      )
    )
  }

  window_width <- window_widths[[1]]
  density_dfs <- df %>%
    group_by(across(all_of(group_cols))) %>%
    group_split() %>%
    map(function(group_df) {
      group_df <- arrange(group_df, start)
      group_label <- paste(group_df[1, group_cols], collapse = " / ")

      positive_steps <- unique(diff(group_df$start))
      positive_steps <- positive_steps[positive_steps > 0]
      if (length(positive_steps) != 1) {
        stop(
          paste(
            "Expected one positive rolling start increment for",
            group_label,
            "but found:",
            paste(sort(positive_steps), collapse = ", ")
          )
        )
      }

      sampling_step <- positive_steps[[1]]
      tile_padding <- window_width - sampling_step
      if (tile_padding < 0 || tile_padding %% 2 != 0) {
        stop(
          paste(
            "Cannot derive midpoint-centered density tiles for",
            group_label,
            "with window width",
            window_width,
            "and sampling step",
            sampling_step
          )
        )
      }

      group_df <- mutate(
        group_df,
        sampling_step = sampling_step,
        density_start = as.integer(start + tile_padding / 2),
        density_end = as.integer(start + tile_padding / 2 + sampling_step)
      )

      if (any((group_df$density_end - group_df$density_start) != group_df$sampling_step)) {
        stop(paste("Derived density tiles have the wrong width for", group_label))
      }
      if (any((group_df$density_start + group_df$density_end) != (group_df$start + group_df$end))) {
        stop(paste("Derived density tiles do not preserve window midpoints for", group_label))
      }

      group_df
    })

  density_df <- bind_rows(density_dfs)
  if (n_distinct(density_df$sampling_step) != 1) {
    stop("Expected a single sampling step across all rolling-stat groups")
  }

  density_df
}

# Setup plot parameters for single genome half
setup_single_plot_params <- function() {
  plot.params <- getDefaultPlotParams(plot.type = 1)
  plot.params$leftmargin <- 0.2
  plot.params$topmargin <- 5
  plot.params$bottommargin <- 5
  plot.params$data1height <- 200
  plot.params$data1inmargin <- 15
  plot.params$ideogramheight <- 30
  return(plot.params)
}

# ============================================================================
# SINGLE GENOME HALF GENERATION FUNCTIONS
# ============================================================================

# Add the regional data tracks to a selected karyoploteR panel
add_regional_panel_tracks <- function(kp, region, region_data, density_region_data, cnv_data, segdup_data,
                                      ser_color, gt_color, max_var_density,
                                      data_panel, config,
                                      add_track_labels = TRUE,
                                      add_base_numbers = TRUE,
                                      var_label_margin = config$plot$regional_var_label_margin,
                                      data_label_margin = config$plot$regional_data_label_margin) {
  segdup_region_data <- clip_and_merge_highlight_data(segdup_data, region)
  cnv_region_data <- clip_data(cnv_data, region)

  kp <- add_highlight_rectangles(
    kp,
    highlight_data = segdup_region_data,
    fill_color = config$plot$segdup_highlight_color,
    r0 = config$plot$layout$data_1_start,
    r1 = config$plot$layout$data_1_end,
    data_panel = data_panel,
    config = config
  )

  kp <- add_highlight_rectangles(
    kp,
    highlight_data = segdup_region_data,
    fill_color = config$plot$segdup_highlight_color,
    r0 = config$plot$layout$data_2_start,
    r1 = 1,
    data_panel = data_panel,
    config = config
  )

  kp <- add_highlight_rectangles(
    kp,
    highlight_data = cnv_region_data,
    fill_color = config$plot$cnv_highlight_color,
    r0 = config$plot$layout$data_1_start,
    r1 = config$plot$layout$data_1_end,
    data_panel = data_panel,
    config = config
  )

  kp <- add_highlight_rectangles(
    kp,
    highlight_data = cnv_region_data,
    fill_color = config$plot$cnv_highlight_color,
    r0 = config$plot$layout$data_2_start,
    r1 = 1,
    data_panel = data_panel,
    config = config
  )
  
  # `kpHeatmap` uses the full interval width, so feed it the condensed 10 kb
  # tiles instead of the original 500 kb rolling windows.
  kp <- kpHeatmap(kp, data = density_region_data,
                  r0 = 0, r1 = config$plot$layout$var_band_stop, 
                  y = clamp(density_region_data$variants_per_1k, 0, max_var_density),
                  ymax = max_var_density, 
                  colors = config$plot$var_density_palette, 
                  data.panel = data_panel, clipping = FALSE, lwd = config$plot$linewidth)
  
  # Add switch error rate area plot
  kp <- kpArea(kp, data = region_data, data.panel = data_panel, 
               r0 = config$plot$layout$data_1_start, 
               r1 = config$plot$layout$data_1_end, 
               y = clamp(region_data$switch_error_rate, 0, config$plot$ser_ymax), 
               ymax = config$plot$ser_ymax, col = ser_color, 
               clipping = FALSE, lwd = config$plot$linewidth)
  
  # Add genotype error rate area plot
  kp <- kpArea(kp, data = region_data, data.panel = data_panel, 
               r0 = config$plot$layout$data_2_start, r1 = 1, 
               y = clamp(region_data$gt_error_rate, 0, config$plot$gt_ymax), 
               ymax = config$plot$gt_ymax, col = gt_color, 
               clipping = FALSE, lwd = config$plot$linewidth)
  
  # Add axes
  kpAxis(kp, numticks = config$plot$num_ticks, data.panel = data_panel, 
         r0 = config$plot$layout$data_1_start, r1 = config$plot$layout$data_1_end, 
         ymin = 0, ymax = config$plot$ser_ymax, 
         cex = (config$plot$percent_size_pts / config$plot$default_pts), 
         labels = paste0(seq(0, config$plot$ser_ymax, length.out = config$plot$num_ticks), "%"), 
         lwd = config$plot$linewidth)
  
  kpAxis(kp, numticks = config$plot$num_ticks, data.panel = data_panel, 
         r0 = config$plot$layout$data_2_start, r1 = 1, 
         ymin = 0, ymax = config$plot$gt_ymax, 
         cex = (config$plot$percent_size_pts / config$plot$default_pts), 
         labels = paste0(seq(0, config$plot$gt_ymax, length.out = config$plot$num_ticks), "%"), 
         lwd = config$plot$linewidth)
  
  if (add_track_labels) {
    kpAddLabels(kp, labels = "# Variants/1kb", 
                r0 = 0, r1 = config$plot$layout$var_band_stop, data.panel = data_panel, 
                cex = (config$plot$label_size_pts / config$plot$default_pts), 
                label.margin = var_label_margin, lwd = config$plot$linewidth)
    
    kpAddLabels(kp, labels = "Switch\nError\nRate", 
                r0 = config$plot$layout$data_1_start, r1 = config$plot$layout$data_1_end, 
                data.panel = data_panel, 
                cex = (config$plot$label_size_pts / config$plot$default_pts), 
                label.margin = data_label_margin, lwd = config$plot$linewidth)
    
    kpAddLabels(kp, labels = "Genotype\nError\nRate", 
                r0 = config$plot$layout$data_2_start, r1 = 1, data.panel = data_panel, 
                cex = (config$plot$label_size_pts / config$plot$default_pts), 
                label.margin = data_label_margin, lwd = config$plot$linewidth)
  }
  
  if (add_base_numbers) {
    kpAddBaseNumbers(kp, tick.dist = 1000000, minor.tick.dist = 200000, 
                     cex = (config$plot$base_number_size_pts / config$plot$default_pts), 
                     lwd = config$plot$linewidth, data.panel = data_panel)
  }
  
  return(kp)
}

# Create single genome ideogram plot (regional)
create_single_genome_regional <- function(region, data, density_data, cytobands, cnv_data, segdup_data, genome,
                                        ser_color, gt_color, max_var_density,
                                        plot_params, config,
                                        var_label_margin = config$plot$regional_var_label_margin,
                                        data_label_margin = config$plot$regional_data_label_margin) {
  region_data <- clip_data(data, region)
  density_region_data <- clip_data(density_data, region)
  
  # Create karyotype plot
  kp <- plotKaryotype(
    genome = genome,
    zoom = region,
    plot.type = 1,
    plot.params = plot_params,
    cytobands = clip_data(cytobands, region),
    lwd = config$plot$linewidth,
    cex = (config$plot$chrom_size_pts / config$plot$default_pts)
  )
  
  kp <- add_regional_panel_tracks(
    kp = kp,
    region = region,
    region_data = region_data,
    density_region_data = density_region_data,
    cnv_data = cnv_data,
    segdup_data = segdup_data,
    ser_color = ser_color,
    gt_color = gt_color,
    max_var_density = max_var_density,
    data_panel = 1,
    config = config,
    add_track_labels = TRUE,
    add_base_numbers = TRUE,
    var_label_margin = var_label_margin,
    data_label_margin = data_label_margin
  )
  
  return(kp)
}

# Create single genome whole-genome plot
create_single_genome_whole <- function(data, density_data, cytobands, cnv_data, genome,
                                     ser_color, max_var_density, plot_params, config) {

  # Create karyotype plot
  kp <- plotKaryotype(genome = genome, plot.type = 1, plot.params = plot_params,
                      cytobands = cytobands)

  # Add CNV rectangles
  kp <- kpRect(kp, data = cnv_data, y0 = 0, y1 = 1, col = "#FFDDDD",
               border = NA, r0 = 0.2, r1 = 1, data.panel = 1)

  # Add variant density heatmap
  kp <- kpHeatmap(kp, data = density_data, r0 = 0, r1 = 0.2,
                  y = clamp(density_data$variants_per_1k, 0, max_var_density),
                  ymax = max_var_density, colors = config$plot$var_density_palette,
                  data.panel = 1, clipping = TRUE)

  # Add switch error rate area plot
  kp <- kpArea(kp, data = data, data.panel = 1,
               r0 = 0.3, r1 = 1,
               y = clamp(data$switch_error_rate, 0, config$plot$ser_ymax),
               ymax = config$plot$ser_ymax, col = ser_color,
               clipping = TRUE)

  return(kp)
}

# ============================================================================
# SVG EXPORT FUNCTIONS
# ============================================================================

# Export a single-genome regional plot to SVG
export_svg_regional <- function(region, data, density_data, cytobands, cnv_data, genome,
                                segdup_data,
                                ser_color, gt_color, max_var_density,
                                plot_params, config, out_file,
                                var_label_margin = config$plot$regional_var_label_margin,
                                data_label_margin = config$plot$regional_data_label_margin,
                                width_inches = config$plot$width_inches,
                                height_inches = config$plot$height_inches) {
  svglite(out_file,
          width    = width_inches,
          height   = height_inches,
          pointsize = config$plot$default_pts,
          system_fonts = list(sans = "Arial"))

  create_single_genome_regional(
    region       = region,
    data         = data,
    density_data = density_data,
    cytobands    = cytobands,
    cnv_data     = cnv_data,
    segdup_data  = segdup_data,
    genome       = genome,
    ser_color    = ser_color,
    gt_color     = gt_color,
    max_var_density = max_var_density,
    plot_params  = plot_params,
    config       = config,
    var_label_margin = var_label_margin,
    data_label_margin = data_label_margin
  )

  dev.off()
}

# Export a single-genome whole-genome plot to SVG
export_svg_whole <- function(data, density_data, cytobands, cnv_data, genome,
                             ser_color, max_var_density, plot_params,
                             config, out_file) {
  svglite(out_file,
          width    = config$plot$width_inches * 2,
          height   = config$plot$height_inches,
          pointsize = config$plot$default_pts,
          system_fonts = list(sans = "Arial"))

  create_single_genome_whole(
    data         = data,
    density_data = density_data,
    cytobands    = cytobands,
    cnv_data     = cnv_data,
    genome       = genome,
    ser_color    = ser_color,
    max_var_density = max_var_density,
    plot_params  = plot_params,
    config       = config
  )

  dev.off()
}

export_svg_figure6_panel <- function(region, data, density_data, cytobands, cnv_data, segdup_data, genome,
                                     ser_color, gt_color, max_var_density,
                                     plot_params, config, out_file) {
  export_svg_regional(
    region = region,
    data = data,
    density_data = density_data,
    cytobands = cytobands,
    cnv_data = cnv_data,
    segdup_data = segdup_data,
    genome = genome,
    ser_color = ser_color,
    gt_color = gt_color,
    max_var_density = max_var_density,
    plot_params = plot_params,
    config = config,
    out_file = out_file,
    width_inches = config$plot$figure6_panel_width_inches,
    height_inches = config$plot$figure6_panel_height_inches
  )

  expand_figure6_source_canvas(out_file)
}

export_svg_colorbar <- function(out_file, max_var_density, config) {
  legend_labels <- format_colorbar_breaks(max_var_density)
  width_pt <- config$plot$figure6_colorbar_width_inches * 72
  height_pt <- config$plot$figure6_colorbar_height_inches * 72

  # Span the axis between the centers of the two Figure 6 panel columns after
  # stitching, while keeping the legend centered in its own SVG canvas.
  bar_width <- config$plot$figure6_colorbar_axis_width_pt
  bar_height <- 7.56
  bar_x <- (width_pt - bar_width) / 2
  bar_y <- 10.80
  tick_y1 <- bar_y + bar_height
  tick_y2 <- 19.80
  tick_label_y <- 26.27
  title_y <- 35.99
  title_cex <- config$plot$chrom_size_pts / config$plot$default_pts
  n_steps <- 100

  svglite(out_file,
          width = config$plot$figure6_colorbar_width_inches,
          height = config$plot$figure6_colorbar_height_inches,
          pointsize = config$plot$default_pts,
          system_fonts = list(sans = "Arial"))

  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, width_pt), ylim = c(height_pt, 0), xaxs = "i", yaxs = "i")

  palette_colors <- colorRampPalette(config$plot$var_density_palette)(n_steps)
  x_breaks <- seq(bar_x, bar_x + bar_width, length.out = n_steps + 1)
  rect(
    xleft = x_breaks[-length(x_breaks)],
    ybottom = bar_y,
    xright = x_breaks[-1],
    ytop = bar_y + bar_height,
    col = palette_colors,
    border = NA,
    xpd = NA
  )

  rect(
    xleft = bar_x,
    ybottom = bar_y,
    xright = bar_x + bar_width,
    ytop = bar_y + bar_height,
    border = "black",
    lwd = 0.75,
    xpd = NA
  )

  tick_positions <- seq(bar_x, bar_x + bar_width, length.out = length(legend_labels))
  segments(tick_positions, tick_y1, tick_positions, tick_y2, lwd = 0.75, xpd = NA)
  text(tick_positions, tick_label_y, labels = legend_labels, xpd = NA)

  text(
    bar_x + (bar_width / 2),
    title_y,
    labels = "Variants/kb",
    cex = title_cex,
    xpd = NA
  )

  dev.off()
}

# ============================================================================
# DATA LOADING AND PREPROCESSING (from refactored script)
# ============================================================================

cat("Loading and preprocessing data...\n")

# Load and preprocess rolling data
rolling_data <- read_parquet(CONFIG$rolling_data_file) %>%
  recenter_rolling_window_bounds()
rolling_data$switch_error_rate <- replace_na(rolling_data$n_switch_errors/rolling_data$n_checked, 0) * 100
rolling_data$gt_error_rate <- replace_na(rolling_data$n_gt_errors/rolling_data$n_gt_checked, 0) * 100
# Keep the corrected rolling-window coordinates for the error-rate areas and
# build a second, condensed coordinate set for the density heatmap.
rolling_density_data <- create_variant_density_tiles(rolling_data)

# Split data by genome
GRCh38_rolling_data <- rolling_data[rolling_data$genome == 'GRCh38', ]
T2T_rolling_data <- rolling_data[rolling_data$genome == 'CHM13v2.0', ]
GRCh38_density_data <- rolling_density_data[rolling_density_data$genome == 'GRCh38', ]
T2T_density_data <- rolling_density_data[rolling_density_data$genome == 'CHM13v2.0', ]

# Process data into GRanges objects
grch38_data <- rolling_df_to_granges_map(GRCh38_rolling_data)
chm13_data <- rolling_df_to_granges_map(T2T_rolling_data)
grch38_density_gr <- rolling_df_to_granges_map(
  GRCh38_density_data,
  start_field = "density_start",
  end_field = "density_end",
  drop_cols = c("start", "end")
)
chm13_density_gr <- rolling_df_to_granges_map(
  T2T_density_data,
  start_field = "density_start",
  end_field = "density_end",
  drop_cols = c("start", "end")
)

# Load CNV data
cnv <- read_tsv(CONFIG$cnv_file)
cnv[cnv$start_grch38 > cnv$end_grch38, c('start_grch38', 'end_grch38')] <- 
  cnv[cnv$start_grch38 > cnv$end_grch38, c('end_grch38', 'start_grch38')]
cnv[cnv$start_chm13 > cnv$end_chm13, c('start_chm13', 'end_chm13')] <- 
  cnv[cnv$start_chm13 > cnv$end_chm13, c('end_chm13', 'start_chm13')]

cnv_chm13 <- makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE, 
                                      start.field = 'start_chm13', 
                                      end.field = 'end_chm13',
                                      strand.field = 'strand_chm13', 
                                      starts.in.df.are.0based = TRUE)

cnv_grch38 <- makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE, 
                                       start.field = 'start_grch38', 
                                       end.field = 'end_grch38',
                                       strand.field = 'strand_grch38', 
                                       starts.in.df.are.0based = TRUE)

# Load segmental duplication BED data
segdup_chm13 <- load_bed_ranges(CONFIG$t2t_segdup_file)
segdup_grch38 <- load_bed_ranges(CONFIG$grch38_segdup_file)

# Load cytoband data
chm13_cytobands <- makeGRangesFromDataFrame(read_tsv(CONFIG$t2t_cytobands_header_file), 
                                            keep.extra.columns = TRUE)
grch38_cytobands <- makeGRangesFromDataFrame(read_tsv(CONFIG$grch38_cytobands_header_file), 
                                             keep.extra.columns = TRUE)

# Extract specific dataset
chm13_3202_vs_HPRC_data <- chm13_data[['1kgp_variation_phased_with_reference_panel_HPRC_samples']]
grch38_3202_vs_HPRC_data <- grch38_data[['1kgp_variation_phased_with_reference_panel_HPRC_samples']]
chm13_3202_vs_HPRC_density <- chm13_density_gr[['1kgp_variation_phased_with_reference_panel_HPRC_samples']]
grch38_3202_vs_HPRC_density <- grch38_density_gr[['1kgp_variation_phased_with_reference_panel_HPRC_samples']]

# Calculate maximum variant density
global_max_var_density <- round(max(max(grch38_3202_vs_HPRC_data$variants_per_1k),
                                   max(chm13_3202_vs_HPRC_data$variants_per_1k)))
figure6_max_var_density <- calculate_figure6_max_var_density(
  grch38_data = grch38_3202_vs_HPRC_density,
  chm13_data = chm13_3202_vs_HPRC_density,
  config = CONFIG
)

# ============================================================================
# MAIN EXECUTION
# ============================================================================

cat("Setting up output directories...\n")

# Create output directory
dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)

# Shared plot parameters for regional plots
regional_plot_params <- setup_single_plot_params()

colorbar_svg <- file.path(CONFIG$output_dir, "colorbar.svg")
export_svg_colorbar(colorbar_svg, figure6_max_var_density, CONFIG)
cat(paste("Saved:", colorbar_svg, "\n"))

cat("Generating regional SVGs...\n")

for (region_name in names(CONFIG$grch38_regions)) {
  cat(paste("Processing region:", region_name, "\n"))

  grch38_region_info <- CONFIG$grch38_regions[[region_name]]
  t2t_region_info    <- CONFIG$chm13_regions[[region_name]]

  grch38_region <- toGRanges(data.frame(
    grch38_region_info$chr,
    grch38_region_info$start - CONFIG$plot$slop,
    grch38_region_info$end   + CONFIG$plot$slop
  ))

  t2t_region <- toGRanges(data.frame(
    t2t_region_info$chr,
    t2t_region_info$start - CONFIG$plot$slop,
    t2t_region_info$end   + CONFIG$plot$slop
  ))

  suffix <- gsub("^grch38_|^chm13_", "", region_name)

  grch38_svg <- file.path(CONFIG$output_dir, paste0("grch38_", suffix, ".svg"))
  export_svg_regional(
    region       = grch38_region,
    data         = grch38_3202_vs_HPRC_data,
    density_data = grch38_3202_vs_HPRC_density,
    cytobands    = grch38_cytobands,
    cnv_data     = cnv_grch38,
    segdup_data  = segdup_grch38,
    genome       = "hg38",
    ser_color    = CONFIG$plot$grch38_ser_color,
    gt_color     = CONFIG$plot$grch38_gt_color,
    max_var_density = global_max_var_density,
    plot_params  = regional_plot_params,
    config       = CONFIG,
    out_file     = grch38_svg
  )
  cat(paste("  Saved:", grch38_svg, "\n"))

  if (region_name %in% CONFIG$figure6_regions) {
    grch38_figure6_svg <- file.path(CONFIG$output_dir, paste0("figure6_grch38_", suffix, ".svg"))
    export_svg_figure6_panel(
      region       = grch38_region,
      data         = grch38_3202_vs_HPRC_data,
      density_data = grch38_3202_vs_HPRC_density,
      cytobands    = grch38_cytobands,
      cnv_data     = cnv_grch38,
      segdup_data  = segdup_grch38,
      genome       = "hg38",
      ser_color    = CONFIG$plot$grch38_ser_color,
      gt_color     = CONFIG$plot$grch38_gt_color,
      max_var_density = figure6_max_var_density,
      plot_params  = regional_plot_params,
      config       = CONFIG,
      out_file     = grch38_figure6_svg
    )
    cat(paste("  Saved:", grch38_figure6_svg, "\n"))
  }

  t2t_svg <- file.path(CONFIG$output_dir, paste0("t2t_", suffix, ".svg"))
  export_svg_regional(
    region       = t2t_region,
    data         = chm13_3202_vs_HPRC_data,
    density_data = chm13_3202_vs_HPRC_density,
    cytobands    = chm13_cytobands,
    cnv_data     = cnv_chm13,
    segdup_data  = segdup_chm13,
    genome       = "hs1",
    ser_color    = CONFIG$plot$chm13_ser_color,
    gt_color     = CONFIG$plot$chm13_gt_color,
    max_var_density = global_max_var_density,
    plot_params  = regional_plot_params,
    config       = CONFIG,
    out_file     = t2t_svg
  )
  cat(paste("  Saved:", t2t_svg, "\n"))

  if (region_name %in% CONFIG$figure6_regions) {
    t2t_figure6_svg <- file.path(CONFIG$output_dir, paste0("figure6_t2t_", suffix, ".svg"))
    export_svg_figure6_panel(
      region       = t2t_region,
      data         = chm13_3202_vs_HPRC_data,
      density_data = chm13_3202_vs_HPRC_density,
      cytobands    = chm13_cytobands,
      cnv_data     = cnv_chm13,
      segdup_data  = segdup_chm13,
      genome       = "hs1",
      ser_color    = CONFIG$plot$chm13_ser_color,
      gt_color     = CONFIG$plot$chm13_gt_color,
      max_var_density = figure6_max_var_density,
      plot_params  = regional_plot_params,
      config       = CONFIG,
      out_file     = t2t_figure6_svg
    )
    cat(paste("  Saved:", t2t_figure6_svg, "\n"))
  }
}

cat("Generating whole-genome SVGs...\n")

whole_plot_params <- getDefaultPlotParams(plot.type = 1)
whole_plot_params$leftmargin <- 0.2

grch38_whole_svg <- file.path(CONFIG$output_dir, "grch38_whole_genome.svg")
export_svg_whole(
  data         = grch38_3202_vs_HPRC_data,
  density_data = grch38_3202_vs_HPRC_density,
  cytobands    = grch38_cytobands,
  cnv_data     = cnv_grch38,
  genome       = "hg38",
  ser_color    = CONFIG$plot$grch38_ser_color,
  max_var_density = global_max_var_density,
  plot_params  = whole_plot_params,
  config       = CONFIG,
  out_file     = grch38_whole_svg
)
cat(paste("Saved:", grch38_whole_svg, "\n"))

t2t_whole_svg <- file.path(CONFIG$output_dir, "t2t_whole_genome.svg")
export_svg_whole(
  data         = chm13_3202_vs_HPRC_data,
  density_data = chm13_3202_vs_HPRC_density,
  cytobands    = chm13_cytobands,
  cnv_data     = cnv_chm13,
  genome       = "hs1",
  ser_color    = CONFIG$plot$chm13_ser_color,
  max_var_density = global_max_var_density,
  plot_params  = whole_plot_params,
  config       = CONFIG,
  out_file     = t2t_whole_svg
)
cat(paste("Saved:", t2t_whole_svg, "\n"))

cat("Figure 6 SVG export completed successfully!\n")
cat("Output files saved to:", CONFIG$output_dir, "\n")
cat("Run 'python stitch_svgs.py --batch", CONFIG$output_dir, "' to combine pairs.\n")
