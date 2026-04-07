# ASV relative abundance scatterplots: Cop_ctr (x) vs Cop_pro (y) by day
# Input:
#   - ASV_table.csv in project root
#   - ASVtable_0.8threshold.csv in project root
# Output:
#   - asv_relative_abundance_scatterplots_cop_ctr_vs_cop_pro_D16.png in results/
#   - asv_relative_abundance_scatterplots_cop_ctr_vs_cop_pro_D21.png in results/
#   - asv_relative_abundance_scatterplots_cop_ctr_vs_cop_pro_D24.png in results/

library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tibble)
library(ggplot2)
library(ggrepel)

results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

count_table_path <- "ASV_table.csv"
taxonomy_table_path <- "ASVtable_0.8threshold.csv"
days_to_plot <- c("D16", "D21", "D24")
base_family <- "Times New Roman"
label_max_n <- 8
label_outlier_percent <- 2
fixed_axis_limit <- 30

day_colors <- c(
  "D16" = "#FFD84A",
  "D21" = "#90C95A",
  "D24" = "#E58FA2"
)

make_pretty_axis_limit <- function(max_value) {
  if (!is.finite(max_value) || max_value <= 0) {
    return(1)
  }

  if (max_value <= 10) {
    step <- 2
  } else if (max_value <= 25) {
    step <- 5
  } else if (max_value <= 50) {
    step <- 10
  } else {
    step <- 20
  }

  ceiling(max_value / step) * step
}

make_axis_breaks <- function(axis_limit) {
  if (!is.finite(axis_limit) || axis_limit <= 0) {
    return(c(0, 1))
  }

  if (axis_limit <= 10) {
    step <- 2
  } else if (axis_limit <= 25) {
    step <- 5
  } else if (axis_limit <= 50) {
    step <- 10
  } else {
    step <- 20
  }

  seq(0, axis_limit, by = step)
}

clean_bom <- function(x) {
  sub("^\ufeff", "", x)
}

detect_delimiter <- function(lines) {
  candidate_delimiters <- c(";", ",", "\t")
  lines <- clean_bom(lines)
  lines <- lines[!is.na(lines) & trimws(lines) != ""]

  if (length(lines) == 0) {
    return(";")
  }

  lines <- head(lines, 10)

  candidate_scores <- map_dbl(candidate_delimiters, function(delim) {
    field_counts <- lengths(strsplit(lines, delim, fixed = TRUE))
    valid_counts <- field_counts[field_counts > 1]

    if (length(valid_counts) == 0) {
      return(-Inf)
    }

    mode_count <- max(tabulate(match(valid_counts, unique(valid_counts))))
    mode_count * 1000 + median(valid_counts)
  })

  candidate_delimiters[which.max(candidate_scores)]
}

find_header_line <- function(lines) {
  cleaned_lines <- clean_bom(lines)
  header_idx <- which(
    str_detect(cleaned_lines, regex("(^|[;,\t])(asv|otu|feature).*(id)?($|[;,\t])", ignore_case = TRUE)) |
      str_detect(cleaned_lines, regex("Cop_(pro|ctr)_D(16|21|24)", ignore_case = TRUE))
  )

  if (length(header_idx) == 0) {
    stop("Could not detect a header row in ", count_table_path, ".")
  }

  header_idx[1]
}

read_count_table <- function(path) {
  raw_lines <- readr::read_lines(path, lazy = FALSE)
  header_line_idx <- find_header_line(raw_lines)
  delimiter <- detect_delimiter(raw_lines[header_line_idx:min(length(raw_lines), header_line_idx + 9)])

  table <- readr::read_delim(
    file = path,
    delim = delimiter,
    skip = header_line_idx - 1,
    show_col_types = FALSE,
    na = c("", "NA", "NaN"),
    trim_ws = TRUE
  )

  names(table) <- clean_bom(names(table))
  table
}

read_threshold_taxonomy <- function(path) {
  raw_lines <- readr::read_lines(path, lazy = FALSE)
  non_empty_lines <- raw_lines[trimws(raw_lines) != ""]

  if (length(non_empty_lines) == 0) {
    return(tibble(ASV_ID = character(), taxonomy_threshold = character()))
  }

  delimiter <- detect_delimiter(non_empty_lines)

  taxonomy_tbl <- readr::read_delim(
    file = path,
    delim = delimiter,
    col_names = FALSE,
    show_col_types = FALSE,
    trim_ws = TRUE
  )

  names(taxonomy_tbl)[1] <- "ASV_ID"
  if (ncol(taxonomy_tbl) < 2) {
    taxonomy_tbl$taxonomy_threshold <- NA_character_
  } else {
    names(taxonomy_tbl)[2] <- "taxonomy_threshold"
  }

  taxonomy_tbl %>%
    transmute(
      ASV_ID = clean_bom(as.character(ASV_ID)),
      taxonomy_threshold = na_if(as.character(taxonomy_threshold), "")
    )
}

is_numeric_like <- function(x) {
  if (is.numeric(x)) {
    return(TRUE)
  }

  x_chr <- trimws(as.character(x))
  x_chr <- x_chr[!is.na(x_chr) & x_chr != ""]

  if (length(x_chr) == 0) {
    return(FALSE)
  }

  all(str_detect(x_chr, "^-?[0-9]+(?:\\.[0-9]+)?$"))
}

detect_asv_column <- function(df) {
  name_scores <- case_when(
    str_detect(names(df), regex("^#?otu\\s*id$", ignore_case = TRUE)) ~ 100,
    str_detect(names(df), regex("\\basv\\b|\\bzotu\\b|\\botu\\b|feature", ignore_case = TRUE)) ~ 50,
    TRUE ~ 0
  )

  value_scores <- map_dbl(df, function(col) {
    values <- as.character(col)
    mean(str_detect(values, regex("^(asv|zotu|otu|feature)\\w*$", ignore_case = TRUE)), na.rm = TRUE)
  })

  combined_scores <- name_scores + value_scores
  best_idx <- which.max(combined_scores)

  if (length(best_idx) == 0 || combined_scores[best_idx] <= 0) {
    stop("Could not detect the ASV ID column.")
  }

  names(df)[best_idx]
}

looks_like_taxonomy <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr <- x_chr[!is.na(x_chr) & x_chr != ""]

  if (length(x_chr) == 0) {
    return(FALSE)
  }

  mean(
    str_detect(x_chr, regex("(^|[,;])\\s*[dpkcofgs]:", ignore_case = TRUE)) |
      str_detect(x_chr, regex("Bacteria|Archaea|Eukaryota", ignore_case = TRUE))
  ) > 0.5
}

has_confidence_scores <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr <- x_chr[!is.na(x_chr) & x_chr != ""]

  if (length(x_chr) == 0) {
    return(FALSE)
  }

  mean(
    str_detect(x_chr, "\\([0-9.]+\\)") |
      str_detect(x_chr, "(^|[;|,])\\s*[dpkcofgs]__[^;|,]*\\([0-9.]+\\)") |
      str_detect(x_chr, "(^|[;|,])\\s*[A-Za-z]+\\([^)]*[0-9.]+[^)]*\\)")
  ) > 0.1
}

detect_taxonomy_column <- function(df, asv_col) {
  candidate_names <- names(df)[
    str_detect(names(df), regex("tax|classif|consensus|lineage|assign", ignore_case = TRUE))
  ]

  candidate_names <- setdiff(candidate_names, asv_col)

  if (length(candidate_names) == 0) {
    return(NA_character_)
  }

  valid_candidates <- candidate_names[
    map_lgl(candidate_names, ~ looks_like_taxonomy(df[[.x]]) && !has_confidence_scores(df[[.x]]))
  ]

  if (length(valid_candidates) == 0) {
    return(NA_character_)
  }

  valid_candidates[1]
}

extract_genus <- function(taxonomy) {
  taxonomy_chr <- as.character(taxonomy)
  genus <- str_match(taxonomy_chr, "(^|[,;])\\s*g:([^,;]+)")[, 3]
  genus <- ifelse(is.na(genus) | genus == "", NA_character_, genus)
  genus
}

format_asv_label <- function(asv_id, taxonomy) {
  display_asv_id <- str_replace(as.character(asv_id), regex("^zotu", ignore_case = TRUE), "ASV")
  genus <- extract_genus(taxonomy)

  if_else(
    !is.na(genus),
    paste(display_asv_id, genus, sep = " - "),
    display_asv_id
  )
}

make_day_summary <- function(df, cop_pro_cols, cop_ctr_cols) {
  df %>%
    transmute(
      ASV_ID = ASV_ID,
      taxonomy = taxonomy,
      mean_Cop_pro = rowMeans(across(all_of(cop_pro_cols)), na.rm = TRUE) * 100,
      mean_Cop_ctr = rowMeans(across(all_of(cop_ctr_cols)), na.rm = TRUE) * 100
    ) %>%
    arrange(desc(mean_Cop_pro + mean_Cop_ctr), ASV_ID)
}

make_comparison_plot <- function(summary_tbl, day, panel_title = NULL) {
  plot_df <- summary_tbl %>%
    mutate(
      x_value = mean_Cop_ctr,
      y_value = mean_Cop_pro,
      combined_abundance = x_value + y_value,
      plot_label = format_asv_label(ASV_ID, taxonomy)
    )

  label_df <- plot_df %>%
    filter(pmax(x_value, y_value) >= label_outlier_percent) %>%
    arrange(desc(combined_abundance), ASV_ID) %>%
    slice_head(n = label_max_n)

  if (nrow(label_df) == 0) {
    label_df <- plot_df %>%
      arrange(desc(combined_abundance), ASV_ID) %>%
      slice_head(n = min(5, n()))
  }

  axis_limit <- fixed_axis_limit
  axis_breaks <- make_axis_breaks(axis_limit)
  point_color <- day_colors[[day]]

  plot <- ggplot(plot_df, aes(x = x_value, y = y_value)) +
    geom_point(
      shape = 21,
      size = 4.2,
      stroke = 0.6,
      fill = point_color,
      color = point_color,
      alpha = 0.85
    ) +
    geom_text_repel(
      data = label_df,
      aes(label = plot_label),
      family = base_family,
      size = 6.2,
      max.overlaps = Inf,
      seed = 42,
      box.padding = 0.6,
      point.padding = 0.5,
      force = 2.5,
      force_pull = 0.5,
      max.time = 3,
      max.iter = 20000,
      segment.color = "grey35",
      segment.size = 0.35,
      segment.alpha = 1,
      min.segment.length = 0
    ) +
    coord_cartesian(xlim = c(0, axis_limit), ylim = c(0, axis_limit), expand = TRUE) +
    scale_x_continuous(breaks = axis_breaks) +
    scale_y_continuous(breaks = axis_breaks) +
    labs(
      title = NULL,
      x = "Mean relative abundance (%) in Cop_ctr",
      y = "Mean relative abundance (%) in Cop_pro"
    ) +
    theme_classic(base_size = 16, base_family = base_family) +
    theme(
      axis.title = element_text(color = "black", size = 22),
      axis.text = element_text(color = "black", size = 22)
    )

  if (!is.null(panel_title)) {
    plot <- plot +
      labs(title = panel_title) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20, color = "black")
      )
  }

  plot
}

count_table_raw <- read_count_table(count_table_path)
asv_col <- detect_asv_column(count_table_raw)
taxonomy_col <- detect_taxonomy_column(count_table_raw, asv_col)
sample_cols <- names(count_table_raw)[map_lgl(count_table_raw, is_numeric_like)]
sample_cols <- setdiff(sample_cols, c(asv_col, taxonomy_col))

threshold_taxonomy <- read_threshold_taxonomy(taxonomy_table_path)

analysis_tbl <- count_table_raw %>%
  rename(ASV_ID = all_of(asv_col)) %>%
  mutate(ASV_ID = as.character(ASV_ID))

if (!is.na(taxonomy_col)) {
  analysis_tbl <- analysis_tbl %>%
    rename(taxonomy = all_of(taxonomy_col))
} else {
  analysis_tbl <- analysis_tbl %>%
    mutate(taxonomy = NA_character_)
}

analysis_tbl <- analysis_tbl %>%
  left_join(threshold_taxonomy, by = "ASV_ID") %>%
  mutate(taxonomy = coalesce(as.character(taxonomy), taxonomy_threshold)) %>%
  select(ASV_ID, taxonomy, all_of(sample_cols)) %>%
  mutate(across(all_of(sample_cols), as.numeric)) %>%
  mutate(
    across(
      all_of(sample_cols),
      ~ {
        column_sum <- sum(.x, na.rm = TRUE)
        if (is.na(column_sum) || column_sum <= 0) {
          rep(0, length(.x))
        } else {
          .x / column_sum
        }
      }
    )
  )

for (day in days_to_plot) {
  cop_pro_cols <- sample_cols[str_detect(sample_cols, paste0("^Cop_pro_", day, "_"))]
  cop_ctr_cols <- sample_cols[str_detect(sample_cols, paste0("^Cop_ctr_", day, "_"))]

  cat("Selected sample columns for", day, "\n")
  cat("  Cop_pro:", if (length(cop_pro_cols) > 0) paste(cop_pro_cols, collapse = ", ") else "None", "\n")
  cat("  Cop_ctr:", if (length(cop_ctr_cols) > 0) paste(cop_ctr_cols, collapse = ", ") else "None", "\n\n")

  if (length(cop_pro_cols) == 0 || length(cop_ctr_cols) == 0) {
    warning("Skipping ", day, " because one or more required sample groups were not detected.")
    next
  }

  summary_tbl <- make_day_summary(
    df = analysis_tbl,
    cop_pro_cols = cop_pro_cols,
    cop_ctr_cols = cop_ctr_cols
  )

  day_plot <- make_comparison_plot(summary_tbl, day)

  output_png <- file.path(
    results_dir,
    paste0("asv_relative_abundance_scatterplots_cop_ctr_vs_cop_pro_", day, ".png")
  )

  ggsave(
    filename = output_png,
    plot = day_plot,
    width = 8,
    height = 7.5,
    dpi = 300,
    bg = "white"
  )

  cat("Saved plot:", output_png, "\n\n")
}
