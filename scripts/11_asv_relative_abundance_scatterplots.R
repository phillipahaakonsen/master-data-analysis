# ASV relative abundance scatterplots by day
# Input:
#   - ASV_table.csv in project root (ASV count table)
#   - ASVtable_0.8threshold.csv in project root (fallback taxonomy table)
# Output:
#   - asv_relative_abundance_scatterplots_D16.png in results/
#   - asv_relative_abundance_scatterplots_D21.png in results/
#   - asv_relative_abundance_scatterplots_D24.png in results/
#   - asv_relative_abundance_summary_D16.csv in results/
#   - asv_relative_abundance_summary_D21.csv in results/
#   - asv_relative_abundance_summary_D24.csv in results/

library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(patchwork)

results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

count_table_path <- "ASV_table.csv"
taxonomy_table_path <- "ASVtable_0.8threshold.csv"
days_to_plot <- c("D16", "D21", "D24")
base_family <- "Times New Roman"
label_max_n <- 8
label_outlier_percent <- 2
x_axis_limit_percent <- 15

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

sample_type_colors <- c(
  "Cop_pro_D16" = "#20B7EA",
  "Cop_ctr_D16" = "#FDC300",
  "Cop_pro_D21" = "#4169E1",
  "Cop_ctr_D21" = "#FF9500",
  "Cop_pro_D24" = "#1624F2",
  "Cop_ctr_D24" = "#FF6200",
  "Prob_D16" = "#9ACD32",
  "Prob_D21" = "#86C400",
  "Prob_D24" = "#7FBF00"
)

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
      str_detect(cleaned_lines, regex("Cop_(pro|ctr)_D(16|21|24)|Prob_D(16|21|24)|Alg_D(16|21|24)", ignore_case = TRUE))
  )

  if (length(header_idx) == 0) {
    stop("Could not detect a header row in ", count_table_path, ".")
  }

  header_idx[1]
}

read_count_table <- function(path) {
  raw_lines <- readr::read_lines(path, lazy = FALSE)
  header_line_idx <- find_header_line(raw_lines)
  header_line <- clean_bom(raw_lines[header_line_idx])
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

make_day_summary <- function(df, day, prob_cols, cop_pro_cols, cop_ctr_cols) {
  df %>%
    transmute(
      ASV_ID = ASV_ID,
      taxonomy = taxonomy,
      mean_Prob = if (length(prob_cols) > 0) rowMeans(across(all_of(prob_cols)), na.rm = TRUE) * 100 else NA_real_,
      mean_Cop_pro = if (length(cop_pro_cols) > 0) rowMeans(across(all_of(cop_pro_cols)), na.rm = TRUE) * 100 else NA_real_,
      mean_Cop_ctr = if (length(cop_ctr_cols) > 0) rowMeans(across(all_of(cop_ctr_cols)), na.rm = TRUE) * 100 else NA_real_
    ) %>%
    arrange(desc(mean_Prob + mean_Cop_pro + mean_Cop_ctr), ASV_ID)
}

make_comparison_plot <- function(summary_tbl, x_col, y_col, x_label, y_label, title_text, x_color, y_color, y_axis_limit) {
  plot_df <- summary_tbl %>%
    mutate(
      x_value = .data[[x_col]],
      y_value = .data[[y_col]],
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
      slice_head(n = min(3, n()))
  }

  if (str_detect(title_text, "^D16\\s*: Prob vs Cop_pro$")) {
    label_df <- label_df %>%
      filter(!ASV_ID %in% c("Zotu1", "Zotu6")) %>%
      bind_rows(
        plot_df %>%
          filter(ASV_ID %in% c("Zotu60", "Zotu70"))
      ) %>%
      distinct(ASV_ID, .keep_all = TRUE) %>%
      arrange(desc(combined_abundance), ASV_ID)
  }

  if (str_detect(title_text, "^D21\\s*:")) {
    label_df <- label_df %>%
      filter(!ASV_ID %in% c("Zotu4", "Zotu8", "Zotu7"))
  }

  x_max_value <- x_axis_limit_percent
  y_max_value <- y_axis_limit
  if (!is.finite(y_max_value) || y_max_value <= 0) {
    y_max_value <- 1
  }

  x_breaks <- make_axis_breaks(x_max_value)
  y_breaks <- make_axis_breaks(y_max_value)
  label_nudge_x <- x_max_value * 0.03
  label_nudge_y <- y_max_value * 0.03
  y_min_value <- if (str_detect(title_text, "^(D16|D21)\\s*:")) -1 else 0

  label_df <- label_df %>%
    mutate(
      label_nudge_x = label_nudge_x,
      label_nudge_y = label_nudge_y
    )

  manual_label_df <- label_df %>%
    slice(0) %>%
    mutate(
      label_x = numeric(),
      label_y = numeric(),
      label_hjust = numeric(),
      label_vjust = numeric()
    )

  if (str_detect(title_text, "^D16\\s*: Prob vs Cop_ctr$")) {
    manual_label_df <- label_df %>%
      filter(ASV_ID %in% c("Zotu17", "Zotu29", "Zotu13")) %>%
      mutate(
        label_x = case_when(
          ASV_ID == "Zotu17" ~ x_value,
          ASV_ID == "Zotu29" ~ x_value + 0.18,
          ASV_ID == "Zotu13" ~ x_value + 0.55,
          TRUE ~ x_value + label_nudge_x
        ),
        label_y = case_when(
          ASV_ID == "Zotu17" ~ 2.2,
          ASV_ID == "Zotu29" ~ -0.55,
          ASV_ID == "Zotu13" ~ y_value,
          TRUE ~ y_value + label_nudge_y
        ),
        label_hjust = case_when(
          ASV_ID == "Zotu17" ~ 0,
          ASV_ID == "Zotu29" ~ 0.5,
          ASV_ID == "Zotu13" ~ 0,
          TRUE ~ 0
        ),
        label_vjust = case_when(
          ASV_ID == "Zotu17" ~ 0,
          ASV_ID == "Zotu29" ~ 1,
          ASV_ID == "Zotu13" ~ 0.5,
          TRUE ~ 0.5
        )
      )

    label_df <- label_df %>%
      filter(!ASV_ID %in% c("Zotu17", "Zotu29", "Zotu13"))
  }

  if (str_detect(title_text, "^D16\\s*: Prob vs Cop_pro$")) {
    manual_label_df <- bind_rows(
      manual_label_df,
      label_df %>%
        filter(ASV_ID %in% c("Zotu60", "Zotu70", "Zotu29", "Zotu17", "Zotu13")) %>%
        mutate(
          label_x = case_when(
            ASV_ID == "Zotu60" ~ x_value,
            ASV_ID == "Zotu70" ~ x_value,
            ASV_ID == "Zotu29" ~ x_value + 0.18,
            ASV_ID == "Zotu17" ~ x_value + 0.45,
            ASV_ID == "Zotu13" ~ x_value + 0.18,
            TRUE ~ x_value + label_nudge_x
          ),
          label_y = case_when(
            ASV_ID == "Zotu60" ~ -0.75,
            ASV_ID == "Zotu70" ~ 2.7,
            ASV_ID == "Zotu29" ~ 0.75,
            ASV_ID == "Zotu17" ~ 2.0,
            ASV_ID == "Zotu13" ~ -0.55,
            TRUE ~ y_value + label_nudge_y
          ),
          label_hjust = case_when(
            ASV_ID == "Zotu17" ~ 0,
            ASV_ID %in% c("Zotu60", "Zotu70", "Zotu29", "Zotu13") ~ 0.5,
            TRUE ~ 0
          ),
          label_vjust = case_when(
            ASV_ID %in% c("Zotu60", "Zotu13") ~ 1,
            ASV_ID == "Zotu29" ~ 0,
            ASV_ID == "Zotu17" ~ 0,
            TRUE ~ 0.5
          )
        )
    )

    label_df <- label_df %>%
      filter(!ASV_ID %in% c("Zotu60", "Zotu70", "Zotu29", "Zotu17", "Zotu13"))
  }

  if (str_detect(title_text, "^D21\\s*: Prob vs Cop_pro$")) {
    manual_label_df <- bind_rows(
      manual_label_df,
      plot_df %>%
        filter(ASV_ID %in% c("Zotu29", "Zotu70", "Zotu17", "Zotu13", "Zotu6")) %>%
        mutate(
          label_x = case_when(
            ASV_ID == "Zotu29" ~ x_value,
            ASV_ID == "Zotu70" ~ x_value + 0.18,
            ASV_ID == "Zotu17" ~ x_value + 0.55,
            ASV_ID == "Zotu13" ~ x_value + 0.28,
            ASV_ID == "Zotu6" ~ x_value + 0.22,
            TRUE ~ x_value + label_nudge_x
          ),
          label_y = case_when(
            ASV_ID == "Zotu29" ~ -0.35,
            ASV_ID == "Zotu70" ~ 1.45,
            ASV_ID == "Zotu17" ~ y_value,
            ASV_ID == "Zotu13" ~ 1.45,
            ASV_ID == "Zotu6" ~ y_value + 0.45,
            TRUE ~ y_value + label_nudge_y
          ),
          label_hjust = case_when(
            ASV_ID == "Zotu29" ~ 0.5,
            ASV_ID == "Zotu70" ~ 0.5,
            ASV_ID == "Zotu17" ~ 0,
            ASV_ID == "Zotu13" ~ 0.5,
            ASV_ID == "Zotu6" ~ 0,
            TRUE ~ 0
          ),
          label_vjust = case_when(
            ASV_ID == "Zotu29" ~ 1,
            ASV_ID == "Zotu70" ~ 0,
            ASV_ID == "Zotu17" ~ 0.5,
            ASV_ID == "Zotu13" ~ 0,
            ASV_ID == "Zotu6" ~ 0.5,
            TRUE ~ 0.5
          )
        )
    )

    label_df <- label_df %>%
      filter(!ASV_ID %in% c("Zotu29", "Zotu70", "Zotu17", "Zotu13", "Zotu6"))
  }

  if (str_detect(title_text, "^D21\\s*: Prob vs Cop_ctr$")) {
    manual_label_df <- bind_rows(
      manual_label_df,
      plot_df %>%
        filter(ASV_ID %in% c("Zotu29", "Zotu70", "Zotu17", "Zotu13")) %>%
        mutate(
          label_x = case_when(
            ASV_ID == "Zotu29" ~ x_value,
            ASV_ID == "Zotu70" ~ x_value + 0.18,
            ASV_ID == "Zotu17" ~ x_value + 0.55,
            ASV_ID == "Zotu13" ~ x_value + 0.28,
            TRUE ~ x_value + label_nudge_x
          ),
          label_y = case_when(
            ASV_ID == "Zotu29" ~ -0.35,
            ASV_ID == "Zotu70" ~ 1.45,
            ASV_ID == "Zotu17" ~ y_value,
            ASV_ID == "Zotu13" ~ 1.45,
            TRUE ~ y_value + label_nudge_y
          ),
          label_hjust = case_when(
            ASV_ID == "Zotu29" ~ 0.5,
            ASV_ID == "Zotu70" ~ 0.5,
            ASV_ID == "Zotu17" ~ 0,
            ASV_ID == "Zotu13" ~ 0.5,
            TRUE ~ 0
          ),
          label_vjust = case_when(
            ASV_ID == "Zotu29" ~ 1,
            ASV_ID == "Zotu70" ~ 0,
            ASV_ID == "Zotu17" ~ 0.5,
            ASV_ID == "Zotu13" ~ 0,
            TRUE ~ 0.5
          )
        )
    )

    label_df <- label_df %>%
      filter(!ASV_ID %in% c("Zotu29", "Zotu70", "Zotu17", "Zotu13"))
  }

  if (str_detect(title_text, "^D24\\s*: Prob vs Cop_pro$")) {
    manual_label_df <- bind_rows(
      manual_label_df,
      label_df %>%
        filter(ASV_ID == "Zotu29") %>%
        mutate(
          label_x = x_value + 0.55,
          label_y = 0.9,
          label_hjust = 0,
          label_vjust = 0.5
        )
    )

    label_df <- label_df %>%
      filter(ASV_ID != "Zotu29")
  }

  ggplot(plot_df, aes(x = x_value, y = y_value)) +
    geom_point(
      shape = 21,
      size = 4.2,
      stroke = 0.6,
      fill = y_color,
      color = y_color,
      alpha = 0.85
    ) +
    geom_segment(
      data = manual_label_df,
      aes(x = x_value, y = y_value, xend = label_x, yend = label_y),
      inherit.aes = FALSE,
      color = "grey35",
      linewidth = 0.35
    ) +
    geom_text(
      data = manual_label_df,
      aes(x = label_x, y = label_y, label = plot_label),
      inherit.aes = FALSE,
      family = base_family,
      size = 6.2,
      hjust = manual_label_df$label_hjust,
      vjust = manual_label_df$label_vjust
    ) +
    geom_text_repel(
      data = label_df,
      aes(label = plot_label),
      size = 6.2,
      family = base_family,
      max.overlaps = Inf,
      seed = 42,
      nudge_x = label_df$label_nudge_x,
      nudge_y = label_df$label_nudge_y,
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
    coord_cartesian(xlim = c(0, x_max_value), ylim = c(y_min_value, y_max_value), expand = TRUE) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks) +
    labs(
      title = NULL,
      x = x_label,
      y = y_label
    ) +
    theme_classic(base_size = 16, base_family = base_family) +
    theme(
      axis.title = element_text(color = "black", size = 22),
      axis.text = element_text(color = "black", size = 22)
    )
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
  mutate(
    taxonomy = coalesce(as.character(taxonomy), taxonomy_threshold)
  ) %>%
  select(ASV_ID, taxonomy, all_of(sample_cols))

analysis_tbl <- analysis_tbl %>%
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

detected_taxonomy_source <- if (!is.na(taxonomy_col)) {
  taxonomy_col
} else if (nrow(threshold_taxonomy) > 0) {
  "ASVtable_0.8threshold.csv::taxonomy_threshold"
} else {
  "None detected"
}

cat("Detected ASV column:", asv_col, "\n")
cat("Detected taxonomy column:", detected_taxonomy_source, "\n")
cat("Detected sample count columns:", paste(sample_cols, collapse = ", "), "\n\n")

for (day in days_to_plot) {
  prob_cols <- sample_cols[str_detect(sample_cols, paste0("^Prob_", day, "_"))]
  cop_pro_cols <- sample_cols[str_detect(sample_cols, paste0("^Cop_pro_", day, "_"))]
  cop_ctr_cols <- sample_cols[str_detect(sample_cols, paste0("^Cop_ctr_", day, "_"))]

  cat("Selected sample columns for", day, "\n")
  cat("  Prob:", if (length(prob_cols) > 0) paste(prob_cols, collapse = ", ") else "None", "\n")
  cat("  Cop_pro:", if (length(cop_pro_cols) > 0) paste(cop_pro_cols, collapse = ", ") else "None", "\n")
  cat("  Cop_ctr:", if (length(cop_ctr_cols) > 0) paste(cop_ctr_cols, collapse = ", ") else "None", "\n\n")

  if (length(prob_cols) == 0 || length(cop_pro_cols) == 0 || length(cop_ctr_cols) == 0) {
    warning("Skipping ", day, " because one or more required sample groups were not detected.")
    next
  }

  summary_tbl <- make_day_summary(
    df = analysis_tbl,
    day = day,
    prob_cols = prob_cols,
    cop_pro_cols = cop_pro_cols,
    cop_ctr_cols = cop_ctr_cols
  )

  day_y_axis_limit <- max(
    c(summary_tbl$mean_Cop_pro, summary_tbl$mean_Cop_ctr),
    na.rm = TRUE
  )

  if (!is.finite(day_y_axis_limit) || day_y_axis_limit <= 0) {
    day_y_axis_limit <- 1
  }

  day_y_axis_limit <- make_pretty_axis_limit(day_y_axis_limit)

  cat("First rows of", day, "summary table:\n")
  print(head(summary_tbl, 6))
  cat("\n")

  plot_prob_vs_cop_pro <- make_comparison_plot(
    summary_tbl = summary_tbl,
    x_col = "mean_Prob",
    y_col = "mean_Cop_pro",
    x_label = "Mean relative abundance (%) in Prob",
    y_label = "Mean relative abundance (%) in Cop_pro",
    title_text = paste(day, ": Prob vs Cop_pro"),
    x_color = sample_type_colors[[paste0("Prob_", day)]],
    y_color = sample_type_colors[[paste0("Cop_pro_", day)]],
    y_axis_limit = day_y_axis_limit
  )

  plot_prob_vs_cop_ctr <- make_comparison_plot(
    summary_tbl = summary_tbl,
    x_col = "mean_Prob",
    y_col = "mean_Cop_ctr",
    x_label = "Mean relative abundance (%) in Prob",
    y_label = "Mean relative abundance (%) in Cop_ctr",
    title_text = paste(day, ": Prob vs Cop_ctr"),
    x_color = sample_type_colors[[paste0("Prob_", day)]],
    y_color = sample_type_colors[[paste0("Cop_ctr_", day)]],
    y_axis_limit = day_y_axis_limit
  )

  combined_plot <- plot_prob_vs_cop_pro + plot_prob_vs_cop_ctr + plot_layout(ncol = 2)

  output_png <- file.path(results_dir, paste0("asv_relative_abundance_scatterplots_", day, ".png"))
  output_csv <- file.path(results_dir, paste0("asv_relative_abundance_summary_", day, ".csv"))

  ggsave(
    filename = output_png,
    plot = combined_plot,
    width = 16,
    height = 7,
    dpi = 300,
    bg = "white"
  )

  readr::write_csv(summary_tbl, output_csv)

  cat("Saved plot:", output_png, "\n")
  cat("Saved summary table:", output_csv, "\n\n")
}
