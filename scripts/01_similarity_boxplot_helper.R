# Similarity boxplot helper

library(dplyr)
library(readr)
library(stringr)
library(tibble)

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

  candidate_scores <- vapply(candidate_delimiters, function(delim) {
    field_counts <- lengths(strsplit(lines, delim, fixed = TRUE))
    valid_counts <- field_counts[field_counts > 1]

    if (length(valid_counts) == 0) {
      return(-Inf)
    }

    mode_count <- max(tabulate(match(valid_counts, unique(valid_counts))))
    mode_count * 1000 + median(valid_counts)
  }, numeric(1))

  candidate_delimiters[which.max(candidate_scores)]
}

read_similarity_matrix <- function(path = "BrayCurtisSimilarities_copepods.csv") {
  raw_lines <- readr::read_lines(path, lazy = FALSE)
  delimiter <- detect_delimiter(raw_lines)

  similarity_tbl <- readr::read_delim(
    file = path,
    delim = delimiter,
    show_col_types = FALSE,
    na = c("", "NA", "NaN"),
    trim_ws = TRUE,
    name_repair = "minimal"
  )

  names(similarity_tbl) <- clean_bom(names(similarity_tbl))
  names(similarity_tbl)[1] <- "sample_id"
  similarity_tbl$sample_id <- clean_bom(as.character(similarity_tbl$sample_id))

  similarity_matrix <- similarity_tbl %>%
    mutate(across(-sample_id, as.numeric)) %>%
    column_to_rownames("sample_id") %>%
    as.matrix()

  colnames(similarity_matrix) <- clean_bom(colnames(similarity_matrix))

  if (!setequal(rownames(similarity_matrix), colnames(similarity_matrix))) {
    stop("Row and column sample IDs do not match in ", path, ".")
  }

  similarity_matrix[rownames(similarity_matrix), rownames(similarity_matrix), drop = FALSE]
}

build_similarity_sample_meta <- function(similarity_matrix) {
  tibble(sample_id = rownames(similarity_matrix)) %>%
    mutate(
      treatment = case_when(
        str_detect(sample_id, "^Cop_pro_") ~ "Cop_pro",
        str_detect(sample_id, "^Cop_ctr_") ~ "Cop_ctr",
        TRUE ~ NA_character_
      ),
      day = str_match(sample_id, "^Cop_(?:pro|ctr)_(D(?:16|21|24))_")[, 2]
    ) %>%
    filter(!is.na(treatment), !is.na(day)) %>%
    mutate(
      treatment = factor(treatment, levels = c("Cop_pro", "Cop_ctr")),
      day = factor(day, levels = c("D16", "D21", "D24"))
    )
}

compute_mean_similarity_values <- function(similarity_matrix, focal_samples, target_samples, exclude_self = FALSE) {
  if (length(focal_samples) == 0 || length(target_samples) == 0) {
    return(tibble())
  }

  mean_rows <- lapply(focal_samples, function(sample_id) {
    current_targets <- target_samples

    if (exclude_self) {
      current_targets <- setdiff(current_targets, sample_id)
    }

    if (length(current_targets) == 0) {
      return(NULL)
    }

    tibble(
      sample_id = sample_id,
      n_compared = length(current_targets),
      similarity = mean(similarity_matrix[sample_id, current_targets], na.rm = TRUE)
    )
  })

  bind_rows(mean_rows)
}

build_day_mean_similarity_data <- function(day_label, similarity_matrix = NULL, sample_meta = NULL) {
  if (is.null(similarity_matrix)) {
    similarity_matrix <- read_similarity_matrix()
  }

  if (is.null(sample_meta)) {
    sample_meta <- build_similarity_sample_meta(similarity_matrix)
  }

  pro_samples <- sample_meta %>%
    filter(day == day_label, treatment == "Cop_pro") %>%
    pull(sample_id)

  ctr_samples <- sample_meta %>%
    filter(day == day_label, treatment == "Cop_ctr") %>%
    pull(sample_id)

  bind_rows(
    compute_mean_similarity_values(
      similarity_matrix = similarity_matrix,
      focal_samples = pro_samples,
      target_samples = pro_samples,
      exclude_self = TRUE
    ) %>%
      mutate(group = "Within_pro"),
    compute_mean_similarity_values(
      similarity_matrix = similarity_matrix,
      focal_samples = ctr_samples,
      target_samples = ctr_samples,
      exclude_self = TRUE
    ) %>%
      mutate(group = "Within_ctr"),
    bind_rows(
      compute_mean_similarity_values(
        similarity_matrix = similarity_matrix,
        focal_samples = pro_samples,
        target_samples = ctr_samples,
        exclude_self = FALSE
      ),
      compute_mean_similarity_values(
        similarity_matrix = similarity_matrix,
        focal_samples = ctr_samples,
        target_samples = pro_samples,
        exclude_self = FALSE
      )
    ) %>%
      mutate(group = "Between")
  ) %>%
    mutate(day = day_label)
}

make_similarity_boxplot <- function(day_label, output_name, box_colors) {
  results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

  data <- build_day_mean_similarity_data(day_label)

  label_map <- c(
    "within_pro" = "Cop_pro",
    "within_ctr" = "Cop_ctr",
    "between" = "Cop_pro vs. Cop_ctr"
  )

  data$group_key <- tolower(data$group)
  data$group_label <- unname(label_map[data$group_key])

  if (any(is.na(data$group_label))) {
    missing_levels <- unique(data$group[is.na(data$group_label)])
    stop(sprintf("Unmapped group labels found: %s", paste(missing_levels, collapse = ", ")))
  }

  data$group_label <- factor(
    data$group_label,
    levels = c(
      "Cop_pro",
      "Cop_ctr",
      "Cop_pro vs. Cop_ctr"
    )
  )

  output_file <- file.path(results_dir, output_name)
  png(filename = output_file, width = 1000, height = 700, res = 120)

  axis_text_cex <- 1.5
  x_axis_text_cex <- 1.45
  axis_title_cex <- 1.8

  par(family = "serif", bty = "l", cex.axis = axis_text_cex, cex.lab = axis_title_cex)
  par(mar = c(5.5, 6.5, 2, 1))
  par(yaxs = "i")

  boxplot(
    similarity ~ group_label,
    data = data,
    at = c(1, 2.4, 3.8),
    xaxt = "n",
    yaxt = "n",
    frame.plot = FALSE,
    xlab = "",
    ylab = "",
    ylim = c(0, 1),
    col = box_colors
  )
  axis(1, at = c(1, 2.4, 3.8), labels = levels(data$group_label), cex.axis = x_axis_text_cex)
  axis(2, at = seq(0, 1, by = 0.2), cex.axis = axis_text_cex, las = 1)
  title(ylab = "Mean Bray-Curtis similarity", cex.lab = axis_title_cex)
  box(bty = "l")

  dev.off()
  cat("Saved boxplot to:", output_file, "\n")
}
