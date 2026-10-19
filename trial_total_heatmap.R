suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

work_dir <- "C:/Users/15280/OneDrive/文档/New project/outputs/Field_trapping_analysis_code_current_input/results"
input_file <- "C:/Users/15280/OneDrive/Desktop/Data availability/Field trapping reveals lower attraction of hybrid-type pheromone signals than parental types/Field attraction data.xlsx"
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
captures <- read_excel(input_file, sheet = "Captures")
trials <- read_excel(input_file, sheet = "Trials")

lure_levels <- c("A", "E", "AE", "EA")

trial_key <- trials %>%
  transmute(
    trial_id = as.character(trial_id),
    trial_start_month = as.character(trial_start_month),
    sampling_site_trial = as.character(city_or_county),
    crop = as.character(crop)
  )

dat <- captures %>%
  mutate(trial_id = as.character(trial_id)) %>%
  filter(record_status == "observed", lure_group %in% lure_levels) %>%
  left_join(trial_key, by = "trial_id") %>%
  mutate(
    sampling_site = coalesce(sampling_site_trial, as.character(sampling_site)),
    lure_group = factor(lure_group, levels = lure_levels)
  )

stopifnot(!any(is.na(dat$trial_start_month)))
stopifnot(!any(is.na(dat$sampling_site)))

trial_labels <- dat %>%
  distinct(trial_id, trial_start_month, sampling_site, crop) %>%
  arrange(trial_start_month, sampling_site, trial_id) %>%
  group_by(trial_start_month, sampling_site) %>%
  mutate(
    duplicate_site_month = n() > 1L,
    crop_clean = if_else(is.na(crop) | trimws(crop) == "", NA_character_, trimws(crop)),
    row_label = if_else(
      duplicate_site_month & !is.na(crop_clean),
      paste0(trial_start_month, " | ", sampling_site, " - ", crop_clean, " (", trial_id, ")"),
      paste0(trial_start_month, " | ", sampling_site, " (", trial_id, ")")
    )
  ) %>%
  ungroup()

trial_totals <- dat %>%
  group_by(trial_id, lure_group) %>%
  summarise(
    total_catch = sum(males_caught, na.rm = TRUE),
    n_records = n(),
    .groups = "drop"
  ) %>%
  complete(
    trial_id = unique(trial_labels$trial_id),
    lure_group = factor(lure_levels, levels = lure_levels),
    fill = list(n_records = 0L)
  ) %>%
  left_join(trial_labels, by = "trial_id")

row_order <- trial_labels %>%
  arrange(trial_start_month, sampling_site, trial_id) %>%
  pull(row_label)

max_total <- max(trial_totals$total_catch, na.rm = TRUE)
min_total <- min(trial_totals$total_catch, na.rm = TRUE)
range_log <- log1p(max_total) - log1p(min_total)

trial_totals <- trial_totals %>%
  mutate(
    row_label = factor(row_label, levels = rev(row_order)),
    cell_label = if_else(is.na(total_catch), "NA", comma(total_catch, accuracy = 1)),
    scaled_log = if_else(
      is.na(total_catch) | range_log == 0,
      0.5,
      (log1p(total_catch) - log1p(min_total)) / range_log
    ),
    label_color = if_else(!is.na(total_catch) & (scaled_log < 0.17 | scaled_log > 0.72),
                          "white", "#17212B")
  )

palette <- c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7",
             "#FDDBC7", "#EF8A62", "#B2182B")

p <- ggplot(trial_totals, aes(lure_group, row_label, fill = total_catch)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = cell_label, color = label_color),
            family = "Times New Roman", fontface = "bold", size = 2.8) +
  scale_color_identity() +
  scale_fill_gradientn(
    colours = palette,
    trans = "log1p",
    limits = c(min_total, max_total),
    breaks = c(0, 1, 10, 100, 500, 1500),
    labels = comma,
    na.value = "#D9D9D9",
    name = "Total male catches",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0,
      barheight = grid::unit(44, "mm"),
      barwidth = grid::unit(4.5, "mm")
    )
  ) +
  scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
  scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
  labs(
    title = "Male catches by sampling time, site and lure group",
    subtitle = "Each row represents one field trial; cells show the observed total catch",
    x = "Lure group",
    y = NULL,
    caption = paste0(
      "Fill colors use a log1p scale; exact totals are printed in each cell. Blank controls were excluded.\n",
      "Missing, lost-trap and installation records were excluded; ",
      "observed zeros were retained."
    )
  ) +
  theme_minimal(base_family = "Times New Roman", base_size = 10.5) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 9.5, hjust = 0, margin = margin(b = 7)),
    axis.title.x = element_text(face = "bold", size = 10.5, margin = margin(t = 7)),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 8.2, color = "black", margin = margin(r = 4)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8.5),
    plot.caption = element_text(size = 8.0, hjust = 0, margin = margin(t = 8)),
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(file.path(work_dir, "Figure_Trial_Time_Site_Lure_Total_Heatmap.tiff"), p,
       width = 7.6, height = 9.0, units = "in", dpi = 1000,
       compression = "lzw", bg = "white", limitsize = FALSE)
ggsave(file.path(work_dir, "Figure_Trial_Time_Site_Lure_Total_Heatmap.png"), p,
       width = 7.6, height = 9.0, units = "in", dpi = 300,
       bg = "white", limitsize = FALSE)
ggsave(file.path(work_dir, "Figure_Trial_Time_Site_Lure_Total_Heatmap.pdf"), p,
       width = 7.6, height = 9.0, units = "in", device = cairo_pdf,
       bg = "white", limitsize = FALSE)

export_table <- trial_totals %>%
  arrange(trial_start_month, sampling_site, trial_id, lure_group) %>%
  transmute(
    trial_id,
    trial_start_month,
    sampling_site,
    crop,
    lure_group = as.character(lure_group),
    total_males_caught = total_catch,
    observed_records = n_records
  )

write.csv(export_table,
          file.path(work_dir, "Trial_Time_Site_Lure_Totals.csv"),
          row.names = FALSE, na = "")

caption <- paste0(
  "Figure X. Total field catches of male Mythimna separata by sampling time, site and lure group. ",
  "Each row represents one field trial and is labelled with its start month, sampling site and trial identifier; ",
  "the crop is also shown when trials shared the same site and start month. ",
  "Cells report the exact total number of males captured with A, E, AE and EA lures. ",
  "Fill colors range from blue (lower catches) to red (higher catches) and use a log1p transformation ",
  "to display the wide range of totals, while the printed values remain untransformed counts. ",
  "Blank controls were excluded. Missing observations, lost traps and installation records were excluded, ",
  "whereas observed zeros were retained."
)
print(p)
cat("Trials:", n_distinct(export_table$trial_id), "\n")
cat("Heatmap cells:", nrow(export_table), "\n")
cat("Exact non-control total:", sum(export_table$total_males_caught, na.rm = TRUE), "\n")
cat("Minimum/maximum cell totals:", min_total, "/", max_total, "\n")
cat("Missing cells:", sum(is.na(export_table$total_males_caught)), "\n")
# 在 RStudio 的 Plots 窗口直接展示

