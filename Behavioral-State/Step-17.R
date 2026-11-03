# ============================================================
# PLOT CLUSTER COMPOSITION THROUGH TIME
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# Keep clusters 1–3 individually; combine 4–7 for display only.
# No proportions are discarded or renormalized.
composition_data <- movement_with_clusters

other_columns <- paste0(
  "cluster_",
  4:max_clusters
)

composition_data$other <- rowSums(
  as.matrix(
    composition_data[, other_columns, drop = FALSE]
  )
)

composition_long <- composition_data %>%
  select(
    track_id,
    system,
    date,
    cluster_1,
    cluster_2,
    cluster_3,
    other
  ) %>%
  pivot_longer(
    cols = c(
      cluster_1,
      cluster_2,
      cluster_3,
      other
    ),
    names_to = "cluster",
    values_to = "proportion"
  ) %>%
  mutate(
    cluster = factor(
      cluster,
      levels = c(
        "cluster_1",
        "cluster_2",
        "cluster_3",
        "other"
      ),
      labels = c(
        "Cluster 1",
        "Cluster 2",
        "Cluster 3",
        "Clusters 4–7"
      )
    ),
    
    system_label = factor(
      system,
      levels = c("Argos", "Kineis"),
      labels = c("Argos", "Kin\u00e9is")
    )
  )

composition_plot <- ggplot(
  composition_long,
  aes(
    x = date,
    y = proportion,
    fill = cluster
  )
) +
  geom_col(
    width = 86400 * 0.95
  ) +
  facet_wrap(
    vars(track_id, system_label),
    ncol = 2,
    scales = "free_x"
  ) +
  scale_fill_manual(
    values = c(
      "Cluster 1" = "#3B8C95",
      "Cluster 2" = "#E3AA42",
      "Cluster 3" = "#8664A5",
      "Clusters 4–7" = "#BDBDBD"
    )
  ) +
  scale_x_datetime(
    date_breaks = "2 weeks",
    date_labels = "%b %d"
  ) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    labels = c("0%", "50%", "100%")
  ) +
  labs(
    x = "Date (UTC)",
    y = "Estimated segment composition",
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

print(composition_plot)