# Percent composition for each shark and system.
three_cluster_summary <- individual_budgets %>%
  transmute(
    track_id,
    system,
    n_movement_rows,
    
    cluster_1_percent = 100 * cluster_1,
    cluster_2_percent = 100 * cluster_2,
    cluster_3_percent = 100 * cluster_3,
    
    remaining_percent =
      100 * (1 - cluster_1 - cluster_2 - cluster_3)
  )

print(
  three_cluster_summary,
  n = Inf,
  width = Inf
)