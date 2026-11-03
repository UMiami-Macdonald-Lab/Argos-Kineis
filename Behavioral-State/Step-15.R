# ============================================================
# CHECK ALL CLUSTERS AND THEIR OBSERVATION-WEIGHTED CONTRIBUTIONS
# ============================================================

library(dplyr)
library(tidyr)

# 1. Show each shark's composition across ALL seven clusters.
# Values are proportions: multiply by 100 for percentages.
print(
  individual_budgets,
  n = Inf,
  width = Inf
)

# 2. Calculate pooled contributions, weighting each analyzed
# daily movement row equally.
#
# This is a complementary summary, not a replacement for
# the tutorial's segment-weighted calculation.

observation_contributions <- movement_with_clusters %>%
  summarise(
    across(
      all_of(cluster_names),
      ~ mean(.x)
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to = "cluster",
    values_to = "mean_observation_prop"
  ) %>%
  arrange(desc(mean_observation_prop)) %>%
  mutate(
    cumulative_prop = cumsum(mean_observation_prop),
    
    retained_over_90_rule =
      lag(cumulative_prop, default = 0) <= 0.90
  )

print(
  observation_contributions,
  n = Inf,
  width = Inf
)